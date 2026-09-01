"""Compact reference encoder and permutation-invariant candidate set transformer."""

from __future__ import annotations

import json
import random
from pathlib import Path

import torch
from torch import nn

from . import GENES, MODALITIES
from .alleles import normalize_allele
from .data import FAMILIES, GLOBAL_FEATURE_NAMES, TOOLS, TOOL_FEATURE_NAMES


class NucleotideEncoder(nn.Module):
    def __init__(self, vocab_size: int, hidden_dim: int = 128, layers: int = 4, heads: int = 4,
                 max_tokens: int = 1024, dropout: float = 0.1):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.token = nn.Embedding(vocab_size, hidden_dim, padding_idx=0)
        self.position = nn.Embedding(max_tokens, hidden_dim)
        layer = nn.TransformerEncoderLayer(hidden_dim, heads, hidden_dim * 4, dropout=dropout,
                                           batch_first=True, norm_first=True)
        self.encoder = nn.TransformerEncoder(layer, layers)
        self.norm = nn.LayerNorm(hidden_dim)
        self.mlm_head = nn.Linear(hidden_dim, vocab_size)
        self.projection = nn.Linear(hidden_dim, hidden_dim)

    def forward(self, tokens: torch.Tensor, mask: torch.Tensor) -> dict[str, torch.Tensor]:
        positions = torch.arange(tokens.shape[1], device=tokens.device).unsqueeze(0)
        hidden = self.token(tokens) + self.position(positions)
        encoded = self.encoder(hidden, src_key_padding_mask=~mask.bool())
        encoded = self.norm(encoded)
        weights = mask.float().unsqueeze(-1)
        pooled = (encoded * weights).sum(1) / weights.sum(1).clamp_min(1.0)
        return {"encoded": encoded, "pooled": self.projection(pooled), "mlm_logits": self.mlm_head(encoded)}


class ReferenceEmbeddingStore:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        payload = torch.load(self.path, map_location="cpu")
        self.embeddings = payload["embeddings"]
        self.reference_sha256 = payload["reference_sha256"]
        self.encoder_sha256 = payload["encoder_sha256"]
        self.embedding_dim = int(payload["embedding_dim"])
        self.metadata = payload.get("metadata", {})
        self.cache_sha256 = __import__("hashlib").sha256(self.path.read_bytes()).hexdigest()

    def lookup(self, allele: str, device: torch.device) -> tuple[torch.Tensor, bool]:
        value = self.embeddings.get(normalize_allele(allele))
        if value is None:
            return torch.zeros((1, self.embedding_dim), dtype=torch.float32, device=device), False
        return value.to(device=device, dtype=torch.float32), True


class AttentiveGroupPool(nn.Module):
    def __init__(self, embedding_dim: int):
        super().__init__()
        self.query = nn.Parameter(torch.zeros(embedding_dim))
        nn.init.normal_(self.query, std=0.02)

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        weights = torch.softmax(values @ self.query, dim=0)
        return (values * weights.unsqueeze(-1)).sum(0)


class RefFormerModel(nn.Module):
    """Scores one list of caller-proposed candidates at a time."""

    def __init__(self, reference_dim: int = 128, hidden_dim: int = 128, layers: int = 2,
                 heads: int = 4, dropout: float = 0.1, no_imgt: bool = False):
        super().__init__()
        self.reference_dim = reference_dim
        self.hidden_dim = hidden_dim
        self.no_imgt = no_imgt
        self.dropout = dropout
        self.group_pool = AttentiveGroupPool(reference_dim)
        self.unknown_allele = nn.Parameter(torch.zeros(reference_dim))
        nn.init.normal_(self.unknown_allele, std=0.02)
        self.pair_project = nn.Sequential(nn.Linear(reference_dim * 3, hidden_dim), nn.GELU(),
                                          nn.LayerNorm(hidden_dim))
        self.global_project = nn.Linear(len(GLOBAL_FEATURE_NAMES), hidden_dim)
        self.tool_numeric = nn.Linear(len(TOOL_FEATURE_NAMES), hidden_dim)
        self.tool_embedding = nn.Embedding(len(TOOLS), hidden_dim)
        family_names = sorted(set(FAMILIES.values()))
        self.family_names = tuple(family_names)
        self.family_embedding = nn.Embedding(len(family_names), hidden_dim)
        self.modality_embedding = nn.Embedding(len(MODALITIES), hidden_dim)
        self.gene_embedding = nn.Embedding(len(GENES), hidden_dim)
        self.candidate_type = nn.Parameter(torch.zeros(hidden_dim))
        self.tool_type = nn.Parameter(torch.zeros(hidden_dim))
        layer = nn.TransformerEncoderLayer(hidden_dim, heads, hidden_dim * 4, dropout=dropout,
                                           batch_first=True, norm_first=True)
        self.set_transformer = nn.TransformerEncoder(layer, layers)
        self.score_head = nn.Sequential(nn.LayerNorm(hidden_dim), nn.Linear(hidden_dim, 1))
        self.oracle_head = nn.Sequential(nn.LayerNorm(hidden_dim), nn.Linear(hidden_dim, 1))

    def _allele(self, allele: str, store: ReferenceEmbeddingStore, device: torch.device) -> torch.Tensor:
        if self.no_imgt:
            return torch.zeros(self.reference_dim, device=device)
        values, resolved = store.lookup(allele, device)
        return self.group_pool(values) if resolved else self.unknown_allele

    def _pair(self, pair, store: ReferenceEmbeddingStore, device: torch.device) -> torch.Tensor:
        if not pair or tuple(pair) == ("", ""):
            return torch.zeros(self.hidden_dim, device=device)
        first = self._allele(pair[0], store, device)
        second = self._allele(pair[1], store, device)
        symmetric = torch.cat((first + second, torch.abs(first - second), first * second), dim=-1)
        return self.pair_project(symmetric)

    def forward_locus(self, locus: dict, store: ReferenceEmbeddingStore,
                      tool_order: list[int] | None = None) -> tuple[torch.Tensor, torch.Tensor]:
        device = next(self.parameters()).device
        modality_idx = MODALITIES.index(locus["modality"])
        gene_idx = GENES.index(locus["gene"])
        context = self.modality_embedding(torch.tensor(modality_idx, device=device)) + \
            self.gene_embedding(torch.tensor(gene_idx, device=device))
        tokens, candidate_positions = [], []
        for candidate in locus["candidates"]:
            candidate_token = self._pair(candidate["pair"], store, device) + \
                self.global_project(torch.tensor(candidate["global_features"], dtype=torch.float32, device=device)) + \
                context + self.candidate_type
            candidate_positions.append(len(tokens))
            tokens.append(candidate_token)
            order = tool_order or list(range(len(candidate["tool_tokens"])))
            for index in order:
                proposal = candidate["tool_tokens"][index]
                tool_idx = TOOLS.index(proposal["tool"])
                family_idx = self.family_names.index(proposal["family"])
                token = self._pair(proposal["pair"], store, device) + \
                    self.tool_numeric(torch.tensor(proposal["features"], dtype=torch.float32, device=device)) + \
                    self.tool_embedding(torch.tensor(tool_idx, device=device)) + \
                    self.family_embedding(torch.tensor(family_idx, device=device)) + context + self.tool_type
                tokens.append(token)
        # A single set-transformer pass makes the candidates interact. With no
        # positional encodings, permuting caller tokens leaves scores invariant.
        encoded = self.set_transformer(torch.stack(tokens).unsqueeze(0))
        candidate_encoded = encoded[0, candidate_positions]
        scores = self.score_head(candidate_encoded).squeeze(-1)
        oracle_logit = self.oracle_head(encoded.mean(1)[0]).squeeze(-1)
        return scores, oracle_logit

    def shuffled_tool_order(self, locus: dict, rng: random.Random) -> list[int]:
        order = list(range(len(locus["candidates"][0]["tool_tokens"])))
        rng.shuffle(order)
        return order

    def config(self) -> dict:
        return {"reference_dim": self.reference_dim, "hidden_dim": self.hidden_dim,
                "layers": len(self.set_transformer.layers),
                "heads": self.set_transformer.layers[0].self_attn.num_heads,
                "dropout": self.dropout,
                "no_imgt": self.no_imgt, "tools": list(TOOLS), "families": list(self.family_names),
                "genes": list(GENES), "modalities": list(MODALITIES)}
