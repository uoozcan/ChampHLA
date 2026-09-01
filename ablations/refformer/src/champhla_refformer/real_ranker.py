"""Batched real-pileup ranker over all pinned IMGT unordered allele pairs."""

from __future__ import annotations

import itertools
import math

import torch
from torch import nn

from . import GENES, MODALITIES
from .model import ReferenceEmbeddingStore


class RealPileupPairRanker(nn.Module):
    """Encode binned pileups and jointly rank every reference-compatible pair."""

    def __init__(self, store: ReferenceEmbeddingStore, hidden_dim: int = 64,
                 layers: int = 2, heads: int = 4, bins: int = 128,
                 dropout: float = 0.1, no_imgt: bool = False):
        super().__init__()
        self.reference_dim = int(store.embedding_dim)
        self.hidden_dim = int(hidden_dim)
        self.layers = int(layers)
        self.heads = int(heads)
        self.bins = int(bins)
        self.dropout = float(dropout)
        self.no_imgt = bool(no_imgt)
        self.adapters = nn.ModuleList([
            nn.Sequential(nn.LayerNorm(88), nn.Linear(88, hidden_dim), nn.GELU())
            for _ in MODALITIES
        ])
        self.modality_embedding = nn.Embedding(len(MODALITIES), hidden_dim)
        self.gene_embedding = nn.Embedding(len(GENES), hidden_dim)
        self.position = nn.Embedding(bins, hidden_dim)
        layer = nn.TransformerEncoderLayer(
            hidden_dim, heads, hidden_dim * 4, dropout=dropout,
            batch_first=True, norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(layer, layers)
        self.pool_query = nn.Parameter(torch.empty(hidden_dim))
        nn.init.normal_(self.pool_query, std=0.02)
        self.evidence_norm = nn.LayerNorm(hidden_dim)
        self.allele_project = nn.Sequential(
            nn.LayerNorm(self.reference_dim), nn.Linear(self.reference_dim, hidden_dim), nn.GELU(),
            nn.LayerNorm(hidden_dim),
        )
        self.pair_project = nn.Sequential(
            nn.Linear(hidden_dim * 3, hidden_dim), nn.GELU(), nn.LayerNorm(hidden_dim),
        )
        self.logit_scale = nn.Parameter(torch.tensor(math.log(8.0)))

        self.allele_names: dict[str, tuple[str, ...]] = {}
        self.pair_names: dict[str, tuple[tuple[str, str], ...]] = {}
        for gene in GENES:
            names = tuple(sorted(a for a in store.embeddings if a.startswith(gene + "*")))
            if not names:
                raise ValueError(f"embedding cache contains no HLA-{gene} alleles")
            means = torch.stack([store.embeddings[a].float().mean(0) for a in names])
            means = nn.functional.normalize(means, dim=-1)
            first, second = zip(*itertools.combinations_with_replacement(range(len(names)), 2))
            self.register_buffer(f"reference_{gene}", means)
            self.register_buffer(f"pair_first_{gene}", torch.tensor(first, dtype=torch.long))
            self.register_buffer(f"pair_second_{gene}", torch.tensor(second, dtype=torch.long))
            self.allele_names[gene] = names
            self.pair_names[gene] = tuple((names[i], names[j]) for i, j in zip(first, second))
            if no_imgt:
                learned = nn.Parameter(torch.empty(len(names), self.reference_dim))
                nn.init.normal_(learned, std=0.05)
                self.register_parameter(f"learned_{gene}", learned)

    def encode(self, binned: torch.Tensor, modality: torch.Tensor, gene: str) -> torch.Tensor:
        if binned.ndim != 3 or binned.shape[1:] != (self.bins, 88):
            raise ValueError(f"binned pileup must be [batch,{self.bins},88], got {tuple(binned.shape)}")
        hidden = torch.zeros((binned.shape[0], self.bins, self.hidden_dim),
                             dtype=binned.dtype, device=binned.device)
        for modality_index, adapter in enumerate(self.adapters):
            selected = modality == modality_index
            if selected.any():
                hidden[selected] = adapter(binned[selected])
        pos = self.position(torch.arange(self.bins, device=binned.device)).unsqueeze(0)
        context = self.modality_embedding(modality).unsqueeze(1) + \
            self.gene_embedding(torch.tensor(GENES.index(gene), device=binned.device)).view(1, 1, -1)
        encoded = self.encoder(hidden + pos + context)
        weights = torch.softmax(encoded @ self.pool_query, dim=1)
        return self.evidence_norm((encoded * weights.unsqueeze(-1)).sum(1))

    def pair_vectors(self, gene: str) -> torch.Tensor:
        raw = getattr(self, f"learned_{gene}") if self.no_imgt else getattr(self, f"reference_{gene}")
        allele = self.allele_project(raw)
        first = allele[getattr(self, f"pair_first_{gene}")]
        second = allele[getattr(self, f"pair_second_{gene}")]
        symmetric = torch.cat((first + second, torch.abs(first - second), first * second), dim=-1)
        return self.pair_project(symmetric)

    def forward(self, binned: torch.Tensor, modality: torch.Tensor, gene: str) -> torch.Tensor:
        evidence = nn.functional.normalize(self.encode(binned, modality, gene), dim=-1)
        pairs = nn.functional.normalize(self.pair_vectors(gene), dim=-1)
        scale = self.logit_scale.clamp(math.log(1.0), math.log(100.0)).exp()
        return scale * (evidence @ pairs.T)

    def pair_index(self, gene: str, pair: tuple[str, str]) -> int:
        canonical = tuple(sorted(pair))
        try:
            return self.pair_names[gene].index(canonical)
        except ValueError as exc:
            raise KeyError(f"pair is absent from HLA-{gene} reference space: {canonical}") from exc

    def config(self) -> dict:
        return {
            "architecture": "RealPileupPairRanker",
            "reference_dim": self.reference_dim,
            "hidden_dim": self.hidden_dim,
            "layers": self.layers,
            "heads": self.heads,
            "bins": self.bins,
            "dropout": self.dropout,
            "no_imgt": self.no_imgt,
            "modalities": list(MODALITIES),
            "genes": list(GENES),
            "unordered_pair_scoring": True,
            "multi_sequence_pooling": "mean_over_cached_compatible_higher_field_embeddings",
        }


def bin_pileup(pileup: torch.Tensor, bins: int) -> torch.Tensor:
    """Preserve both average evidence and sharp local peaks in fixed bins."""
    if pileup.ndim != 2 or pileup.shape[1] != 44:
        raise ValueError(f"pileup must be [positions,44], got {tuple(pileup.shape)}")
    values = pileup.T.unsqueeze(0)
    average = nn.functional.adaptive_avg_pool1d(values, bins)
    maximum = nn.functional.adaptive_max_pool1d(values, bins)
    return torch.cat((average, maximum), dim=1).squeeze(0).T.contiguous()
