"""Standalone pileup-to-reference ranker activated only after reranker failure."""

from __future__ import annotations

import itertools

import torch
from torch import nn

from . import MODALITIES
from .model import AttentiveGroupPool, ReferenceEmbeddingStore


class PileupEvidenceEncoder(nn.Module):
    """Encode pileup rows with a separate input adapter for each modality."""

    def __init__(self, input_dim=44, hidden_dim=128, layers=4, heads=4, dropout=0.1,
                 max_positions=8192):
        super().__init__()
        self.input_dim = int(input_dim)
        self.hidden_dim = int(hidden_dim)
        self.layers = int(layers)
        self.heads = int(heads)
        self.dropout = float(dropout)
        self.max_positions = int(max_positions)
        self.adapters = nn.ModuleDict({
            modality: nn.Sequential(
                nn.LayerNorm(input_dim),
                nn.Linear(input_dim, hidden_dim),
                nn.GELU(),
            )
            for modality in MODALITIES
        })
        layer = nn.TransformerEncoderLayer(
            hidden_dim, heads, hidden_dim * 4, dropout=dropout,
            batch_first=True, norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(layer, layers)
        self.modality = nn.Embedding(len(MODALITIES), hidden_dim)
        self.position = nn.Embedding(max_positions, hidden_dim)
        self.norm = nn.LayerNorm(hidden_dim)

    def forward(self, pileup: torch.Tensor, mask: torch.Tensor, modality: str):
        if modality not in MODALITIES:
            raise ValueError(f"unsupported modality: {modality}")
        if pileup.ndim != 3 or pileup.shape[-1] != self.input_dim:
            raise ValueError(
                f"pileup must have shape [batch, positions, {self.input_dim}], got {tuple(pileup.shape)}"
            )
        if mask.shape != pileup.shape[:2]:
            raise ValueError(f"mask shape {tuple(mask.shape)} does not match pileup {tuple(pileup.shape[:2])}")
        if pileup.shape[1] > self.max_positions:
            raise ValueError(f"pileup length {pileup.shape[1]} exceeds {self.max_positions} positions")
        modality_index = torch.tensor(MODALITIES.index(modality), device=pileup.device)
        positions = torch.arange(pileup.shape[1], device=pileup.device)
        hidden = self.adapters[modality](pileup) + self.modality(modality_index) + self.position(positions)
        encoded = self.encoder(hidden, src_key_padding_mask=~mask.bool())
        weights = mask.float().unsqueeze(-1)
        return self.norm((encoded * weights).sum(1) / weights.sum(1).clamp_min(1.0))

    def config(self) -> dict:
        return {
            "input_dim": self.input_dim,
            "hidden_dim": self.hidden_dim,
            "layers": self.layers,
            "heads": self.heads,
            "dropout": self.dropout,
            "max_positions": self.max_positions,
            "modalities": list(MODALITIES),
            "modality_specific_adapters": True,
        }


class StandaloneReferenceRanker(nn.Module):
    """Two-stage allele retrieval followed by symmetric genotype reranking."""

    def __init__(self, reference_dim=128, hidden_dim=128, evidence_layers=4,
                 evidence_heads=4, dropout=0.1, max_positions=8192):
        super().__init__()
        self.reference_dim = int(reference_dim)
        self.hidden_dim = int(hidden_dim)
        self.evidence = PileupEvidenceEncoder(
            hidden_dim=hidden_dim, layers=evidence_layers,
            heads=evidence_heads, dropout=dropout, max_positions=max_positions,
        )
        self.group_pool = AttentiveGroupPool(reference_dim)
        self.allele_project = nn.Linear(reference_dim, hidden_dim)
        self.pair = nn.Sequential(
            nn.Linear(hidden_dim * 4, hidden_dim), nn.GELU(),
            nn.LayerNorm(hidden_dim), nn.Linear(hidden_dim, 1),
        )
        self.log_similarity_scale = nn.Parameter(torch.tensor(0.0))

    def allele_vector(self, allele: str, store: ReferenceEmbeddingStore,
                      device: torch.device) -> torch.Tensor:
        values, resolved = store.lookup(allele, device)
        if not resolved:
            raise KeyError(f"allele is absent from the pinned embedding cache: {allele}")
        return self.allele_project(self.group_pool(values))

    @staticmethod
    def symmetric_pair_features(evidence: torch.Tensor, first: torch.Tensor,
                                second: torch.Tensor) -> torch.Tensor:
        return torch.cat((evidence, first + second, torch.abs(first - second), first * second), dim=-1)

    def score_pairs(self, evidence: torch.Tensor, pairs: list[tuple[str, str]],
                    store: ReferenceEmbeddingStore) -> torch.Tensor:
        """Score unique canonical pairs for one encoded locus."""
        if evidence.ndim != 1:
            raise ValueError(f"evidence must be a vector, got {tuple(evidence.shape)}")
        canonical = [tuple(sorted(pair)) for pair in pairs]
        if len(set(canonical)) != len(canonical):
            raise ValueError("candidate list contains duplicate unordered allele pairs")
        vectors, similarities = [], []
        for first_name, second_name in canonical:
            first = self.allele_vector(first_name, store, evidence.device)
            second = self.allele_vector(second_name, store, evidence.device)
            vectors.append(self.symmetric_pair_features(evidence, first, second))
            pair_mean = 0.5 * (first + second)
            similarities.append(nn.functional.cosine_similarity(evidence, pair_mean, dim=0))
        if not vectors:
            raise ValueError("at least one candidate pair is required")
        learned = self.pair(torch.stack(vectors)).squeeze(-1)
        scale = self.log_similarity_scale.clamp(-2.0, 4.0).exp()
        return learned + scale * torch.stack(similarities)

    def forward_candidates(self, pileup: torch.Tensor, mask: torch.Tensor, modality: str,
                           pairs: list[tuple[str, str]], store: ReferenceEmbeddingStore) -> torch.Tensor:
        evidence = self.evidence(pileup, mask, modality)
        if evidence.shape[0] != 1:
            raise ValueError("forward_candidates currently accepts one locus at a time")
        return self.score_pairs(evidence.squeeze(0), pairs, store)

    def retrieve(self, pileup, mask, modality: str, gene: str, store: ReferenceEmbeddingStore,
                 top_k_alleles: int = 20):
        device = next(self.parameters()).device
        evidence = self.evidence(pileup.to(device), mask.to(device), modality).squeeze(0)
        allele_vectors, allele_names = [], []
        for allele in sorted(store.embeddings):
            if allele.split("*", 1)[0] != gene:
                continue
            allele_names.append(allele)
            allele_vectors.append(self.allele_vector(allele, store, device))
        if not allele_names:
            raise ValueError(f"embedding cache contains no alleles for HLA-{gene}")
        matrix = torch.stack(allele_vectors)
        single_scores = matrix @ evidence
        selected = torch.topk(single_scores, min(top_k_alleles, len(allele_names))).indices.tolist()
        pairs = [
            tuple(sorted((allele_names[first], allele_names[second])))
            for first, second in itertools.combinations_with_replacement(selected, 2)
        ]
        scores = self.score_pairs(evidence, pairs, store)
        order = torch.argsort(scores, descending=True)
        return [(pairs[index], float(scores[index].detach().cpu())) for index in order.tolist()]

    def config(self) -> dict:
        return {
            "architecture": "StandaloneReferenceRanker",
            "reference_dim": self.reference_dim,
            "hidden_dim": self.hidden_dim,
            "evidence_encoder": self.evidence.config(),
            "pair_representation": ["evidence", "sum", "absolute_difference", "elementwise_product"],
            "unordered_pair_scoring": True,
            "direct_evidence_reference_similarity": True,
        }
