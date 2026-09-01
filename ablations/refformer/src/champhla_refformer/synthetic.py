"""Synthetic IMGT-derived pileups for label-free standalone-ranker pretraining."""

from __future__ import annotations

import hashlib
import math
import random
from dataclasses import dataclass

import torch

from . import MODALITIES


BASE_COLUMNS = {"A": 5, "C": 6, "G": 7, "T": 8, "N": 9}
REF_COLUMNS = {"A": 1, "C": 2, "G": 3, "T": 4}


@dataclass(frozen=True)
class SyntheticPileupConfig:
    max_positions: int = 384
    input_dim: int = 44
    mean_depth_wgs: int = 30
    mean_depth_wes: int = 80
    mean_depth_rnaseq: int = 50
    dropout_wgs: float = 0.02
    dropout_wes: float = 0.08
    dropout_rnaseq: float = 0.15
    base_error_rate: float = 0.005


def stable_pair_split(pair: tuple[str, str], seed: int, validation_fraction: float = 0.15) -> str:
    canonical = "|".join(sorted(pair))
    digest = hashlib.sha256(f"{seed}:{canonical}".encode()).hexdigest()
    return "validation" if int(digest[:8], 16) / 0xFFFFFFFF < validation_fraction else "train"


def _depth_and_dropout(modality: str, config: SyntheticPileupConfig) -> tuple[int, float]:
    if modality not in MODALITIES:
        raise ValueError(f"unsupported modality: {modality}")
    return (
        int(getattr(config, f"mean_depth_{modality}")),
        float(getattr(config, f"dropout_{modality}")),
    )


def sequence_pair_to_pileup(sequence1: str, sequence2: str, modality: str,
                            generator: torch.Generator,
                            config: SyntheticPileupConfig | None = None) -> tuple[torch.Tensor, torch.Tensor]:
    """Create a deterministic diploid pseudo-pileup without using subject truth."""
    config = config or SyntheticPileupConfig()
    sequence1 = sequence1.upper()
    sequence2 = sequence2.upper()
    usable = min(len(sequence1), len(sequence2))
    if usable == 0:
        raise ValueError("cannot simulate a pileup from an empty reference sequence")
    positions = min(config.max_positions, usable)
    indices = torch.linspace(0, usable - 1, positions).round().long().tolist()
    mean_depth, dropout = _depth_and_dropout(modality, config)
    rows = torch.zeros((positions, config.input_dim), dtype=torch.float32)
    mask = torch.ones(positions, dtype=torch.bool)
    for row_index, sequence_index in enumerate(indices):
        if float(torch.rand((), generator=generator)) < dropout:
            mask[row_index] = False
            continue
        depth_scale = 0.5 + float(torch.rand((), generator=generator))
        depth = max(1, int(round(mean_depth * depth_scale)))
        first = sequence1[sequence_index] if sequence1[sequence_index] in "ACGT" else None
        second = sequence2[sequence_index] if sequence2[sequence_index] in "ACGT" else None
        if first is None and second is None:
            mask[row_index] = False
            continue
        reference_base = first or second
        first_fraction = (0.45 + 0.10 * float(torch.rand((), generator=generator))) if (
            first is not None and second is not None
        ) else (1.0 if first is not None else 0.0)
        fractions = {base: 0.0 for base in BASE_COLUMNS}
        if first is not None:
            fractions[first] += first_fraction
        if second is not None:
            fractions[second] += 1.0 - first_fraction
        error = min(config.base_error_rate, 0.2)
        for base in fractions:
            fractions[base] *= 1.0 - error
        fractions["N"] += error
        rows[row_index, 0] = min(1.0, math.log1p(depth) / math.log1p(200.0))
        rows[row_index, REF_COLUMNS[reference_base]] = 1.0
        for base, column in BASE_COLUMNS.items():
            rows[row_index, column] = fractions[base]
        rows[row_index, 10] = fractions.get(reference_base, 0.0)
        rows[row_index, 14] = sum(fractions.values())
        rows[row_index, 15] = 0.58 + 0.08 * float(torch.rand((), generator=generator))
        rows[row_index, 16] = 0.9
        rows[row_index, 17] = row_index / max(1, positions - 1)
        rows[row_index, 18] = 1.0
        rows[row_index, 19] = min(1.0, depth / 100.0)
        rows[row_index, 20] = abs(first_fraction - (1.0 - first_fraction))
        rows[row_index, 21] = 1.0 if first is not None and first == second else 0.0
        rows[row_index, 22] = 1.0 if first is not None and second is not None and first != second else 0.0
    if not bool(mask.any()):
        mask[0] = True
    return rows.unsqueeze(0), mask.unsqueeze(0)


def sample_candidate_pairs(alleles: list[str], truth: tuple[str, str], count: int,
                           rng: random.Random) -> list[tuple[str, str]]:
    """Return one true pair plus unique hard and random negative diploid pairs."""
    truth = tuple(sorted(truth))
    if count < 2:
        raise ValueError("candidate count must be at least two")
    if len(alleles) < 2:
        raise ValueError("at least two alleles are required")
    negatives = []
    seen = {truth}
    for anchor in truth:
        shuffled = list(alleles)
        rng.shuffle(shuffled)
        for other in shuffled:
            candidate = tuple(sorted((anchor, other)))
            if candidate not in seen:
                negatives.append(candidate)
                seen.add(candidate)
                break
    attempts = 0
    while len(negatives) < count - 1 and attempts < count * 100:
        candidate = tuple(sorted((rng.choice(alleles), rng.choice(alleles))))
        if candidate not in seen:
            negatives.append(candidate)
            seen.add(candidate)
        attempts += 1
    rng.shuffle(negatives)
    return [truth, *negatives[:count - 1]]
