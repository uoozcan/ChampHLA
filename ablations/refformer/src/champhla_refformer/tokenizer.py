"""Deterministic overlapping nucleotide k-mer tokenizer."""

from __future__ import annotations

import itertools
from dataclasses import dataclass


SPECIAL = ("<PAD>", "<MASK>", "<UNK>", "<CLS>")
DNA = "ACGT"


@dataclass(frozen=True)
class KmerTokenizer:
    k: int = 6
    max_tokens: int = 1024

    def __post_init__(self):
        vocab = list(SPECIAL) + ["".join(chars) for chars in itertools.product(DNA, repeat=self.k)]
        object.__setattr__(self, "vocab", tuple(vocab))
        object.__setattr__(self, "index", {token: idx for idx, token in enumerate(vocab)})

    @property
    def pad_id(self) -> int:
        return self.index["<PAD>"]

    @property
    def mask_id(self) -> int:
        return self.index["<MASK>"]

    def encode(self, sequence: str) -> tuple[list[int], list[int]]:
        sequence = "".join(base if base in DNA else "N" for base in sequence.upper())
        tokens = [self.index["<CLS>"]]
        for index in range(max(0, len(sequence) - self.k + 1)):
            kmer = sequence[index:index + self.k]
            tokens.append(self.index.get(kmer, self.index["<UNK>"]))
            if len(tokens) >= self.max_tokens:
                break
        mask = [1] * len(tokens)
        if len(tokens) < self.max_tokens:
            padding = self.max_tokens - len(tokens)
            tokens.extend([self.pad_id] * padding)
            mask.extend([0] * padding)
        return tokens, mask

    def manifest(self) -> dict:
        return {"type": "overlapping_kmer", "k": self.k, "max_tokens": self.max_tokens,
                "vocab_size": len(self.vocab), "special_tokens": list(SPECIAL)}

