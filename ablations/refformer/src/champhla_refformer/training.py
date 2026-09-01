"""Training, calibration, prediction, and checksummed model bundles."""

from __future__ import annotations

import hashlib
import json
import math
import random
from copy import deepcopy
from pathlib import Path

import torch
from torch import nn

from . import SCHEMA_VERSION
from .model import RefFormerModel, ReferenceEmbeddingStore


def seed_everything(seed: int):
    random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def labels(locus: dict, device) -> torch.Tensor:
    return torch.tensor([float(candidate.get("label", 0)) for candidate in locus["candidates"]],
                        dtype=torch.float32, device=device)


def accuracy(model: RefFormerModel, loci: list[dict], store: ReferenceEmbeddingStore) -> float:
    if not loci:
        return 0.0
    model.eval()
    correct = 0
    with torch.no_grad():
        for locus in loci:
            scores, _ = model.forward_locus(locus, store)
            selected = int(scores.argmax().item())
            correct += int(locus["candidates"][selected].get("label", 0) == 1)
    return correct / len(loci)


def train_model(train_loci: list[dict], val_loci: list[dict], store: ReferenceEmbeddingStore,
                config: dict, seed: int, no_imgt: bool = False) -> tuple[RefFormerModel, dict]:
    seed_everything(seed)
    device = torch.device("cuda" if torch.cuda.is_available() and config.get("device", "auto") != "cpu" else "cpu")
    model = RefFormerModel(reference_dim=store.embedding_dim,
                           hidden_dim=int(config.get("hidden_dim", 128)),
                           layers=int(config.get("set_layers", 2)), heads=int(config.get("heads", 4)),
                           dropout=float(config.get("dropout", 0.1)), no_imgt=no_imgt).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=float(config.get("learning_rate", 3e-4)),
                                  weight_decay=float(config.get("weight_decay", 1e-4)))
    bce = nn.BCEWithLogitsLoss()
    rng = random.Random(seed)
    best_state, best_accuracy, history = None, -1.0, []
    patience = int(config.get("patience", 8))
    stale = 0
    for epoch in range(1, int(config.get("epochs", 40)) + 1):
        model.train()
        order = list(range(len(train_loci)))
        rng.shuffle(order)
        losses = []
        optimizer.zero_grad(set_to_none=True)
        accumulation = max(1, int(config.get("gradient_accumulation", 8)))
        for step, index in enumerate(order, 1):
            locus = train_loci[index]
            tool_order = model.shuffled_tool_order(locus, rng)
            scores, oracle_logit = model.forward_locus(locus, store, tool_order)
            target = labels(locus, device)
            oracle = torch.tensor(float(locus.get("candidate_set_oracle", 0)), device=device)
            loss = float(config.get("oracle_loss_weight", 0.2)) * bce(oracle_logit, oracle)
            if target.sum() > 0:
                truth_index = target.argmax().long().unsqueeze(0)
                loss = loss + nn.functional.cross_entropy(scores.unsqueeze(0), truth_index)
            (loss / accumulation).backward()
            losses.append(float(loss.detach().cpu()))
            if step % accumulation == 0 or step == len(order):
                nn.utils.clip_grad_norm_(model.parameters(), float(config.get("clip_grad_norm", 1.0)))
                optimizer.step()
                optimizer.zero_grad(set_to_none=True)
        val_accuracy = accuracy(model, val_loci or train_loci, store)
        record = {"epoch": epoch, "train_loss": sum(losses) / max(1, len(losses)),
                  "validation_top_call_accuracy": val_accuracy}
        history.append(record)
        if val_accuracy > best_accuracy + 1e-12:
            best_accuracy, best_state, stale = val_accuracy, deepcopy(model.state_dict()), 0
        else:
            stale += 1
        if stale >= patience:
            break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, {"best_validation_accuracy": best_accuracy, "epochs_completed": len(history),
                   "history": history, "device": str(device), "seed": seed, "no_imgt": no_imgt}


def fit_temperature(raw: list[tuple[list[float], int]]) -> float:
    if not raw:
        return 1.0
    best = (float("inf"), 1.0)
    for step in range(20, 401):
        temperature = step / 100.0
        loss = 0.0
        for scores, truth_index in raw:
            values = torch.tensor(scores, dtype=torch.float64) / temperature
            loss -= float(torch.log_softmax(values, 0)[truth_index])
        best = min(best, (loss / len(raw), temperature))
    return best[1]


def calibrated_probabilities(scores: torch.Tensor, temperature: float) -> torch.Tensor:
    return torch.softmax(scores / max(1e-6, temperature), dim=0)


def predict_loci(model: RefFormerModel, loci: list[dict], store: ReferenceEmbeddingStore,
                 temperature: float = 1.0, threshold: float | None = None,
                 model_sha256: str = "") -> tuple[list[dict], list[dict]]:
    calls, audit = [], []
    model.eval()
    with torch.no_grad():
        for locus in loci:
            scores, oracle_logit = model.forward_locus(locus, store)
            probabilities = calibrated_probabilities(scores, temperature)
            ranking = sorted(range(len(locus["candidates"])), key=lambda idx: (-float(probabilities[idx]),
                                                                                  locus["candidates"][idx]["pair_key"]))
            winner = ranking[0]
            confidence = float(probabilities[winner])
            selected = locus["candidates"][winner]
            abstain = threshold is not None and confidence < threshold
            call = {"sample": locus["sample"], "modality": locus["modality"], "gene": locus["gene"],
                    "method": "ChampHLA-RefFormer", "allele1": selected["pair"][0],
                    "allele2": selected["pair"][1], "confidence": confidence,
                    "call_status": "no_call" if abstain else "called", "is_callable": "0" if abstain else "1",
                    "abstention_reason": "below_refformer_confidence" if abstain else "",
                    "candidate_count": len(ranking), "candidate_set_oracle": locus.get("candidate_set_oracle", ""),
                    "oracle_probability": float(torch.sigmoid(oracle_logit)), "model_sha256": model_sha256,
                    "imgt_sha256": store.reference_sha256,
                    "embedding_encoder_sha256": store.encoder_sha256}
            if "label" in selected:
                call["is_correct"] = str(int(selected.get("label", 0) == 1 and not abstain))
                call["is_correct_2field"] = call["is_correct"]
            calls.append(call)
            for rank, index in enumerate(ranking, 1):
                candidate = locus["candidates"][index]
                audit.append({"sample": locus["sample"], "modality": locus["modality"], "gene": locus["gene"],
                              "candidate_pair": candidate["pair_key"], "rank": rank,
                              "candidate_probability": float(probabilities[index]),
                              "raw_score": float(scores[index]), "selected": int(index == winner),
                              "supporting_callers": ",".join(token["tool"] for token in candidate["tool_tokens"]
                                                               if token["pair"] == candidate["pair"]),
                              "reference_flags_json": json.dumps(candidate["reference_flags"], sort_keys=True),
                              "feature_mask_json": json.dumps({"global": candidate["global_features"],
                                                                 "tools": [token["features"] for token in candidate["tool_tokens"]]},
                                                                separators=(",", ":")),
                              "label": candidate.get("label", ""), "model_sha256": model_sha256,
                              "imgt_sha256": store.reference_sha256,
                              "embedding_encoder_sha256": store.encoder_sha256})
    return calls, audit


def coverage_thresholds(calls: list[dict]) -> dict:
    values = sorted(float(row["confidence"]) for row in calls)
    output = {}
    for coverage in (0.95, 0.90, 0.80):
        if not values:
            threshold = 1.0
        else:
            retained = max(1, int(math.ceil(coverage * len(values))))
            threshold = values[max(0, len(values) - retained)]
        output[str(int(coverage * 100))] = threshold
    return output


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_bundle(directory: str | Path, model: RefFormerModel, context: dict, store: ReferenceEmbeddingStore,
                temperature: float, thresholds: dict, training_summary: dict, extra: dict | None = None) -> dict:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    weights = directory / "model.pt"
    torch.save({"state_dict": model.state_dict(), "model_config": model.config()}, weights)
    manifest = {"schema_version": SCHEMA_VERSION, "weights": weights.name,
                "weights_sha256": file_sha256(weights), "reference_sha256": store.reference_sha256,
                "architecture": model.config(),
                "embedding_encoder_sha256": store.encoder_sha256, "temperature": temperature,
                "embedding_cache_sha256": store.cache_sha256,
                "tokenizer": store.metadata.get("tokenizer", {}),
                "reference_manifest": store.metadata.get("reference_manifest", {}),
                "coverage_thresholds": thresholds, "runtime_context": context,
                "training_summary": training_summary, "runtime_default_changed": False,
                "public_name": "ChampHLA-RefFormer", **(extra or {})}
    unsigned = json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()
    manifest["bundle_sha256"] = hashlib.sha256(unsigned).hexdigest()
    (directory / "model_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def load_bundle(directory: str | Path, store: ReferenceEmbeddingStore, device: str = "cpu"):
    directory = Path(directory)
    manifest = json.loads((directory / "model_manifest.json").read_text())
    unsigned = dict(manifest)
    checksum = unsigned.pop("bundle_sha256")
    observed = hashlib.sha256(json.dumps(unsigned, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if checksum != observed or file_sha256(directory / manifest["weights"]) != manifest["weights_sha256"]:
        raise ValueError("RefFormer model bundle checksum mismatch")
    if manifest["reference_sha256"] != store.reference_sha256:
        raise ValueError("RefFormer/IMGT reference checksum mismatch")
    if manifest.get("embedding_cache_sha256") != store.cache_sha256:
        raise ValueError("RefFormer embedding cache checksum mismatch")
    payload = torch.load(directory / manifest["weights"], map_location=device)
    config = payload["model_config"]
    model = RefFormerModel(reference_dim=int(config["reference_dim"]), hidden_dim=int(config["hidden_dim"]),
                           layers=int(config["layers"]), heads=int(config["heads"]),
                           dropout=float(config.get("dropout", 0.1)),
                           no_imgt=bool(config.get("no_imgt", False)))
    model.load_state_dict(payload["state_dict"])
    model.to(torch.device(device))
    return model, manifest
