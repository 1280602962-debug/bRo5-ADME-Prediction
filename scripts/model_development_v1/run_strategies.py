"""Fixed-parameter G2/G3 refits on the same evaluation-C test labels."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from model_lib import FIXED, OUT, load_store
from run_final import predict_split

PARTS = OUT / "final_parts"


def pool_for(dev: pd.DataFrame, strategy: str, outer: int) -> pd.DataFrame:
    part = dev[dev.assigned_fold.ne(outer)].copy()
    if strategy == "G2":
        part = part[part.cohort.astype(str).str.startswith("nonmacro")]
    elif strategy == "G3":
        part = part[part.cohort.eq("nonmacro_500_700")]
    else:
        raise KeyError(strategy)
    return part


def main() -> None:
    store = load_store()
    dev = store["dev"]
    dev = dev[dev.g1].copy()
    dev["assigned_fold"] = dev.assigned_fold.astype(int)
    test = pd.read_csv(Path(__file__).resolve().parents[2] / "data" / "splits_v1" / "test_labels.csv.gz")
    frames = []
    for strategy in ["G2", "G3"]:
        for outer in range(5):
            train = pool_for(dev, strategy, outer)
            evaluation = test[test.assigned_fold.eq(outer)].copy()
            for model, feature in [("RF", "F05"), ("RF", "F04")]:
                key = f"C_{strategy}_o{outer}_fixed_{model}_{feature}"
                path = PARTS / f"{key}.csv"
                if path.exists():
                    frames.append(pd.read_csv(path))
                    continue
                pred = predict_split(store, train, evaluation, model, feature, FIXED[model], key)
                pred["evaluation"] = "C"
                pred["outer_fold"] = outer
                pred["parameter_source"] = f"fixed_{strategy}"
                pred["exploratory"] = False
                pred.to_csv(path, index=False)
                frames.append(pred)
                print(key, "n_train", int(pred.n_train_labels.iloc[0]), "train_mae", round(float(pred.train_mae_label.iloc[0]), 4), flush=True)
    extra = pd.concat(frames, ignore_index=True)
    current = pd.read_csv(OUT / "oof_predictions.csv.gz")
    merged = pd.concat([current, extra], ignore_index=True)
    merged = merged.drop_duplicates(["evaluation", "outer_fold", "model", "feature", "parameter_source", "label_id"])
    merged.to_csv(OUT / "oof_predictions.csv.gz", index=False)
    print("STRATEGY_DONE", len(extra))


if __name__ == "__main__":
    main()
