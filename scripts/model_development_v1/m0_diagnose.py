"""N-only baselines on the frozen Stage B donors, plus distribution diagnostics."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ["STAGEB_NJOBS"] = "4"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "stage_b"))
from score_stage_b import bootstrap_weights, components, regression_stats, source_table  # noqa: E402
from stageb_lib import derive_seed  # noqa: E402
from train_stage_b import predict_model, training_frame  # noqa: E402

SPLIT = REPO / "data" / "splits_v1"
OUT = REPO / "runs" / "model_development_v1"
FIG = OUT / "figures"
SEEDS = [17, 29, 42, 71, 101, 137, 271, 314, 577, 911]
TL_ARMS = ["TL-N0", "TL-FN", "TL-JN", "TL-B"]
QUANTILES = [0.1, 0.25, 0.5, 0.75, 0.9]


def summarize(values: pd.Series) -> dict:
    array = values.to_numpy(np.float64)
    row = {
        "n": int(array.size),
        "mean": float(np.mean(array)),
        "std": float(np.std(array, ddof=1)) if array.size > 1 else np.nan,
    }
    for q in QUANTILES:
        row[f"p{int(q * 100)}"] = float(np.quantile(array, q))
    return row


def load_c02() -> pd.DataFrame:
    frames = [pd.read_csv(path) for path in sorted((REPO / "runs" / "stage_c" / "c02" / "pred").glob("*.csv"))]
    pred = pd.concat(frames, ignore_index=True)
    return pred[pred.arm.isin(TL_ARMS)].copy()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    features = np.load(REPO / "runs" / "stage_b" / "features_dev.npz", allow_pickle=True)
    index = {str(gid): i for i, gid in enumerate(features["graph_id"])}
    desc = features["desc"].astype(np.float64)
    fp = features["fp"]
    usable = pd.read_csv(SPLIT / "outer_members.csv.gz")
    usable = usable[usable.role.isin(["test", "train_eligible"])].copy()
    usable["assigned_fold"] = usable.assigned_fold.astype(int)
    test = pd.read_csv(SPLIT / "test_labels.csv.gz")
    test["assigned_fold"] = test.assigned_fold.astype(int)
    members = pd.read_csv(SPLIT / "sampling_members.csv.gz")
    donors = members[members.block.eq("B01") & members.arm.eq("BN") & members.component.eq("N")]
    rows = []
    dist_rows = []
    for fold in range(5):
        evaluation = test[test.assigned_fold.eq(fold)].copy()
        dist_rows.append({"record": "test_labels", "fold": fold, "seed": -1, "component": "TEST", **summarize(evaluation.y)})
        for seed in SEEDS:
            selected = donors[donors.fold.eq(fold) & donors.seed.eq(seed)][["Graph_ID", "component"]]
            if selected.Graph_ID.nunique() != 50:
                raise RuntimeError(f"fold {fold} seed {seed} has {selected.Graph_ID.nunique()} N donors")
            train = training_frame(usable, selected, fold, "B01", {}, set())
            for model in ["Median", "Ridge", "RF"]:
                rf_seed = derive_seed("Caco-2:AB", fold, seed, "N-only", model)
                y_pred = predict_model(model, train, evaluation, desc, fp, index, rf_seed)
                part = evaluation[["label_id", "Graph_ID", "Study_key", "y"]].copy()
                part["y_pred"] = y_pred
                part["fold"] = fold
                part["seed"] = seed
                part["model"] = f"N-only {model}"
                rows.append(part)
    pred = pd.concat(rows, ignore_index=True)
    pred = pred.rename(columns={"Graph_ID": "graph_id", "Study_key": "study_id"})
    pred["abs_error"] = (pred.y - pred.y_pred).abs()
    pred["signed_error"] = pred.y_pred - pred.y
    pred.to_csv(OUT / "m0_predictions.csv.gz", index=False)

    c02 = load_c02()
    c02["abs_error"] = (c02.y - c02.y_pred).abs()
    c02["signed_error"] = c02.y_pred - c02.y
    c02["model"] = c02.arm
    combined = pd.concat(
        [
            pred[["model", "fold", "seed", "label_id", "graph_id", "study_id", "y", "y_pred", "abs_error", "signed_error"]],
            c02[["model", "fold", "seed", "label_id", "graph_id", "study_id", "y", "y_pred", "abs_error", "signed_error"]],
        ],
        ignore_index=True,
    )
    study_component = components(test)
    studies, weights = bootstrap_weights(study_component, np.random.default_rng(271828))
    metric_rows = []
    seed_frames = []
    for model, part in combined.groupby("model"):
        seed_mean = part.groupby("label_id", as_index=False).agg(
            abs_error=("abs_error", "mean"),
            signed_error=("signed_error", "mean"),
            y=("y", "first"),
            y_pred=("y_pred", "mean"),
            study_id=("study_id", "first"),
            graph_id=("graph_id", "first"),
        )
        seed_mean["model"] = model
        seed_frames.append(seed_mean)
        per_seed = [regression_stats(chunk) for _, chunk in part.groupby("seed")]
        stat = pd.DataFrame(per_seed)
        y = seed_mean.y.to_numpy()
        yhat = seed_mean.y_pred.to_numpy()
        metric_rows.append(
            {
                "model": model,
                "mae_source": float(source_table(seed_mean, "abs_error").mean()),
                "signed_source": float(source_table(seed_mean, "signed_error").mean()),
                "rmse_seed_mean": float(stat.rmse.mean()),
                "r2_seed_mean": float(stat.r2.mean()),
                "spearman_seed_mean": float(stat.spearman.mean()),
                "pred_std_over_y_std": float(np.std(yhat, ddof=1) / np.std(y, ddof=1)),
                "n_labels": int(seed_mean.label_id.nunique()),
            }
        )
    metrics = pd.DataFrame(metric_rows).sort_values("mae_source")
    seed_mean = pd.concat(seed_frames, ignore_index=True)
    draws = {}
    points = {}
    for model, chunk in seed_mean.groupby("model"):
        vector = source_table(chunk, "abs_error").loc[studies].to_numpy(np.float64)
        draws[model] = (weights * vector).sum(1) / weights.sum(1)
        points[model] = float(vector.mean())
    delta_rows = []
    for left, right in [
        ("N-only Median", "TL-N0"),
        ("N-only Ridge", "TL-N0"),
        ("N-only RF", "TL-N0"),
        ("TL-N0", "TL-FN"),
        ("TL-N0", "TL-JN"),
        ("TL-N0", "TL-B"),
    ]:
        draw = draws[left] - draws[right]
        low, high = np.quantile(draw, [0.025, 0.975])
        delta_rows.append(
            {
                "left": left,
                "right": right,
                "estimate": points[left] - points[right],
                "ci_low": float(low),
                "ci_high": float(high),
                "interval_excludes_zero": bool(low > 0 or high < 0),
                "n_studies": int(len(study_component)),
                "n_components": int(study_component.max()) + 1,
            }
        )
    deltas = pd.DataFrame(delta_rows)
    for component, arm in [("B", "B"), ("S", "BS"), ("M", "BM"), ("N", "BN")]:
        selected = members[members.block.eq("B01") & members.arm.eq(arm) & members.component.eq(component)]
        for (fold, seed), part in selected.groupby(["fold", "seed"]):
            graphs = set(part.Graph_ID)
            labels = usable[usable.Graph_ID.isin(graphs) & usable.assigned_fold.ne(int(fold))]
            dist_rows.append({"record": "train_labels", "fold": int(fold), "seed": int(seed), "component": component, **summarize(labels.y)})
    for model, part in combined.groupby("model"):
        for (fold, seed), chunk in part.groupby(["fold", "seed"]):
            row = {"record": "prediction", "fold": int(fold), "seed": int(seed), "component": model, **summarize(chunk.y_pred)}
            row["pred_std_over_y_std"] = row["std"] / float(np.std(chunk.y, ddof=1))
            dist_rows.append(row)
    bins = pd.cut(seed_mean.y, [-np.inf, -6.5, -5.5, np.inf], labels=["le_-6.5", "-6.5_to_-5.5", "gt_-5.5"])
    regional = (
        seed_mean.assign(region=bins)
        .groupby(["model", "region"], observed=True)
        .agg(n=("label_id", "size"), signed_mean=("signed_error", "mean"), mae=("abs_error", "mean"))
        .reset_index()
    )
    metrics.to_csv(OUT / "m0_metrics.csv", index=False)
    deltas.to_csv(OUT / "m0_paired_deltas.csv", index=False)
    pd.DataFrame(dist_rows).to_csv(OUT / "m0_distributions.csv", index=False)
    regional.to_csv(OUT / "m0_regional_error.csv", index=False)
    order = ["N-only Median", "N-only Ridge", "N-only RF", "TL-N0", "TL-FN", "TL-JN", "TL-B"]
    plot_df = metrics.set_index("model").loc[order]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.bar(range(len(plot_df)), plot_df.mae_source, color="#4C78A8")
    ax.set_xticks(range(len(plot_df)))
    ax.set_xticklabels(order, rotation=30, ha="right")
    ax.set_ylabel("Source-equal MAE")
    ax.set_title("N-only baselines and C02 transfer models")
    fig.tight_layout()
    fig.savefig(FIG / "m0_mae_comparison.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4), sharex=True, sharey=True)
    for ax, model in zip(axes, ["N-only Median", "TL-N0"]):
        part = seed_mean[seed_mean.model.eq(model)]
        ax.scatter(part.y, part.y_pred, s=12, alpha=0.7, color="#4C78A8")
        ax.plot([-8, -3], [-8, -3], color="#333333", lw=0.8)
        ax.set_title(model)
        ax.set_xlabel("Experimental log10(Papp [cm/s])")
        ax.set_ylabel("Predicted")
    fig.tight_layout()
    fig.savefig(FIG / "m0_pred_vs_obs.png", dpi=160)
    plt.close(fig)
    (OUT / "m0_summary.json").write_text(json.dumps({"metrics": metrics.to_dict("records"), "deltas": deltas.to_dict("records")}, indent=2))
    print(metrics.to_string(index=False))
    print(deltas.to_string(index=False))


if __name__ == "__main__":
    main()
