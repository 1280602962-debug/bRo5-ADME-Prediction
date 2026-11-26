"""Figures for the completed M5 tables."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

OUT = Path(__file__).resolve().parents[2] / "runs" / "model_development_m5"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def seed_figure() -> None:
    frame = pd.read_csv(OUT / "seed_sensitivity.csv")
    frame = frame[frame.summary.isin(["one_seed", "seed_mean_absolute_error"])]
    labels = {"R3-R0": "ensemble − F04", "R3-R1": "ensemble − F05", "R3-R2": "ensemble − static"}
    figure, axis = plt.subplots(figsize=(7.2, 4.2))
    for index, (contrast, part) in enumerate(frame.groupby("contrast", sort=False)):
        for row_index, row in enumerate(part.itertuples(index=False)):
            y = index + row_index * 0.16
            axis.plot([row.ci95_low, row.ci95_high], [y, y], color="#4C78A8", lw=1.5)
            axis.plot(row.estimate, y, "o", color="#4C78A8")
        axis.text(-0.055, index + 0.24, labels.get(contrast, contrast), ha="right", va="center", fontsize=9)
    axis.axvline(0, color="#333333", lw=0.8)
    axis.set_yticks([])
    axis.set_xlabel("Paired source-MAE difference, candidate minus reference")
    axis.set_title("Shared-seed sensitivity; intervals are not the M4 family")
    figure.tight_layout()
    figure.savefig(FIG / "M5_seed_sensitivity.png", dpi=200)
    plt.close(figure)


def group_figure() -> None:
    path = OUT / "feature_group_ablation.csv"
    if not path.exists():
        return
    frame = pd.read_csv(path)
    figure, axis = plt.subplots(figsize=(6.4, 3.6))
    for index, row in enumerate(frame.itertuples(index=False)):
        axis.plot([row.ci95_low, row.ci95_high], [index, index], color="#F58518", lw=1.5)
        axis.plot(row.estimate, index, "o", color="#F58518")
        axis.text(row.ci95_high + 0.001, index, row.contrast, va="center", fontsize=9)
    axis.axvline(0, color="#333333", lw=0.8)
    axis.set_yticks([])
    axis.set_xlabel("Exploratory paired source-MAE difference versus F04")
    figure.tight_layout()
    figure.savefig(FIG / "M5_group_ablation.png", dpi=200)
    plt.close(figure)


def reliability_figure() -> None:
    path = OUT / "reliability_metrics.csv"
    if not path.exists():
        return
    frame = pd.read_csv(path)
    kept = frame[frame.method.isin(["F05_full_coverage", "ridge_crossfit", "nearest_neighbor", "tree_dispersion", "f04_f05_gap"])]
    figure, axis = plt.subplots(figsize=(6.6, 4.0))
    for method, part in kept.groupby("method"):
        axis.plot(part.realized_molecule_coverage, part.mae_source, marker="o", label=method)
    axis.set_xlabel("Realized molecule coverage")
    axis.set_ylabel("Source MAE on the retained subset")
    axis.set_title("Subset error is not full-coverage accuracy")
    axis.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(FIG / "M5_reliability_coverage.png", dpi=200)
    plt.close(figure)


if __name__ == "__main__":
    seed_figure()
    group_figure()
    reliability_figure()
