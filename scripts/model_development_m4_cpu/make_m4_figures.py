"""Figures for completed M4 modules. Solvent figures are omitted when M4.4 did not run."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

OUT = Path(__file__).resolve().parents[2] / "runs" / "model_development_m4_cpu"
FIG = OUT / "figures"


def save(fig, name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / name, dpi=160)
    plt.close(fig)


def figure_repeats() -> None:
    table = pd.read_csv(OUT / "measurement_variability_summary.csv")
    table = table[table.cohort.ne("ALL_DEVELOPMENT")]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    x = np.arange(len(table))
    ax.bar(x - 0.18, table.database_reentry_groups, width=0.36, label="Database re-entry", color="#4C78A8")
    ax.bar(x + 0.18, table.independent_repeat_groups, width=0.36, label="Independent repeat", color="#F58518")
    ax.set_xticks(x, table.cohort, rotation=15, ha="right")
    ax.set_ylabel("Multi-record groups")
    ax.set_title("Within-study repeats by cohort")
    ax.legend(frameon=False)
    save(fig, "M4-1_repeatability.png")


def figure_coverage() -> None:
    summary = pd.read_csv(OUT / "conformer_generation_summary.csv")
    eligibility = pd.read_csv(OUT / "conformer_eligibility.csv")
    counts = summary.sampling_status.value_counts()
    excluded = int((~eligibility.eligible_main.astype(str).str.lower().eq("true")).sum())
    labels = list(counts.index) + ["excluded"]
    values = list(counts.to_numpy()) + [excluded]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.bar(labels, values, color="#54A24B")
    ax.set_ylabel("Graphs")
    ax.set_title("Conformer sampling status")
    ax.tick_params(axis="x", rotation=15)
    save(fig, "M4-2_feature_coverage.png")


def figure_geometry() -> None:
    ensemble = pd.read_csv(OUT / "l1_ensemble_features.csv.gz")
    groups = pd.read_csv(Path(__file__).resolve().parents[2] / "data" / "modeling_v1" / "molecule_groups.csv.gz")
    groups["Graph_ID"] = groups.Graph_ID.astype(str)
    ensemble["Graph_ID"] = ensemble.Graph_ID.astype(str)
    merged = ensemble.merge(groups[["Graph_ID", "MW", "cLogP", "HBD", "HBA"]], on="Graph_ID", how="left")
    # TPSA is not in molecule_groups; use eligibility.
    eligibility = pd.read_csv(OUT / "conformer_eligibility.csv")
    eligibility["Graph_ID"] = eligibility.Graph_ID.astype(str)
    merged = merged.merge(eligibility[["Graph_ID", "TPSA"]], on="Graph_ID", how="left")
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))
    pairs = [
        ("MW", "ens_Rgyr_median", "Radius of gyration"),
        ("TPSA", "ens_PolarSASA_median", "Polar SASA"),
        ("HBD", "ens_IMHB_mean", "Intramolecular H-bonds"),
    ]
    for ax, (xname, yname, title) in zip(axes, pairs):
        ax.scatter(merged[xname], merged[yname], s=8, alpha=0.35, color="#4C78A8", linewidths=0)
        ax.set_xlabel(xname)
        ax.set_ylabel(title)
    save(fig, "M4-3_geometry_vs_2d.png")


def figure_ablation() -> None:
    metrics = pd.read_csv(OUT / "m4_metrics.csv")
    part = metrics[metrics.evaluation.eq("C") & metrics.model.eq("RF") & metrics.parameter_source.eq("fixed_rf")]
    keep = ["P0_ensemble", "P1_ensemble", "P2_ion", "P3_static", "P3_ensemble", "P4_ensemble", "P4_both", "P5_both"]
    part = part[part.arm.isin(keep)].sort_values("mae_source")
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    ax.barh(part.arm, part.mae_source, color="#4C78A8")
    ax.set_xlabel("Source-equal MAE")
    ax.set_title("Matched evaluation C")
    save(fig, "M4-4_ablation.png")


def figure_deltas() -> None:
    deltas = pd.read_csv(OUT / "m4_paired_deltas.csv")
    estimated = deltas[deltas.status.eq("estimated")].copy()
    if estimated.empty:
        return
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    y = np.arange(len(estimated))
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.errorbar(
        estimated.estimate,
        y,
        xerr=[estimated.estimate - estimated.ci95_low, estimated.ci95_high - estimated.estimate],
        fmt="o",
        color="#4C78A8",
        label="95%",
    )
    if "ci_bonferroni_low" in estimated.columns:
        primary = estimated[estimated.role.eq("primary")]
        ax.errorbar(
            primary.estimate,
            np.arange(len(estimated))[estimated.role.eq("primary").to_numpy()],
            xerr=[primary.estimate - primary.ci_bonferroni_low, primary.ci_bonferroni_high - primary.estimate],
            fmt="none",
            ecolor="#F58518",
            label="Bonferroni",
        )
    ax.set_yticks(y, estimated.contrast)
    ax.set_xlabel("Source-equal MAE difference")
    ax.legend(frameon=False)
    save(fig, "M4-5_paired_deltas.png")


def figure_conformer_image() -> None:
    sdf_dir = OUT / "representative_conformers"
    paths = sorted(sdf_dir.glob("*.sdf")) if sdf_dir.exists() else []
    if not paths:
        return
    mols = []
    legends = []
    for path in paths[:4]:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
        mol = next((item for item in supplier if item is not None), None)
        if mol is None:
            continue
        mols.append(mol)
        legends.append(path.stem)
    if not mols:
        return
    image = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=2, subImgSize=(280, 220), returnPNG=True)
    FIG.mkdir(parents=True, exist_ok=True)
    (FIG / "M4-7_representative_conformers.png").write_bytes(image)


def figure_remaining() -> None:
    coverage = pd.read_csv(OUT / "m4_coverage_metrics.csv")
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    ax.plot(coverage.outer_fold, coverage.test_labels, marker="o", label="All test labels")
    ax.plot(coverage.outer_fold, coverage.ensemble_test_labels, marker="o", label="Ensemble available")
    ax.plot(coverage.outer_fold, coverage.fallback_uses_enhanced, marker="o", label="Fallback uses 3D")
    ax.set_xlabel("Outer fold")
    ax.set_ylabel("Test labels")
    ax.set_title("Coverage of the pre-specified fallback")
    ax.legend(frameon=False)
    save(fig, "M4-8_coverage_and_gaps.png")


def main() -> None:
    figure_repeats()
    figure_coverage()
    figure_geometry()
    figure_ablation()
    figure_deltas()
    figure_conformer_image()
    figure_remaining()
    print("figures", sorted(path.name for path in FIG.glob("*.png")))


if __name__ == "__main__":
    main()
