"""Rescore the archived C02 predictions. Does not refit D-MPNN or alter C01/C02 outputs."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stage_b"))
from score_stage_b import (  # noqa: E402
    bootstrap_weights,
    components,
    graph_mean,
    regression_stats,
    source_table,
)

REPO = Path(__file__).resolve().parents[2]
SPLIT = REPO / "data" / "splits_v1"
RUN = REPO / "runs" / "stage_c"
PRED = RUN / "c02" / "pred"
SEEDS = [17, 29, 42, 71, 101, 137, 271, 314, 577, 911]
ARMS = ["TL-B", "TL-N0", "TL-JS", "TL-JM", "TL-JN", "TL-FS", "TL-FM", "TL-FN"]
CONTRASTS = [
    ("TL01_FN_minus_JN", "TL-FN", "TL-JN", True),
    ("FN_minus_B", "TL-FN", "TL-B", True),
    ("FN_minus_N0", "TL-FN", "TL-N0", True),
    ("FM_minus_JM", "TL-FM", "TL-JM", True),
    ("FS_minus_JS", "TL-FS", "TL-JS", True),
    ("JN_minus_B", "TL-JN", "TL-B", False),
    ("JM_minus_B", "TL-JM", "TL-B", False),
    ("JS_minus_B", "TL-JS", "TL-B", False),
    ("N0_minus_B", "TL-N0", "TL-B", False),
    ("JN_minus_N0", "TL-JN", "TL-N0", False),
    ("JM_minus_N0", "TL-JM", "TL-N0", False),
    ("JS_minus_N0", "TL-JS", "TL-N0", False),
]
DONOR = {
    "TL-B": "B",
    "TL-N0": "N",
    "TL-FN": "N",
    "TL-JN": "N",
    "TL-FS": "S",
    "TL-JS": "S",
    "TL-FM": "M",
    "TL-JM": "M",
}
ARM_MEMBER = {"B": ("B", "B"), "S": ("BS", "S"), "M": ("BM", "M"), "N": ("BN", "N")}


def load_predictions() -> pd.DataFrame:
    frames = [pd.read_csv(path) for path in sorted(PRED.glob("*.csv"))]
    pred = pd.concat(frames, ignore_index=True)
    pred["abs_error"] = (pred.y - pred.y_pred).abs()
    pred["signed_error"] = pred.y_pred - pred.y
    return pred


def coverage(pred: pd.DataFrame, test: pd.DataFrame) -> dict:
    expected = set(test.label_id)
    problems = []
    if len(pred) != 42400:
        problems.append(f"rows {len(pred)} != 42400")
    for arm, part in pred.groupby("arm"):
        if len(part) != 5300:
            problems.append(f"{arm} rows {len(part)} != 5300")
        for seed, chunk in part.groupby("seed"):
            ids = set(chunk.label_id)
            if ids != expected:
                problems.append(f"{arm} seed {seed} label mismatch missing {len(expected - ids)} extra {len(ids - expected)}")
            if chunk.label_id.duplicated().any():
                problems.append(f"{arm} seed {seed} duplicate labels")
    merged = pred[["label_id", "y", "study_id", "graph_id"]].drop_duplicates("label_id").merge(
        test[["label_id", "y", "Study_key", "Graph_ID", "corrected_task", "analysis_set"]],
        on="label_id",
        how="outer",
        indicator=True,
    )
    if not (merged._merge == "both").all():
        problems.append("prediction keys do not match the frozen test labels")
    both = merged[merged._merge.eq("both")]
    if not np.allclose(both.y_x, both.y_y, atol=1e-8, equal_nan=False):
        problems.append("prediction y differs from the frozen test y")
    if not both.study_id.eq(both.Study_key).all():
        problems.append("study_id differs from Study_key")
    if not both.graph_id.eq(both.Graph_ID).all():
        problems.append("graph_id differs from Graph_ID")
    tasks = sorted(test.corrected_task.unique())
    sets = sorted(test.analysis_set.unique())
    if tasks != ["Caco-2:AB"] or sets != ["caco2_ab_main"]:
        problems.append(f"unexpected endpoint {tasks} {sets}")
    holdout = set(pd.read_csv(SPLIT / "holdout_members.csv").Graph_ID)
    if set(pred.graph_id) & holdout:
        problems.append("holdout graph entered C02 predictions")
    return {
        "n_rows": int(len(pred)),
        "n_files": int(len(list(PRED.glob('*.csv')))),
        "tasks": tasks,
        "analysis_sets": sets,
        "problems": problems,
        "fold_counts": {int(k): int(v) for k, v in test.groupby("assigned_fold").size().items()},
    }


def training_audit() -> pd.DataFrame:
    rows = []
    path = RUN / "c02" / "training_records.jsonl"
    for line in path.read_text().splitlines():
        if line.strip():
            rows.append(json.loads(line))
    frame = pd.DataFrame(rows)
    summary = []
    for arm, part in frame.groupby("arm"):
        before = part.get("encoder_hash_before")
        after = part.get("encoder_hash_after")
        changed = 0
        if before is not None and after is not None:
            mask = before.notna() & after.notna()
            changed = int((before[mask] != after[mask]).sum())
        summary.append(
            {
                "arm": arm,
                "n": int(len(part)),
                "encoder_frozen": sorted(part.encoder_frozen.dropna().unique().tolist()),
                "encoder_trainable": sorted({None if pd.isna(v) else int(v) for v in part.encoder_trainable}) if "encoder_trainable" in part else [],
                "steps": sorted(int(v) for v in part.steps.unique()),
                "scaler": sorted(part.scaler.fillna("none").astype(str).unique()) if "scaler" in part else [],
                "n_reused": int(part.reused_from.notna().sum()) if "reused_from" in part else 0,
                "unique_hash_before": int(before.nunique(dropna=True)) if before is not None else 0,
                "unique_hash_after": int(after.nunique(dropna=True)) if after is not None else 0,
                "n_hash_changed": changed,
            }
        )
    return frame, pd.DataFrame(summary)


def nearest_tanimoto(pred: pd.DataFrame) -> pd.DataFrame:
    store = np.load(REPO / "runs" / "stage_b" / "features_dev.npz", allow_pickle=True)
    ids = [str(v) for v in store["graph_id"]]
    index = {gid: i for i, gid in enumerate(ids)}
    fp = store["fp"]
    members = pd.read_csv(SPLIT / "sampling_members.csv.gz")
    members = members[members.block.eq("B01")]
    test = pd.read_csv(SPLIT / "test_labels.csv.gz")
    rows = []
    for arm, donor in DONOR.items():
        member_arm, component = ARM_MEMBER[donor]
        selected = members[members.arm.eq(member_arm) & members.component.eq(component)]
        for (fold, seed), part in selected.groupby(["fold", "seed"]):
            train_ids = [gid for gid in part.Graph_ID.unique() if gid in index]
            test_ids = [gid for gid in test.loc[test.assigned_fold.eq(fold), "Graph_ID"].unique() if gid in index]
            train_fp = fp[[index[gid] for gid in train_ids]].astype(np.float32)
            query_fp = fp[[index[gid] for gid in test_ids]].astype(np.float32)
            intersection = query_fp @ train_fp.T
            denom = query_fp.sum(1, keepdims=True) + train_fp.sum(1) - intersection
            nearest = (intersection / np.maximum(denom, 1e-12)).max(1)
            rows.append({"arm": arm, "donor": donor, "fold": int(fold), "seed": int(seed), "nearest_mean": float(nearest.mean())})
    detail = pd.DataFrame(rows)
    out = detail.groupby(["arm", "donor"], as_index=False).nearest_mean.mean()
    return out


def diagnostics(seed_mean: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    for arm, part in seed_mean.groupby("arm"):
        y = part.y.to_numpy()
        yhat = part.y_pred.to_numpy()
        rows.append(
            {
                "arm": arm,
                "y_std": float(np.std(y, ddof=1)),
                "pred_std": float(np.std(yhat, ddof=1)),
                "pred_std_over_y_std": float(np.std(yhat, ddof=1) / np.std(y, ddof=1)),
                "y_p10": float(np.quantile(y, 0.1)),
                "y_p90": float(np.quantile(y, 0.9)),
                "pred_p10": float(np.quantile(yhat, 0.1)),
                "pred_p90": float(np.quantile(yhat, 0.9)),
                "signed_source": float(source_table(part, "signed_error").mean()),
            }
        )
    bins = pd.cut(seed_mean.y, [-np.inf, -6.5, -5.5, np.inf], labels=["le_-6.5", "-6.5_to_-5.5", "gt_-5.5"])
    regional = (
        seed_mean.assign(region=bins)
        .groupby(["arm", "region"], observed=True)
        .agg(n=("label_id", "size"), signed_mean=("signed_error", "mean"), mae=("abs_error", "mean"))
        .reset_index()
    )
    wide = seed_mean.pivot(index="label_id", columns="arm", values="abs_error")
    corr = wide.corr()
    corr_rows = corr.reset_index().melt(id_vars="arm", var_name="other", value_name="abs_error_pearson")
    return pd.DataFrame(rows), regional, corr_rows


def write_report(coverage_info: dict, pooled: pd.DataFrame, deltas: pd.DataFrame, audit: pd.DataFrame, compression: pd.DataFrame, published_delta: pd.DataFrame) -> None:
    def fmt(frame: pd.DataFrame) -> str:
        return "```\n" + frame.to_csv(index=False) + "```"

    n0 = audit[audit.arm.eq("TL-N0")].iloc[0]
    text = f"""# 阶段 C 统计归档

日期：2026-11-25。父提交 `25cbde3`。本文件只归档已经完成的 C02 预测，不重新训练 D-MPNN，也不补跑缺失的 C01。

预测终点是 Caco-2 A→B 的 log10(Papp [cm/s])。冻结测试表的 `corrected_task` 为 Caco-2:AB，`analysis_set` 为 caco2_ab_main。这里的 y 不是 log10 活性，也不是摩尔浓度。

## 完成状态

C02 八个策略已经执行并入库：预测文件 {coverage_info['n_files']} 个，预测行 {coverage_info['n_rows']}。每个策略 5300 行，对应 530 条正式测试标签和 10 个种子。覆盖检查问题数：{len(coverage_info['problems'])}。

C01 只有部分预测有记录。TL-B 的 50 条训练记录标明复用集群上的 C01 D-MPNN B 臂预测。仓库中没有 C01 XGBoost，也没有 C01 的 S、M、N 和 ALL_DEV。阶段 C 的完成范围是 C02 的统计归档和独立复核。缺少训练源码，因此本归档不表示这些网络训练已经被独立复现。

TL01 仍是 TL-FN 减 TL-JN，检验家族仍是 FTRANSFER。P01 不在本阶段修改。阴性结果保留。

## TL-N0 的现有证据

TL-N0 共 {int(n0.n)} 条记录。`encoder_frozen` 为 {n0.encoder_frozen}，`encoder_trainable` 为 {n0.encoder_trainable}，优化步数为 {n0.steps}，scaler 为 {n0.scaler}。训练前哈希有 {int(n0.unique_hash_before)} 个不同值，训练后哈希有 {int(n0.unique_hash_after)} 个不同值，哈希发生变化的记录为 {int(n0.n_hash_changed)}。这些是随机初始化并更新编码器的现有支持证据。冻结头策略的编码器哈希在训练前后不变。没有阶段 C 训练源码，不能把这一证据写成已经独立复现了网络训练。

## 重算指标

来源等权 MAE 使用阶段 B 已复核的聚合顺序：先对种子平均绝对误差，再在研究内对分子图平均，然后对研究等权。RMSE、R²、Spearman 和 Pearson 在每个种子的完整 530 条预测上计算，再对种子取均值。有符号误差是预测值减实验值。

{fmt(pooled)}

与 `25cbde3` 已发布汇总表的来源等权 MAE 差值：

{fmt(published_delta)}

## 配对差值

Bootstrap 为 10000 次，生成器种子 271828，重抽样研究–SC 连通组件，沿用阶段 B 的权重。估计值是左侧来源等权 MAE 减右侧。区间包含 0 不是等价。

{fmt(deltas)}

## 预测诊断

预测标准差与实验标准差之比、分区域有符号误差和策略间绝对误差相关写在 `prediction_diagnostics.csv`、`regional_signed_error.csv` 和 `error_correlation.csv`。训练供体与测试图的最近邻 Tanimoto 使用阶段 B 冻结 ECFP4 和 B01 供体，写在 `nearest_tanimoto.csv`。这些诊断不修改 C02 主比较。

{fmt(compression)}
"""
    (RUN / "STAGE_C_REPORT.md").write_text(text)


def main() -> None:
    test = pd.read_csv(SPLIT / "test_labels.csv.gz")
    pred = load_predictions()
    coverage_info = coverage(pred, test)
    if coverage_info["problems"]:
        raise RuntimeError("; ".join(coverage_info["problems"]))
    pred.to_csv(RUN / "c02" / "oof_predictions.csv.gz", index=False)
    study_component = components(test)
    studies, weights = bootstrap_weights(study_component, np.random.default_rng(271828))
    pooled_rows = []
    seed_frames = []
    for arm, part in pred.groupby("arm"):
        seed_mean = part.groupby("label_id", as_index=False).agg(
            abs_error=("abs_error", "mean"),
            signed_error=("signed_error", "mean"),
            y=("y", "first"),
            y_pred=("y_pred", "mean"),
            study_id=("study_id", "first"),
            graph_id=("graph_id", "first"),
            fold=("fold", "first"),
        )
        seed_mean["arm"] = arm
        seed_frames.append(seed_mean)
        per_seed = []
        for _, chunk in part.groupby("seed"):
            stat = regression_stats(chunk)
            stat["pearson"] = float(pearsonr(chunk.y, chunk.y_pred).statistic)
            per_seed.append(stat)
        stat = pd.DataFrame(per_seed)
        pooled_rows.append(
            {
                "block": "C02",
                "model": "D-MPNN",
                "arm": arm,
                "mae_source": float(source_table(seed_mean, "abs_error").mean()),
                "mae_graph": graph_mean(seed_mean, "abs_error"),
                "mae_label": float(seed_mean.abs_error.mean()),
                "signed_source": float(source_table(seed_mean, "signed_error").mean()),
                "rmse_seed_mean": float(stat.rmse.mean()),
                "r2_seed_mean": float(stat.r2.mean()),
                "spearman_seed_mean": float(stat.spearman.mean()),
                "pearson_seed_mean": float(stat.pearson.mean()),
                "n_labels": int(seed_mean.label_id.nunique()),
                "n_graphs": int(seed_mean.graph_id.nunique()),
                "n_studies": int(seed_mean.study_id.nunique()),
                "n_rows": int(len(part)),
            }
        )
    pooled = pd.DataFrame(pooled_rows).sort_values("arm")
    seed_mean = pd.concat(seed_frames, ignore_index=True)
    point = {row.arm: row.mae_source for row in pooled.itertuples()}
    arm_draws = {}
    for arm, chunk in seed_mean.groupby("arm"):
        vector = source_table(chunk, "abs_error").loc[studies].to_numpy(np.float64)
        arm_draws[arm] = (weights * vector).sum(1) / weights.sum(1)
    n_comp = int(study_component.max()) + 1
    largest = int(study_component.value_counts().iloc[0])
    delta_rows = []
    for name, left, right, registered in CONTRASTS:
        draw = arm_draws[left] - arm_draws[right]
        low, high = np.quantile(draw, [0.025, 0.975])
        delta_rows.append(
            {
                "block": "C02",
                "model": "D-MPNN",
                "comparison": name,
                "left": left,
                "right": right,
                "weighting": "source_equal_seed_mean_abs_error",
                "estimate": float(point[left] - point[right]),
                "ci_low": float(low),
                "ci_high": float(high),
                "n_studies": int(len(study_component)),
                "n_components": n_comp,
                "max_component_studies": largest,
                "interval_excludes_zero": bool(low > 0 or high < 0),
                "in_c3_registry": registered,
                "family": "FTRANSFER" if name == "TL01_FN_minus_JN" else "DESCRIPTIVE",
            }
        )
    deltas = pd.DataFrame(delta_rows)
    records, audit = training_audit()
    compression, regional, corr = diagnostics(seed_mean)
    nearest = nearest_tanimoto(pred)
    published = pd.read_csv(RUN / "c02" / "metrics_pooled.csv")
    published_delta = pooled.merge(published[["arm", "mae_source"]], on="arm", suffixes=("_recomputed", "_published"))
    published_delta["mae_source_difference"] = published_delta.mae_source_recomputed - published_delta.mae_source_published
    pooled.to_csv(RUN / "c02" / "metrics_recomputed.csv", index=False)
    deltas.to_csv(RUN / "c02" / "paired_deltas_recomputed.csv", index=False)
    audit.to_csv(RUN / "c02" / "training_audit.csv", index=False)
    compression.to_csv(RUN / "c02" / "prediction_diagnostics.csv", index=False)
    regional.to_csv(RUN / "c02" / "regional_signed_error.csv", index=False)
    corr.to_csv(RUN / "c02" / "error_correlation.csv", index=False)
    nearest.to_csv(RUN / "c02" / "nearest_tanimoto.csv", index=False)
    published_delta.to_csv(RUN / "c02" / "published_comparison.csv", index=False)
    write_report(coverage_info, pooled, deltas, audit, compression, published_delta[["arm", "mae_source_recomputed", "mae_source_published", "mae_source_difference"]])
    print("C02_SCORE_DONE", len(pred), "max_published_abs_diff", float(published_delta.mae_source_difference.abs().max()))


if __name__ == "__main__":
    main()
