#!/usr/bin/env python3
"""Redraw the three NEXTrun figures from sealed artifacts.

Reads artifacts/v1/curves.npz and artifacts/v1/result.json. Does not rewrite
the report or copy artifacts.
"""
import json
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
ART = ROOT / "artifacts" / "v1"
OUT = ROOT / "figures"
CURVES = "artifacts/v1/curves.npz"
RESULT = "artifacts/v1/result.json"

NAVY, ORANGE, GREEN, RED, GRAY = "#1f4e79", "#e07b39", "#3a9a5b", "#b23a48", "#6b7280"
METHOD_COLOR = {"B": NAVY, "F": GRAY, "G": ORANGE, "A": GREEN}
EDGES = (160, 180, 210, 260)
SENSITIVITY = [
    "half_step-B", "half_step-G", "half_step-A",
    "index20-B", "index20-G", "index20-A",
    "A_full-A", "A_clip-A", "A_order-A", "A_inf-A",
    "G_L-B", "rolling_F-F", "B_df-B",
]


def use_cjk():
    from matplotlib import font_manager

    listed = subprocess.check_output(["fc-list", "Noto Sans CJK SC", "file"], text=True)
    paths = [line.split(":")[0].strip() for line in listed.splitlines() if line.strip()]
    if not paths:
        raise SystemExit("Noto Sans CJK SC is not installed")
    for path in paths:
        font_manager.fontManager.addfont(path)
    plt.rcParams.update({
        "font.family": "Noto Sans CJK SC",
        "font.weight": "regular",
        "axes.titleweight": "regular",
        "axes.labelweight": "regular",
        "axes.unicode_minus": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "#333333",
        "text.color": "#222222",
        "axes.labelcolor": "#222222",
        "xtick.color": "#222222",
        "ytick.color": "#222222",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 8,
    })


def strip_png_text(path):
    data = Path(path).read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise SystemExit(f"not a png: {path.name}")
    out = bytearray(data[:8])
    i = 8
    while i + 8 <= len(data):
        length = int.from_bytes(data[i:i + 4], "big")
        ctype = data[i + 4:i + 8]
        end = i + 12 + length
        if end > len(data):
            raise SystemExit(f"truncated png: {path.name}")
        if ctype not in (b"tEXt", b"iTXt", b"zTXt"):
            out += data[i:end]
        i = end
        if ctype == b"IEND":
            break
    Path(path).write_bytes(out)
    width = int.from_bytes(out[16:20], "big")
    if width != 1600:
        raise SystemExit(f"{path.name} width {width} != 1600")


def save_png(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    strip_png_text(path)


def yerr(mean, ci):
    lo, hi = ci
    if not (lo <= mean <= hi):
        raise SystemExit(f"interval does not cover the mean: {lo}, {mean}, {hi}")
    return np.array([[mean - lo], [hi - mean]])


def summary_metric(result, condition, method, window, metric):
    return result["summary"][f"main/{condition}/{method}/{window}"][metric]


def fig_curves(curves):
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.6), sharex=True)
    for ci, (condition, title) in enumerate((
        ("constant", "常噪声 Constant noise"),
        ("shift", "突变噪声 Shift noise"),
    )):
        ax = axes[ci]
        for mi, method in enumerate(("B", "F", "G", "A")):
            line = curves[ci, mi, :, 0]
            ax.plot(np.arange(800), line, color=METHOD_COLOR[method], alpha=0.18, lw=0.6)
            tail = np.convolve(line, np.ones(100) / 100, mode="valid")
            ax.plot(np.arange(99, 800), tail, color=METHOD_COLOR[method], lw=1.4, label=method)
        for edge in EDGES:
            ax.axvline(edge, color="#222222", ls=":", alpha=0.35, lw=0.7)
        ax.axhline(0.95, color=RED, ls="--", lw=0.9, label="名义 0.95 Nominal" if ci == 0 else None)
        ax.set_ylim(0.4, 1.015)
        ax.set_ylabel("覆盖率 Coverage")
        ax.set_title(title + "：淡线为逐步均值，实线为 100 点尾随\nFaint = per-step mean; solid = trailing 100")
        ax.grid(axis="y", color="#e5e7eb")
        ax.set_axisbelow(True)
    axes[0].legend(ncol=5, frameon=False, fontsize=7, loc="lower right")
    axes[1].set_xlabel("从零计的测试索引 Zero-based test index\n竖线 Vertical lines: 160 / 180 / 210 / 260。纵轴自 0.4 起 Axis starts at 0.4.")
    fig.suptitle("逐步覆盖 Stepwise coverage", fontsize=12, y=0.98)
    fig.text(0.01, 0.005, f"数据来源 Source: {CURVES}  (array main, coverage channel)", fontsize=7, color=GRAY)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.16, hspace=0.42)
    save_png(fig, OUT / "coverage_curves.png")


def fig_costs(result):
    fig, axes = plt.subplots(2, 3, figsize=(8.0, 5.4))
    metrics = (
        ("coverage", "覆盖率 Coverage"),
        ("width", "宽度 Width"),
        ("winkler", "Winkler"),
    )
    rows = (
        ("constant", "all", "常噪声 Constant"),
        ("shift", "high", "突变 Shift"),
    )
    for row, (condition, window, row_name) in enumerate(rows):
        for col, (metric, metric_name) in enumerate(metrics):
            ax = axes[row, col]
            for k, method in enumerate(("B", "F", "G", "A")):
                item = summary_metric(result, condition, method, window, metric)
                if not item.get("estimable", True):
                    raise SystemExit(f"not estimable: {condition} {method} {window} {metric}")
                mean = item["mean"]
                ax.errorbar(k, mean, yerr=yerr(mean, item["ci95"]), fmt="o",
                            color=METHOD_COLOR[method], capsize=3, ms=5)
            ax.set_xticks(range(4), ["B", "F", "G", "A"])
            ax.set_title(f"{row_name} · {metric_name}", fontsize=8)
            if metric == "coverage":
                ax.axhline(0.95, color=RED, ls="--", lw=0.8)
            if col == 0:
                ax.set_ylabel("均值与 95% 区间\nMean and 95% CI", fontsize=7)
            ax.grid(axis="y", color="#e5e7eb")
            ax.set_axisbelow(True)
    fig.suptitle("覆盖与代价 Coverage and cost", fontsize=12)
    fig.text(0.01, 0.005, f"数据来源 Source: {RESULT}    误差线为配对 ID bootstrap 95% 百分位区间 "
             "Error bars are paired-ID bootstrap 95% percentile intervals.", fontsize=6.5, color=GRAY)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.86, bottom=0.10, hspace=0.42, wspace=0.38)
    save_png(fig, OUT / "coverage_costs.png")


def fig_sensitivity(result):
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 6.5), sharey=True)
    titles = {"early50": "前 50 点 Early 50", "high": "高噪声全段 High-noise window"}
    for ax, window in zip(axes, ("early50", "high")):
        for k, name in enumerate(SENSITIVITY):
            item = result["comparisons"][f"{name}/shift/{window}"]["coverage"]
            if not item.get("estimable", True):
                raise SystemExit(f"not estimable: {name} {window}")
            mean = 100 * item["mean"]
            lo, hi = (100 * x for x in item["ci95"])
            ax.errorbar(mean, k, xerr=yerr(mean, (lo, hi)), fmt="o", color=NAVY, capsize=2, ms=4)
        ax.axvline(0, color=GRAY, ls="--", lw=0.8)
        ax.set_xlabel("配对覆盖变化（百分点）\nPaired coverage change (pp)")
        ax.set_title(titles[window] + "（描述性 descriptive）", fontsize=9)
        ax.grid(axis="x", color="#e5e7eb")
        ax.set_axisbelow(True)
    axes[0].set_yticks(range(len(SENSITIVITY)), SENSITIVITY, fontsize=7)
    axes[0].invert_yaxis()
    fig.suptitle("预定敏感性 Prespecified sensitivity", fontsize=12)
    fig.text(0.01, 0.005, f"数据来源 Source: {RESULT}    突变条件 shift condition。不同变体含不同机制，不替换主结果。 "
             "Variants are not alternative primary results.", fontsize=6.5, color=GRAY)
    fig.subplots_adjust(left=0.16, right=0.98, top=0.88, bottom=0.14, wspace=0.12)
    save_png(fig, OUT / "sensitivity.png")


def main():
    use_cjk()
    curves_path = ART / "curves.npz"
    result_path = ART / "result.json"
    if not curves_path.is_file() or not result_path.is_file():
        raise SystemExit(f"missing {CURVES} or {RESULT}")
    curves = np.load(curves_path)["main"]
    if curves.shape[0] < 2 or curves.shape[1] < 4 or curves.shape[2] != 800:
        raise SystemExit(f"unexpected curves shape {curves.shape}")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    fig_curves(curves)
    fig_costs(result)
    fig_sensitivity(result)
    print("wrote figures/coverage_curves.png")
    print("wrote figures/coverage_costs.png")
    print("wrote figures/sensitivity.png")


if __name__ == "__main__":
    main()
