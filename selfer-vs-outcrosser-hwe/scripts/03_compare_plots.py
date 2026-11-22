#!/usr/bin/env python3
"""
03_compare_plots.py - compare a selfing and an outcrossing dataset.

Reads the two TSV files from 02_run_hwe.py and writes:
  figures/figure1_selfer_vs_outcrosser.png   (4 panels)
  results/comparison_summary.csv             (one row per dataset)

Example
-------
  python scripts/03_compare_plots.py \
      --selfer results/athaliana_hwe.tsv --selfer-label "A. thaliana (selfer)" \
      --outcrosser results/populus_hwe.tsv --outcrosser-label "P. trichocarpa (outcrosser)"
"""

import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# colour-blind-safe pair (Okabe-Ito)
COL = {"outcrosser": "#0072B2", "selfer": "#D55E00"}


def qq(ax, p, color, label):
    p = np.sort(np.clip(np.asarray(p, float), 1e-300, 1.0))
    m = len(p)
    expected = (np.arange(1, m + 1) - 0.5) / m
    step = max(1, m // 3000)          # thin the points so the plot stays light
    ax.scatter(-np.log10(expected[::step]), -np.log10(p[::step]), s=6, color=color, alpha=0.7,
               label=label, linewidths=0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selfer", required=True)
    ap.add_argument("--outcrosser", required=True)
    ap.add_argument("--selfer-label", default="Selfer")
    ap.add_argument("--outcrosser-label", default="Outcrosser")
    ap.add_argument("--figure", default="figures/figure1_selfer_vs_outcrosser.png")
    ap.add_argument("--table", default="results/comparison_summary.csv")
    ap.add_argument("--title", default="", help="optional figure title (e.g. 'SYNTHETIC TEST DATA')")
    args = ap.parse_args()

    d = {"selfer": pd.read_csv(args.selfer, sep="\t"), "outcrosser": pd.read_csv(args.outcrosser, sep="\t")}
    lab = {"selfer": args.selfer_label, "outcrosser": args.outcrosser_label}
    order = ["outcrosser", "selfer"]

    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
    fig.subplots_adjust(hspace=0.33, wspace=0.27)

    # A: F_IS distribution
    ax = axs[0, 0]
    bins = np.linspace(-1, 1, 41)
    for k in order:
        ax.hist(np.clip(d[k]["f_is"], -1, 1), bins=bins, density=True, alpha=0.65, color=COL[k], label=lab[k])
        ax.axvline(1 - d[k]["h_obs"].sum() / d[k]["h_exp"].sum(), color=COL[k], ls="--", lw=1.5)
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_xlabel("F_IS per SNP  (dashed line = overall F_IS)")
    ax.set_ylabel("density")
    ax.set_title("A. Inbreeding coefficient", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)

    # B: observed vs expected heterozygosity
    ax = axs[0, 1]
    for k in order:
        s = d[k].sample(n=min(4000, len(d[k])), random_state=1)
        ax.scatter(s["h_exp"], s["h_obs"], s=5, alpha=0.35, color=COL[k], label=lab[k], linewidths=0)
    ax.plot([0, 0.5], [0, 0.5], color="black", lw=1)
    ax.set_xlim(0, 0.52)
    ax.set_ylim(0, 0.52)
    ax.set_xlabel("expected heterozygosity (2pq)")
    ax.set_ylabel("observed heterozygosity")
    ax.set_title("B. Observed vs expected heterozygosity", loc="left", fontsize=11, fontweight="bold")
    ax.text(0.22, 0.095, "points on the line = Hardy-Weinberg", fontsize=8)
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    # C: QQ plot of exact-test p-values
    ax = axs[1, 0]
    for k in order:
        qq(ax, d[k]["p_exact"], COL[k], lab[k])
    x_top = max(-np.log10(0.5 / len(d[name])) for name in order) * 1.05   # expected values only go this far
    ax.plot([0, x_top], [0, x_top], color="black", lw=1)
    ax.set_xlim(-0.1, x_top)
    ax.set_xlabel("expected  -log10(p)  if HWE holds")
    ax.set_ylabel("observed  -log10(p)  (exact test)")
    ax.set_title("C. How far from HWE? (QQ plot)", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)

    # D: minor allele frequency spectrum
    ax = axs[1, 1]
    bins = np.linspace(0, 0.5, 26)
    for k in order:
        ax.hist(d[k]["maf"], bins=bins, density=True, alpha=0.65, color=COL[k], label=lab[k])
    ax.set_xlabel("minor allele frequency")
    ax.set_ylabel("density")
    ax.set_title("D. Allele-frequency spectrum", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)

    for a in axs.ravel():
        a.spines[["top", "right"]].set_visible(False)
    if args.title:
        fig.suptitle(args.title, fontsize=12, color="#9a1c1c", fontweight="bold")
    os.makedirs(os.path.dirname(os.path.abspath(args.figure)), exist_ok=True)
    fig.savefig(args.figure, dpi=200, bbox_inches="tight")

    # Summary table
    rows = []
    for k in order:
        t = d[k]
        m = len(t)
        rows.append({
            "dataset": lab[k], "snps_tested": m,
            "mean_H_obs": round(t["h_obs"].mean(), 4), "mean_H_exp": round(t["h_exp"].mean(), 4),
            "overall_F_IS": round(1 - t["h_obs"].sum() / t["h_exp"].sum(), 4),
            "pct_fail_HWE_raw_p05_exact": round(100 * (t["p_exact"] < 0.05).mean(), 2),
            "pct_fail_HWE_FDR05": round(100 * (t["q_bh_exact"] < 0.05).mean(), 2),
            "pct_fail_HWE_bonferroni": round(100 * t["bonferroni_sig_exact"].mean(), 2),
            "pct_chi2_valid_expected_ge5": round(100 * t["chi2_valid"].mean(), 2),
        })
    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(os.path.abspath(args.table)), exist_ok=True)
    out.to_csv(args.table, index=False)
    print(out.to_string(index=False))
    print(f"\nfigure -> {args.figure}\ntable  -> {args.table}")


if __name__ == "__main__":
    main()
