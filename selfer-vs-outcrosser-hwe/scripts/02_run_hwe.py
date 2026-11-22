#!/usr/bin/env python3
"""
02_run_hwe.py - per-SNP Hardy-Weinberg tests and F_IS for one VCF.

Run it once per species, then compare the two outputs with 03_compare_plots.py.

Example
-------
  python scripts/02_run_hwe.py --vcf data/Panel4SNPv3.vcf --label "P. trichocarpa (outcrosser)" \
         --out results/populus_hwe.tsv --maf 0.05 --max-missing 0.1 --subsample 5000 --seed 1
"""

import argparse
import json
import os
import random
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hwe_lib as hl  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="input VCF (.vcf or .vcf.gz)")
    ap.add_argument("--label", required=True, help="dataset name used in outputs and plots")
    ap.add_argument("--out", required=True, help="output TSV path (a .summary.json is written next to it)")
    ap.add_argument("--samples", help="text file with one sample ID per line; analyse only these samples")
    ap.add_argument("--maf", type=float, default=0.05, help="minimum minor allele frequency (default 0.05)")
    ap.add_argument("--max-missing", type=float, default=0.10,
                    help="maximum fraction of missing genotypes per SNP (default 0.10)")
    ap.add_argument("--min-called", type=int, default=20, help="minimum genotyped individuals per SNP (default 20)")
    ap.add_argument("--thin", type=int, default=1,
                    help="keep every Nth SNP that passes filters (cheap way to reduce linkage; default 1 = keep all)")
    ap.add_argument("--subsample", type=int, default=0,
                    help="randomly keep this many SNPs after filtering (use the same number for both species)")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    keep = None
    if args.samples:
        with open(args.samples) as fh:
            keep = [ln.strip() for ln in fh if ln.strip()]

    n_seen = n_pass = 0
    rows = []
    n_samples = 0
    for rec in hl.iter_vcf_counts(args.vcf, keep_samples=keep):
        n_seen += 1
        n_samples = rec["n_samples"]
        n_called = rec["n0"] + rec["n1"] + rec["n2"]
        if n_called < args.min_called:
            continue
        if rec["n_missing"] / rec["n_samples"] > args.max_missing:
            continue
        st = hl.site_stats(rec["n0"], rec["n1"], rec["n2"])
        if st is None or st["maf"] < args.maf:
            continue
        n_pass += 1
        if args.thin > 1 and (n_pass % args.thin) != 0:
            continue
        rows.append({**{k: rec[k] for k in ("chrom", "pos", "id", "ref", "alt", "n0", "n1", "n2")}, **st})

    if args.subsample and len(rows) > args.subsample:
        random.Random(args.seed).shuffle(rows)
        rows = rows[: args.subsample]
        def _chrom_key(r: dict):
            c = re.sub(r"^chr", "", str(r["chrom"]), flags=re.IGNORECASE)
            return (0, int(c), r["pos"]) if c.isdigit() else (1, str(r["chrom"]), r["pos"])

        rows.sort(key=_chrom_key)
    if not rows:
        sys.exit("No SNPs passed the filters - check the VCF, sample list and thresholds.")

    # Multiple-testing correction across all SNPs tested in this dataset
    m = len(rows)
    q_exact = hl.benjamini_hochberg([r["p_exact"] for r in rows])
    q_chi2 = hl.benjamini_hochberg([r["p_chi2"] for r in rows])
    for r, qe, qc in zip(rows, q_exact, q_chi2):
        r["q_bh_exact"] = float(qe)
        r["q_bh_chi2"] = float(qc)
        r["bonferroni_sig_exact"] = bool(r["p_exact"] < 0.05 / m)
        r["chi2_valid"] = bool(r["min_expected"] >= 5)

    cols = ["chrom", "pos", "id", "ref", "alt", "n0", "n1", "n2", "n_called", "p_ref", "maf",
            "h_obs", "h_exp", "f_is", "exp_homref", "exp_het", "exp_homalt", "chi2", "p_chi2",
            "p_exact", "min_expected", "chi2_valid", "q_bh_chi2", "q_bh_exact", "bonferroni_sig_exact"]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

    h_obs = np.array([r["h_obs"] for r in rows])
    h_exp = np.array([r["h_exp"] for r in rows])
    p_chi = np.array([r["p_chi2"] for r in rows])
    p_ex = np.array([r["p_exact"] for r in rows])
    summary = {
        "label": args.label,
        "vcf": os.path.basename(args.vcf),
        "individuals_analysed": n_samples,
        "snps_read": n_seen,
        "snps_passing_filters": n_pass,
        "snps_tested": m,
        "maf_min": args.maf,
        "mean_h_obs": float(h_obs.mean()),
        "mean_h_exp": float(h_exp.mean()),
        "f_is_overall": float(1 - h_obs.sum() / h_exp.sum()),
        "f_is_mean_per_snp": float(np.mean([r["f_is"] for r in rows])),
        "frac_sig_chi2_p05_raw": float((p_chi < 0.05).mean()),
        "frac_sig_exact_p05_raw": float((p_ex < 0.05).mean()),
        "frac_sig_exact_fdr05": float((np.array([r["q_bh_exact"] for r in rows]) < 0.05).mean()),
        "frac_sig_exact_bonferroni": float(np.mean([r["bonferroni_sig_exact"] for r in rows])),
        "frac_chi2_valid_expected_ge5": float(np.mean([r["chi2_valid"] for r in rows])),
    }
    with open(os.path.splitext(args.out)[0] + ".summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"{args.label}: {m} SNPs tested in {n_samples} individuals")
    print(f"  mean H_obs = {summary['mean_h_obs']:.3f}   mean H_exp = {summary['mean_h_exp']:.3f}   "
          f"overall F_IS = {summary['f_is_overall']:.3f}")
    print(f"  SNPs failing HWE (exact, raw p<0.05): {100*summary['frac_sig_exact_p05_raw']:.1f}%  "
          f"(5% expected by chance if HWE held)")
    print(f"  after FDR 5%: {100*summary['frac_sig_exact_fdr05']:.1f}%   after Bonferroni: "
          f"{100*summary['frac_sig_exact_bonferroni']:.1f}%")
    print(f"  results -> {args.out}")


if __name__ == "__main__":
    main()
