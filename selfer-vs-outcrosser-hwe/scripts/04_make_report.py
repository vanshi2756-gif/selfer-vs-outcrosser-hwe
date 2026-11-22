#!/usr/bin/env python3
"""
04_make_report.py - fill the numeric parts of the report from the real results.

Reads the two *.summary.json files written by 02_run_hwe.py and writes docs/REPORT.md
containing: data and filter details, a results table, the figure, and a one-paragraph
statement of what the numbers are. Everything numeric comes from your results files.

It deliberately does NOT write the Discussion for you. Explaining why the numbers look
the way they do (structure, ascertainment bias, genotyping error) has to be your own
reading of your own results. The Discussion section contains prompts to answer.

Example
-------
  python scripts/04_make_report.py \
      --selfer results/athaliana_hwe.summary.json --outcrosser results/populus_hwe.summary.json \
      --figure ../figures/figure1_selfer_vs_outcrosser.png --out docs/REPORT.md
"""

import argparse
import json
import os


def pct(x: float) -> str:
    return f"{100 * x:.1f}%"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selfer", required=True, help="summary.json for the selfing species")
    ap.add_argument("--outcrosser", required=True, help="summary.json for the outcrossing species")
    ap.add_argument("--figure", default="../figures/figure1_selfer_vs_outcrosser.png",
                    help="figure path as it should appear inside the report (relative to the report file)")
    ap.add_argument("--out", default="docs/REPORT.md")
    ap.add_argument("--synthetic", action="store_true", help="add a banner saying the data are simulated")
    args = ap.parse_args()

    s = json.load(open(args.selfer))
    o = json.load(open(args.outcrosser))

    banner = ""
    if args.synthetic:
        banner = ("> **SYNTHETIC TEST DATA.** These numbers come from simulated genotypes used to check the "
                  "pipeline. They are not a biological result.\n\n")

    def row(label, key, fmt):
        return f"| {label} | {fmt(o[key])} | {fmt(s[key])} |"

    table = "\n".join([
        f"| | {o['label']} | {s['label']} |",
        "|---|---|---|",
        row("Individuals analysed", "individuals_analysed", str),
        row("SNPs tested", "snps_tested", str),
        row("Mean observed heterozygosity", "mean_h_obs", lambda v: f"{v:.3f}"),
        row("Mean expected heterozygosity (2pq)", "mean_h_exp", lambda v: f"{v:.3f}"),
        row("Overall F_IS", "f_is_overall", lambda v: f"{v:.3f}"),
        row("SNPs failing HWE, exact test p < 0.05 (raw)", "frac_sig_exact_p05_raw", pct),
        row("SNPs failing HWE after FDR 5%", "frac_sig_exact_fdr05", pct),
        row("SNPs failing HWE after Bonferroni", "frac_sig_exact_bonferroni", pct),
        row("SNPs where chi-square is valid (all expected counts >= 5)", "frac_chi2_valid_expected_ge5", pct),
    ])

    direction = ("higher in the selfing species, as predicted"
                 if s["f_is_overall"] > o["f_is_overall"] else
                 "NOT higher in the selfing species, which is the opposite of the prediction")

    text = f"""# Mating system and Hardy-Weinberg equilibrium: a selfer vs. an outcrosser

*[Your name, college, date]*

{banner}## 1. Question and prediction
How does mating system change genotype frequencies? Selfing should leave almost no heterozygotes
(observed heterozygosity far below 2pq, F_IS near 1), whereas random mating gives genotype
frequencies p^2, 2pq, q^2 and F_IS near 0. [Add 2-3 lines in your own words on Hardy-Weinberg and F_IS.]

## 2. Data and filters
* Selfer: {s['label']}, from [source, accession or DOI]. File: `{s['vcf']}`.
* Outcrosser: {o['label']}, from [source, accession or DOI]. File: `{o['vcf']}`.
* Individuals analysed: {s['individuals_analysed']} (selfer) and {o['individuals_analysed']} (outcrosser).
  [State which individuals were chosen and why, for example one region or one genetic group.]
* Biallelic SNPs with minor allele frequency >= {s['maf_min']} ; {s['snps_passing_filters']} (selfer) and
  {o['snps_passing_filters']} (outcrosser) passed the filters out of {s['snps_read']} and {o['snps_read']} read;
  {s['snps_tested']} and {o['snps_tested']} were tested after any subsampling.
* [State the marker types, e.g. SNP array vs. whole-genome sequencing, and what that implies.]

## 3. Methods
Per-SNP chi-square test of Hardy-Weinberg proportions (1 d.f.), exact test (Wigginton et al. 2005),
F_IS = 1 - H_obs/H_exp, and Benjamini-Hochberg FDR and Bonferroni correction across SNPs.
Code and unit tests: [repository link].

## 4. Results
![Selfer vs outcrosser]({args.figure})

{table}

Overall F_IS was {s['f_is_overall']:.3f} in the selfer and {o['f_is_overall']:.3f} in the outcrosser, {direction}.
After FDR correction, {pct(s['frac_sig_exact_fdr05'])} of SNPs deviated from Hardy-Weinberg in the selfer and
{pct(o['frac_sig_exact_fdr05'])} in the outcrosser.

## 5. Discussion
Answer these using your own results:
* Do the results match the prediction? By how much?
* What else could produce the heterozygote deficit or excess you see: population structure
  (Wahlund effect), SNP-array ascertainment, genotyping error, selection, linkage?
* How did you handle structure, and what would you do with more time?

## 6. Limitations
Reanalysis of public data, not new experimental work; population structure (Wahlund effect); linkage disequilibrium across neighbouring SNPs; genotyping error; different marker types in the two species (SNP array vs. WGS); [anything specific you found].

## 7. References
McKown et al. (2022) Zenodo 10.5281/zenodo.7129889 (check the record for the associated paper);
1001 Genomes Consortium (2016) Cell 166:481-491; Wigginton et al. (2005) Am J Hum Genet 76:887-893.
"""
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(text)
    print(f"report -> {args.out}")


if __name__ == "__main__":
    main()
