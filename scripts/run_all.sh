#!/usr/bin/env bash
# run_all.sh - run the whole real-data analysis in one go (on your own computer).
#
# Usage:
#   bash scripts/run_all.sh data/athaliana_subset.vcf.gz data/Panel4SNPv3.vcf [N_SNPS]
#         ^ selfer VCF                                     ^ outcrosser VCF      ^ default 3000
#
# Produces: results/*.tsv, results/*.summary.json, results/comparison_summary.csv,
#           figures/figure1_selfer_vs_outcrosser.png, docs/REPORT.md
set -euo pipefail
cd "$(dirname "$0")/.."

SELFER_VCF="${1:?give the A. thaliana VCF path}"
OUTC_VCF="${2:?give the Populus VCF path}"
N="${3:-3000}"
for f in "$SELFER_VCF" "$OUTC_VCF"; do [ -f "$f" ] || { echo "missing file: $f" >&2; exit 1; }; done

python3 tests/test_hwe.py            # stop early if the code itself is broken

python3 scripts/02_run_hwe.py --vcf "$OUTC_VCF"   --label "P. trichocarpa (outcrosser)" \
        --out results/populus_hwe.tsv   --maf 0.05 --max-missing 0.1 --subsample "$N" --seed 1
python3 scripts/02_run_hwe.py --vcf "$SELFER_VCF" --label "A. thaliana (selfer)" \
        --out results/athaliana_hwe.tsv --maf 0.05 --max-missing 0.1 --subsample "$N" --seed 1

python3 scripts/03_compare_plots.py \
        --selfer results/athaliana_hwe.tsv   --selfer-label "A. thaliana (selfer)" \
        --outcrosser results/populus_hwe.tsv --outcrosser-label "P. trichocarpa (outcrosser)"

python3 scripts/04_make_report.py \
        --selfer results/athaliana_hwe.summary.json --outcrosser results/populus_hwe.summary.json

echo
echo "Done. Open docs/REPORT.md and fill the [bracketed] parts and the Discussion."
echo "Check results/*.summary.json: if the selfer shows F_IS far from 1 or the outcrosser far from 0,"
echo "look at population structure before interpreting (see README, caveat 1)."
