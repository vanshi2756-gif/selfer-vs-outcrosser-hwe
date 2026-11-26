#!/usr/bin/env bash
# 01_fetch_data.sh - download the two public datasets. Run this on YOUR computer.
# (Raw VCF files are not stored in the repo. The analysis that was run on them is in
#  results/ and docs/REPORT.md.)
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data

# ---------------------------------------------------------------------------
# 1) OUTCROSSER: Populus trichocarpa (black cottonwood), SNP-array genotypes
#    Zenodo record 7129889, DOI 10.5281/zenodo.7129889, licence CC0.
#    Panel4SNPv3.vcf is ~46 MB; the record lists its md5 as below.
# ---------------------------------------------------------------------------
curl -L -o data/Panel4SNPv3.vcf "https://zenodo.org/records/7129889/files/Panel4SNPv3.vcf?download=1"
EXPECTED_MD5="3c26eba2aa931c716ff9dd41e0d86907"
if command -v md5sum >/dev/null 2>&1; then
    echo "${EXPECTED_MD5}  data/Panel4SNPv3.vcf" | md5sum -c -
elif command -v md5 >/dev/null 2>&1; then
    ACTUAL_MD5="$(md5 -r data/Panel4SNPv3.vcf | awk '{print $1}')"
    [ "${ACTUAL_MD5}" = "${EXPECTED_MD5}" ] || { echo "MD5 mismatch for data/Panel4SNPv3.vcf" >&2; exit 1; }
    echo "data/Panel4SNPv3.vcf: OK"
fi

# Sample metadata (check these for sampling locations so you can stratify by region):
curl -L -o data/NEW_IDs.xlsx          "https://zenodo.org/records/7129889/files/NEW_IDs.xlsx?download=1"
curl -L -o data/NEW_408_accession.xlsx "https://zenodo.org/records/7129889/files/NEW_408_accession.xlsx?download=1"

# ---------------------------------------------------------------------------
# 2) SELFER: Arabidopsis thaliana, 1001 Genomes (1,135 accessions) - do this by hand
#    The full VCF is huge, so use the project's subset tool instead:
#        https://tools.1001genomes.org/vcfsubset
#    Its four steps are: select strains -> select loci -> select format -> download.
#      * Strains: pick ONE genetically coherent group (e.g. accessions from a single
#        country or admixture group - see https://tools.1001genomes.org for the
#        ADMIXTURE group viewer). About 100-300 accessions is plenty.
#      * Loci: keep the region small, e.g. Chr1:1..2000000. The tool writes a record for
#        every base position, so large regions give multi-GB files and the server can
#        fail part-way (the file then ends with an HTML error page). Check that the
#        last line of the file is a complete VCF record before using it.
#      * Format: VCF.
#    Save the result as data/athaliana_subset.vcf(.gz).
#    The reported run used the 110 accessions from the "Spain" filter, Chr1:1..2,000,000.
# ---------------------------------------------------------------------------
echo
echo "Populus data downloaded. Now fetch the A. thaliana subset by hand (see comments above)."
