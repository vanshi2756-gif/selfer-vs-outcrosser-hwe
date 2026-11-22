# Mating system and Hardy-Weinberg equilibrium: a selfer vs. an outcrosser

*Author, college, date. Fill in after running the real analysis. Replace every [bracket].*

## 1. Question and prediction (4-5 lines)
How does mating system change genotype frequencies? Predict heterozygote deficit and F_IS near 1 in
*A. thaliana* (selfing) and F_IS near 0 in *P. trichocarpa* (outcrossing). Briefly explain Hardy-Weinberg
expectations (p², 2pq, q²) and F_IS.

## 2. Data
* Species, source and accession (DOI or URL), number of individuals, marker type.
* Filters used (MAF ≥ [0.05], missing ≤ [10%], [N] SNPs after subsampling, seed [1]).
* Which individuals were analysed and why (one coherent group / region).

## 3. Methods
Per-SNP chi-square test (1 d.f.), exact test (Wigginton et al. 2005), F_IS, Benjamini-Hochberg FDR and
Bonferroni. Code: [repo link]. Mention the unit tests.

## 4. Results
Insert `figures/figure1_selfer_vs_outcrosser.png` and `results/comparison_summary.csv`.
Report for each species: mean H_obs and H_exp, overall F_IS, and % of SNPs failing HWE (raw, FDR, Bonferroni).

## 5. Discussion
* Do the results match the prediction? By how much?
* Why a deviation cannot be attributed to selfing alone (structure, ascertainment, genotyping error, selection).
* What the Wahlund effect would do here and how you handled it.

## 6. Limitations (be explicit)
Public data reanalysis, not new experimental work; population structure (Wahlund effect); linkage disequilibrium across neighbouring SNPs; genotyping error; different marker types (SNP array vs. WGS); [anything specific you found].

## 7. References
McKown et al. (2022) Zenodo 10.5281/zenodo.7129889; 1001 Genomes Consortium (2016) Cell 166:481;
Wigginton et al. (2005) Am J Hum Genet 76:887.
