#!/usr/bin/env python3
"""
simulate_data.py - SYNTHETIC genotypes for testing the pipeline. NOT real data.

Because the true inbreeding coefficient F is set by us, we can check that the
pipeline recovers it (F_IS ~ F). Never present these results as biology.

Model
-----
  * M biallelic SNPs, N diploid individuals.
  * Each individual is a homozygote "by descent" with probability F (selfing-like),
    otherwise its two alleles are drawn independently (random mating).
  * Optional population structure: --pops K subpopulations with differentiation
    --fst (Balding-Nichols). Pooling them without stratifying creates a spurious
    heterozygote deficit (the Wahlund effect).

Example
-------
  python scripts/simulate_data.py --out data/sim_selfer.vcf --f 0.95 --n 200 --snps 3000
"""

import argparse

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=200, help="individuals")
    ap.add_argument("--snps", type=int, default=3000)
    ap.add_argument("--f", type=float, default=0.0, help="inbreeding coefficient F (0 = random mating)")
    ap.add_argument("--pops", type=int, default=1, help="number of subpopulations")
    ap.add_argument("--fst", type=float, default=0.0, help="differentiation between subpopulations")
    ap.add_argument("--missing", type=float, default=0.02, help="fraction of missing genotypes")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    n, m, F = args.n, args.snps, args.f
    pop_of = np.arange(n) % args.pops

    p0 = np.clip(rng.beta(0.6, 0.6, size=m), 0.06, 0.94)          # ancestral ALT frequency
    if args.pops > 1 and args.fst > 0:
        a = p0 * (1 - args.fst) / args.fst
        b = (1 - p0) * (1 - args.fst) / args.fst
        p_pop = np.stack([rng.beta(a, b) for _ in range(args.pops)])   # pops x snps
    else:
        p_pop = np.tile(p0, (args.pops, 1))
    p = p_pop[pop_of]                                                # n x snps

    a1 = rng.random((n, m)) < p                                      # allele 1 is ALT?
    a2 = rng.random((n, m)) < p
    ibd = rng.random((n, m)) < F                                     # homozygous by descent
    a2 = np.where(ibd, a1, a2)
    geno = a1.astype(int) + a2.astype(int)                           # number of ALT alleles
    miss = rng.random((n, m)) < args.missing

    gt_str = np.array(["0/0", "0/1", "1/1"])
    with open(args.out, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=simulate_data.py SYNTHETIC TEST DATA - not real\n")
        fh.write(f"##simulation=F={F};pops={args.pops};fst={args.fst};seed={args.seed}\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + "\t".join(f"S{i+1:03d}" for i in range(n)) + "\n")
        for j in range(m):
            col = np.where(miss[:, j], "./.", gt_str[geno[:, j]])
            fh.write(f"chr1\t{j+1}\tsnp{j+1}\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(col) + "\n")
    print(f"wrote {args.out}: {n} individuals x {m} SNPs, F={F}, pops={args.pops}, fst={args.fst}")


if __name__ == "__main__":
    main()
