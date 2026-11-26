"""
hwe_lib.py - core functions for Hardy-Weinberg and F_IS analysis of biallelic SNPs.

Everything here is plain Python + NumPy so it is easy to read and check against
textbook problems (Genetics syllabus, Unit 7: allele/genotype frequencies,
Hardy-Weinberg law, chi-square test).

Notation used throughout
------------------------
  n0, n1, n2 : observed counts of genotypes REF/REF, REF/ALT, ALT/ALT
  p, q       : allele frequencies of REF and ALT (p + q = 1)
  Expected genotype counts under Hardy-Weinberg equilibrium (HWE):
                 REF/REF = n * p^2,   REF/ALT = n * 2pq,   ALT/ALT = n * q^2
  H_obs      : observed heterozygosity  = n1 / n
  H_exp      : expected heterozygosity  = 2pq
  F_IS       : inbreeding coefficient   = 1 - H_obs / H_exp
               ~0 for random mating, -> 1 for strong inbreeding (e.g. selfing),
               negative when there is an excess of heterozygotes.
"""

from __future__ import annotations

import gzip
import math
from typing import Iterator, Optional, Sequence

import numpy as np

# Genotype strings -> number of ALT alleles. Anything else is treated as missing.
GT_MAP = {
    "0/0": 0, "0|0": 0,
    "0/1": 1, "1/0": 1, "0|1": 1, "1|0": 1,
    "1/1": 2, "1|1": 2,
}


# --------------------------------------------------------------------------- #
# VCF reading
# --------------------------------------------------------------------------- #
def open_text(path: str):
    """Open a plain or gzip/bgzip-compressed text file (detected by magic bytes)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_sample_ids(path: str) -> list[str]:
    """Return the sample IDs from the #CHROM header line of a VCF."""
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise ValueError(f"No #CHROM header line found in {path}")


def iter_vcf_counts(path: str, keep_samples: Optional[Sequence[str]] = None) -> Iterator[dict]:
    """
    Stream a VCF and yield genotype counts for every biallelic SNP.

    Only records with a single-base REF and a single-base ALT are used
    (indels and multi-allelic sites are skipped). Genotypes are read from the GT
    field. Haploid or partly missing calls count as missing.

    keep_samples : optional list of sample IDs. Use it to restrict the analysis
                   to one genetically coherent group (see README: Wahlund effect).
    """
    samples = None
    idx: list[int] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                if keep_samples is not None:
                    keep = set(keep_samples)
                    idx = [i for i, s in enumerate(samples) if s in keep]
                    if not idx:
                        raise ValueError("None of the requested samples were found in the VCF")
                else:
                    idx = list(range(len(samples)))
                continue
            if samples is None:
                raise ValueError("Data line found before the #CHROM header")

            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            ref, alt = f[3], f[4]
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                continue  # not a biallelic SNP
            fmt = f[8].split(":")
            if "GT" not in fmt:
                continue
            g = fmt.index("GT")
            cols = f[9:]

            n0 = n1 = n2 = nm = 0
            for i in idx:
                c = cols[i]
                gt = c[:3] if g == 0 else c.split(":")[g]
                v = GT_MAP.get(gt, -1)
                if v == 0:
                    n0 += 1
                elif v == 1:
                    n1 += 1
                elif v == 2:
                    n2 += 1
                else:
                    nm += 1
            yield {
                "chrom": f[0], "pos": int(f[1]), "id": f[2], "ref": ref, "alt": alt,
                "n0": n0, "n1": n1, "n2": n2, "n_missing": nm, "n_samples": len(idx),
            }


# --------------------------------------------------------------------------- #
# Hardy-Weinberg tests
# --------------------------------------------------------------------------- #
def hwe_exact(n_het: int, n_hom1: int, n_hom2: int) -> float:
    """
    Exact test of Hardy-Weinberg equilibrium (two-sided), following the
    algorithm of Wigginton, Cutler & Abecasis (2005), Am J Hum Genet 76:887-893.

    Unlike the chi-square test it stays valid when some expected counts are small,
    which is common for rare alleles. Returns the p-value.
    """
    hom_c = max(n_hom1, n_hom2)   # common homozygote
    hom_r = min(n_hom1, n_hom2)   # rare homozygote
    rare = 2 * hom_r + n_het      # copies of the rare allele
    n = n_het + hom_c + hom_r     # individuals
    if n == 0 or rare == 0:
        return 1.0

    probs = [0.0] * (rare + 1)
    mid = rare * (2 * n - rare) // (2 * n)   # most likely number of heterozygotes
    if (mid % 2) != (rare % 2):
        mid += 1

    cur_hets = mid
    cur_homr = (rare - mid) // 2
    cur_homc = n - cur_hets - cur_homr
    probs[mid] = 1.0
    total = 1.0

    # walk down from the mode
    while cur_hets > 1:
        probs[cur_hets - 2] = probs[cur_hets] * cur_hets * (cur_hets - 1.0) / (
            4.0 * (cur_homr + 1.0) * (cur_homc + 1.0))
        total += probs[cur_hets - 2]
        cur_homr += 1
        cur_homc += 1
        cur_hets -= 2

    # walk up from the mode
    cur_hets = mid
    cur_homr = (rare - mid) // 2
    cur_homc = n - cur_hets - cur_homr
    while cur_hets <= rare - 2:
        probs[cur_hets + 2] = probs[cur_hets] * 4.0 * cur_homr * cur_homc / (
            (cur_hets + 2.0) * (cur_hets + 1.0))
        total += probs[cur_hets + 2]
        cur_homr -= 1
        cur_homc -= 1
        cur_hets += 2

    p_obs = probs[n_het]
    p_val = sum(x for x in probs if x <= p_obs) / total
    return min(1.0, p_val)


def hwe_exact_bruteforce(n_het: int, n_hom1: int, n_hom2: int) -> float:
    """
    Slow but transparent version of the exact test, used only to check hwe_exact().
    P(h hets | n individuals, r rare-allele copies) =
        n! r! (2n-r)! 2^h / ( ((r-h)/2)! h! ((2n-r-h)/2)! (2n)! )
    """
    hom_r = min(n_hom1, n_hom2)
    n = n_het + n_hom1 + n_hom2
    r = 2 * hom_r + n_het
    lg = math.lgamma

    def logp(h: int) -> float:
        return (lg(n + 1) + lg(r + 1) + lg(2 * n - r + 1) + h * math.log(2)
                - lg((r - h) / 2 + 1) - lg(h + 1) - lg((2 * n - r - h) / 2 + 1) - lg(2 * n + 1))

    hs = range(r % 2, min(r, 2 * n - r) + 1, 2)
    probs = {h: math.exp(logp(h)) for h in hs}
    p_obs = probs[n_het]
    return min(1.0, sum(v for v in probs.values() if v <= p_obs * (1 + 1e-9)))


def site_stats(n0: int, n1: int, n2: int) -> Optional[dict]:
    """
    Allele frequencies, expected counts, chi-square HWE test (1 d.f.), exact test
    and F_IS for one SNP. Returns None for monomorphic sites.
    """
    n = n0 + n1 + n2
    if n == 0:
        return None
    p = (2 * n0 + n1) / (2 * n)
    q = 1.0 - p
    if p <= 0.0 or q <= 0.0:
        return None

    e0, e1, e2 = n * p * p, n * 2 * p * q, n * q * q
    if e0 <= 0.0 or e1 <= 0.0 or e2 <= 0.0:
        chi2 = float("nan")
        p_chi2 = float("nan")
    else:
        chi2 = (n0 - e0) ** 2 / e0 + (n1 - e1) ** 2 / e1 + (n2 - e2) ** 2 / e2
        p_chi2 = math.erfc(math.sqrt(chi2 / 2.0))   # survival function of chi-square, 1 d.f.
    h_obs = n1 / n
    h_exp = 2 * p * q
    return {
        "n_called": n, "p_ref": p, "maf": min(p, q),
        "h_obs": h_obs, "h_exp": h_exp, "f_is": 1.0 - h_obs / h_exp,
        "exp_homref": e0, "exp_het": e1, "exp_homalt": e2,
        "chi2": chi2, "p_chi2": p_chi2,
        "p_exact": hwe_exact(n1, n0, n2),
        "min_expected": min(e0, e1, e2),
    }


# --------------------------------------------------------------------------- #
# Multiple testing
# --------------------------------------------------------------------------- #
def benjamini_hochberg(pvals: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (false discovery rate 'q-values')."""
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    return out
