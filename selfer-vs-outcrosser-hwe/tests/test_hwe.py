"""
Tests for hwe_lib. Run with:  python tests/test_hwe.py   (or: pytest tests/)
"""
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import hwe_lib as hl  # noqa: E402


def test_exact_matches_bruteforce():
    cases = [(10, 40, 50), (0, 30, 70), (50, 25, 25), (3, 1, 96), (1, 0, 99), (40, 30, 30), (200, 150, 150)]
    for het, h1, h2 in cases:
        a, b = hl.hwe_exact(het, h1, h2), hl.hwe_exact_bruteforce(het, h1, h2)
        assert abs(a - b) < 1e-9, (het, h1, h2, a, b)


def test_perfect_hwe_gives_p_one():
    # p = 0.5, n = 400 -> expected 100 / 200 / 100
    st = hl.site_stats(100, 200, 100)
    assert abs(st["chi2"]) < 1e-12 and abs(st["p_chi2"] - 1.0) < 1e-12
    assert abs(st["f_is"]) < 1e-12 and abs(st["p_ref"] - 0.5) < 1e-12


def test_chi2_pvalue_matches_scipy():
    from scipy.stats import chi2
    st = hl.site_stats(60, 20, 20)
    assert abs(st["p_chi2"] - chi2.sf(st["chi2"], 1)) < 1e-12


def test_textbook_example():
    # Worked by hand: 50 AA, 20 Aa, 30 aa  -> p = 0.6, q = 0.4, n = 100
    st = hl.site_stats(50, 20, 30)
    assert abs(st["p_ref"] - 0.6) < 1e-12
    assert abs(st["exp_homref"] - 36) < 1e-9 and abs(st["exp_het"] - 48) < 1e-9 and abs(st["exp_homalt"] - 16) < 1e-9
    assert abs(st["chi2"] - (14**2 / 36 + 28**2 / 48 + 14**2 / 16)) < 1e-9
    assert abs(st["f_is"] - (1 - 0.20 / 0.48)) < 1e-9


def test_benjamini_hochberg():
    q = hl.benjamini_hochberg([0.01, 0.04, 0.03, 0.005])
    assert np.allclose(sorted(q), sorted([0.02, 0.04, 0.04, 0.02]))


def test_pipeline_recovers_simulated_F():
    with tempfile.TemporaryDirectory() as tmp:
        for f_true, lo, hi in [(0.0, -0.05, 0.05), (0.9, 0.84, 0.96)]:
            vcf = os.path.join(tmp, f"sim_{f_true}.vcf")
            out = os.path.join(tmp, f"res_{f_true}.tsv")
            subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "simulate_data.py"), "--out", vcf,
                            "--f", str(f_true), "--n", "200", "--snps", "1500"], check=True, capture_output=True)
            subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "02_run_hwe.py"), "--vcf", vcf,
                            "--label", "t", "--out", out], check=True, capture_output=True)
            import json
            s = json.load(open(os.path.splitext(out)[0] + ".summary.json"))
            assert lo < s["f_is_overall"] < hi, (f_true, s["f_is_overall"])


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
