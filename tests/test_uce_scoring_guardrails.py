from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from types import SimpleNamespace

from main_assembler import (
    Build_Uce_Guardrails,
    Calculate_Read_Threading,
    Compute_Uce_Iteration,
    Load_Uce_Reference_Manifest,
    Pass_Uce_Guardrails,
    Score_Contig,
)


def guard_args(**overrides):
    values = {
        "max_contig_length": 5000,
        "min_read_density": 0.003,
        "density_check_min_length": 1000,
        "max_depth_cv": 0,
        "max_depth_ratio": 0,
        "max_unsupported_fraction": 0,
    }
    values.update(overrides)
    return values


class UceScoringGuardrailTests(unittest.TestCase):
    def test_long_low_density_contig_is_rejected(self):
        self.assertFalse(
            Pass_Uce_Guardrails(
                contig_len=3000,
                read_density=0.001,
                depth_cv=0.2,
                max_depth_ratio=2,
                max_unsupported_span=500,
                guardrails=guard_args(),
            )
        )

    def test_short_low_density_contig_is_not_rejected_by_density_guardrail(self):
        self.assertTrue(
            Pass_Uce_Guardrails(
                contig_len=500,
                read_density=0.001,
                depth_cv=0.2,
                max_depth_ratio=2,
                max_unsupported_span=100,
                guardrails=guard_args(),
            )
        )

    def test_long_unsupported_contig_is_rejected_when_enabled(self):
        self.assertFalse(
            Pass_Uce_Guardrails(
                contig_len=3000,
                read_density=0.01,
                depth_cv=0.2,
                max_depth_ratio=2,
                max_unsupported_span=2700,
                guardrails=guard_args(max_unsupported_fraction=0.8),
            )
        )

    def test_uce_score_penalizes_weak_read_density(self):
        weak_long = ["A" * 3000, 1, 0, 10, 3, 2900, 0.5, 0.001, 0.967, 2, 0.1, 1.0, 0.2, 2000, 1.0]
        dense_short = ["A" * 800, 1, 0, 10, 8, 700, 0.5, 0.01, 0.875, 2, 0.1, 1.0, 0.9, 50, 2.0]

        self.assertGreater(Score_Contig(dense_short, "uce"), Score_Contig(weak_long, "uce"))

    def test_read_threading_summarizes_supported_bases_and_gaps(self):
        coverage_fraction, max_unsupported, median_depth = Calculate_Read_Threading(
            "AAAACCCCGGGG",
            4,
            {"AAAA": 2, "GGGG": 4},
        )

        self.assertAlmostEqual(coverage_fraction, 8 / 12)
        self.assertEqual(max_unsupported, 4)
        self.assertEqual(median_depth, 3)

    def test_dynamic_search_reduces_weak_loci_only(self):
        args = SimpleNamespace(
            assembly_mode="uce",
            uce_dynamic_search=True,
            uce_min_search_depth=512,
            uce_manifest={},
        )

        self.assertEqual(Compute_Uce_Iteration(4096, args, 4, 40, "uce-1"), 1024)
        self.assertEqual(Compute_Uce_Iteration(4096, args, 60, 1000, "uce-1"), 4096)

    def test_manifest_overrides_guardrails_and_dynamic_search(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest_path = Path(tmp) / "manifest.tsv"
            manifest_path.write_text(
                "locus\tmax_contig_length\tmin_read_density\tmin_search_depth\n"
                "uce-1\t1200\t0.01\t2048\n"
            )
            manifest = Load_Uce_Reference_Manifest(str(manifest_path))

        args = SimpleNamespace(
            uce_max_contig_length=5000,
            uce_min_read_density=0.003,
            uce_density_check_min_length=1000,
            uce_max_depth_cv=0,
            uce_max_depth_ratio=0,
            uce_max_unsupported_fraction=0,
            assembly_mode="uce",
            uce_dynamic_search=True,
            uce_min_search_depth=512,
            uce_manifest=manifest,
        )

        guardrails = Build_Uce_Guardrails(args, "uce-1")
        self.assertEqual(guardrails["max_contig_length"], 1200)
        self.assertEqual(guardrails["min_read_density"], 0.01)
        self.assertEqual(Compute_Uce_Iteration(4096, args, 4, 40, "uce-1"), 2048)


if __name__ == "__main__":
    unittest.main()
