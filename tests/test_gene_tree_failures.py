from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from unix_command import run_gene_tree_job, write_failed_gene_trees


class GeneTreeFailureTests(unittest.TestCase):
    def test_gene_tree_failure_is_reported_without_raising(self):
        def fail(_gene):
            raise RuntimeError("bad alignment")

        gene, alignment, tree_path, error = run_gene_tree_job("uce-1", "/tmp/aligned", fail)

        self.assertEqual(gene, "uce-1")
        self.assertEqual(alignment, "/tmp/aligned/uce-1.fasta")
        self.assertIsNone(tree_path)
        self.assertIn("bad alignment", error)

    def test_gene_tree_success_returns_tree_path(self):
        gene, alignment, tree_path, error = run_gene_tree_job("uce-2", "/tmp/aligned", lambda _gene: "/tmp/tree")

        self.assertEqual(gene, "uce-2")
        self.assertEqual(alignment, "/tmp/aligned/uce-2.fasta")
        self.assertEqual(tree_path, "/tmp/tree")
        self.assertEqual(error, "")

    def test_failed_gene_tree_report_is_removed_when_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "failed_gene_trees.tsv"
            report.write_text("old\n")

            write_failed_gene_trees(tmp, [])

            self.assertFalse(report.exists())


if __name__ == "__main__":
    unittest.main()
