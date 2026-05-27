from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from unix_command import run_gene_tree_job


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


if __name__ == "__main__":
    unittest.main()
