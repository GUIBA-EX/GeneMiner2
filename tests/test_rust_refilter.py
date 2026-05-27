from pathlib import Path
import filecmp
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class RustRefilterGoldenTests(unittest.TestCase):
    @unittest.skipIf(shutil.which("cargo") is None, "cargo is not installed")
    def test_rust_refilter_matches_python_keep_linked_mates(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            ref_dir = work / "ref"
            pe_dir = work / "filtered_pe"
            py_out = work / "py_out"
            rs_out = work / "rs_out"
            ref_dir.mkdir()
            pe_dir.mkdir()

            (ref_dir / "gene.fasta").write_text(">gene_ref\nACGTACGTACGTACGTACGT\n")
            (pe_dir / "gene_1.fq").write_text(
                "@read1/1\nACGTACGTACGT\n+\nFFFFFFFFFFFF\n"
                "@read2/1\nNNNNNNNNNNNN\n+\nFFFFFFFFFFFF\n"
            )
            (pe_dir / "gene_2.fq").write_text(
                "@read1/2\nNNNNNNNNNNNN\n+\nFFFFFFFFFFFF\n"
                "@read2/2\nNNNNNNNNNNNN\n+\nFFFFFFFFFFFF\n"
            )

            common_args = [
                "-r",
                str(ref_dir),
                "-qd",
                str(pe_dir),
                "-kf",
                "5",
                "--min-depth",
                "0",
                "--max-depth",
                "999999",
                "--max-size",
                "999999",
                "--keep-linked-mates",
                "-p",
                "1",
            ]

            subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "main_refilter_new.py"), "-o", str(py_out), *common_args],
                check=True,
                cwd=ROOT,
            )
            subprocess.run(
                [
                    "cargo",
                    "run",
                    "--quiet",
                    "--manifest-path",
                    str(ROOT / "rust" / "main_refilter_new" / "Cargo.toml"),
                    "--",
                    "-o",
                    str(rs_out),
                    *common_args,
                ],
                check=True,
                cwd=ROOT,
            )

            self.assertTrue(filecmp.cmp(py_out / "gene.fq", rs_out / "gene.fq", shallow=False))


if __name__ == "__main__":
    unittest.main()
