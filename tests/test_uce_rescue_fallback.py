from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from unix_command import read_density_or_blank, rescue_density_decreased


class UceRescueFallbackTests(unittest.TestCase):
    def test_density_decrease_triggers_fallback(self):
        before = {'selected_contig_length': '326', 'read_count': '14'}
        after = {'selected_contig_length': '14426', 'read_count': '13'}
        self.assertTrue(rescue_density_decreased(before, after))

    def test_density_increase_keeps_rescue(self):
        before = {'selected_contig_length': '532', 'read_count': '10'}
        after = {'selected_contig_length': '1476', 'read_count': '38'}
        self.assertFalse(rescue_density_decreased(before, after))

    def test_blank_density_does_not_trigger_fallback(self):
        self.assertEqual(read_density_or_blank({'selected_contig_length': '0', 'read_count': '5'}), '')
        self.assertFalse(rescue_density_decreased({}, {'selected_contig_length': '10', 'read_count': '1'}))


if __name__ == "__main__":
    unittest.main()
