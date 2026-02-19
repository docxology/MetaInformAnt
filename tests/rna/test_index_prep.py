"""
Zero-Mock Unit Test for IndexComplexityManager.

Verifies that the IndexComplexityManager correctly filters:
1. Low-complexity / non-coding sequences (XR_)
2. Short fragments (< 200bp)
3. Ribosomal RNA (rRNA)
"""

import gzip
import shutil
import tempfile
import unittest
from pathlib import Path

from metainformant.rna.amalgkit.index_prep import IndexComplexityManager

class TestIndexComplexityManager(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        self.input_fasta = Path(self.test_dir) / "test_input.fna"
        self.output_fasta = Path(self.test_dir) / "test_output.fna"

    def tearDown(self):
        # Clean up
        shutil.rmtree(self.test_dir)

    def test_filter_fasta(self):
        """Test filtering logic with a dummy FASTA file."""
        # Create dummy FASTA content
        content = [
            ">valid_transcript_1 len=250",
            "A" * 250,  # Keep (Valid length, no bad patterns)
            ">XR_bad_transcript_1 len=300",
            "G" * 300,  # discard (XR_ prefix)
            ">short_transcript_1 len=50",
            "C" * 50,   # discard (Short < 200bp)
            ">ribosomal_rna_1",
            "T" * 250,  # discard (rRNA in header)
            ">NR_bad_transcript_2 len=250",
            "A" * 250,  # discard (NR_ prefix)
            ">valid_transcript_2 len=200",
            "A" * 200,  # Keep (Exactly 200bp)
        ]
        
        # Write to file
        with open(self.input_fasta, "w") as f:
            f.write("\n".join(content))

        # Initialize manager
        manager = IndexComplexityManager(min_length=200, exclude_patterns=["XR_", "NR_", "rRNA", "ribosomal"])
        
        # Run filtering
        stats = manager.filter_fasta(self.input_fasta, self.output_fasta)

        # Verify stats
        self.assertEqual(stats["total"], 6, "Should process 6 entries")
        self.assertEqual(stats["kept"], 2, "Should keep 2 entries (valid_1, valid_2)")
        self.assertEqual(stats["removed_pattern"], 3, "Should remove 3 patterns (XR, NR, rRNA)")
        self.assertEqual(stats["removed_length"], 1, "Should remove 1 short entry")

        # Verify output file content
        with open(self.output_fasta, "r") as f:
            output_content = f.read()
            
        self.assertIn(">valid_transcript_1", output_content)
        self.assertIn(">valid_transcript_2", output_content)
        self.assertNotIn(">XR_", output_content)
        self.assertNotIn(">short_", output_content)
        self.assertNotIn(">ribosomal_", output_content)
        self.assertNotIn(">NR_", output_content)

if __name__ == "__main__":
    unittest.main()
