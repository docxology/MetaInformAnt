import os
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

# Add the scripts/rna directory to the path so we can import generate_orthologs
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../scripts/rna")))
import generate_orthologs


class TestOrthologGeneration(unittest.TestCase):

    def test_get_taxons_from_configs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)

            # Create a minimal regular amalgkit config
            with open(tmp_path / "amalgkit_good.yaml", "w") as f:
                yaml.dump({"taxon_id": 12345, "species_list": ["Good_Species"]}, f)

            # Create a config with string taxon
            with open(tmp_path / "amalgkit_string.yaml", "w") as f:
                yaml.dump({"taxon_id": "67890"}, f)

            # Create an ignored test config
            with open(tmp_path / "amalgkit_test.yaml", "w") as f:
                yaml.dump({"taxon_id": 99999}, f)

            # Create an ignored template config
            with open(tmp_path / "amalgkit_template.yaml", "w") as f:
                yaml.dump({"taxon_id": 10101}, f)

            taxons = generate_orthologs.get_taxons_from_configs(str(tmp_path), skip_tests=True)

            # We expect only the valid configs, sorted
            self.assertEqual(len(taxons), 2)
            self.assertIn("12345", taxons)
            self.assertIn("67890", taxons)
            self.assertNotIn("99999", taxons)
            self.assertNotIn("10101", taxons)
            self.assertEqual(taxons, ["12345", "67890"])

    def test_get_taxons_from_configs_no_skip(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)

            with open(tmp_path / "amalgkit_test.yaml", "w") as f:
                yaml.dump({"taxon_id": 99999}, f)

            taxons = generate_orthologs.get_taxons_from_configs(str(tmp_path), skip_tests=False)
            self.assertEqual(len(taxons), 1)
            self.assertIn("99999", taxons)


if __name__ == "__main__":
    unittest.main()
