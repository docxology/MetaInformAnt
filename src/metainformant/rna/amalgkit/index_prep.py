"""
Index Preparation and Complexity Management.

This module provides utilities to manage the complexity of transcriptome indices
for pseudoaligners like kallisto. It handles filtering of:
1. Low-complexity / non-coding sequences (XR_, NR_)
2. Short fragments
3. Ribosomal RNA (rRNA)
4. Duplicate sequences (future)
"""

import gzip
import logging
import shutil
from pathlib import Path
from typing import List, Optional, Union

# Use internal logging
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

class IndexComplexityManager:
    """Manages transcriptome complexity by filtering problematic sequences."""

    def __init__(self, min_length: int = 200, exclude_patterns: Optional[List[str]] = None):
        """
        Initialize the manager.

        Args:
            min_length: Minimum sequence length to keep (default: 200).
            exclude_patterns: List of header substrings to exclude (default: ["XR_", "NR_", "rRNA", "ribosomal"]).
        """
        self.min_length = min_length
        self.exclude_patterns = exclude_patterns or ["XR_", "NR_", "rRNA", "ribosomal RNA"]

    def filter_fasta(self, input_path: Union[str, Path], output_path: Union[str, Path]) -> dict:
        """
        Filter a FASTA file (gzipped or plain) based on complexity rules.

        Args:
            input_path: Path to input FASTA (can be .gz).
            output_path: Path to write filtered FASTA.

        Returns:
            Dictionary with statistics: {
                "total": int,
                "kept": int,
                "removed_pattern": int,
                "removed_length": int
            }
        """
        input_path = Path(input_path)
        output_path = Path(output_path)
        
        stats = {
            "total": 0,
            "kept": 0,
            "removed_pattern": 0,
            "removed_length": 0
        }

        # Determine open function based on extension
        is_gzip = input_path.suffix == ".gz"
        open_func = gzip.open if is_gzip else open
        mode = "rt" if is_gzip else "r"

        try:
            with open_func(input_path, mode) as f_in, open(output_path, "w") as f_out:
                header = None
                sequence = []

                def process_entry(h, s):
                    stats["total"] += 1
                    s_str = "".join(s)

                    # Check 1: Exclude patterns
                    for pattern in self.exclude_patterns:
                        if pattern in h:
                            stats["removed_pattern"] += 1
                            return

                    # Check 2: Length
                    if len(s_str) < self.min_length:
                        stats["removed_length"] += 1
                        return

                    # Keep
                    f_out.write(f"{h}\n{s_str}\n")
                    stats["kept"] += 1

                for line in f_in:
                    line = line.strip()
                    if not line:
                        continue
                        
                    if line.startswith(">"):
                        if header:
                            process_entry(header, sequence)
                        header = line
                        sequence = []
                    else:
                        sequence.append(line)
                
                # Process last entry
                if header:
                    process_entry(header, sequence)

            logger.info(f"Filtering complete for {input_path.name}: {stats}")
            return stats

        except Exception as e:
            logger.error(f"Failed to filter FASTA {input_path}: {e}")
            raise
