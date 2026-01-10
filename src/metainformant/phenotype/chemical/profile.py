from typing import Dict, List, Any
import math
from .compound import Compound
from metainformant.core.utils.errors import ValidationError

class ChemicalProfile:
    """
    Represents a chemical profile (e.g., from GC-MS).
    
    Attributes:
        sample_id (str): Unique identifier for the sample.
        compounds (Dict[Compound, float]): Mapping of Compound to abundance/concentration.
        metadata (Dict[str, Any]): Experimental conditions.
    """
    
    def __init__(self, sample_id: str, compounds: Dict[Compound, float], metadata: Dict[str, Any] = None):
        self.sample_id = sample_id
        self.compounds = compounds
        self.metadata = metadata or {}
        
    def normalize(self, method: str = "total_ion_current") -> 'ChemicalProfile':
        """
        Return a new normalized profile.
        
        Args:
            method: 'total_ion_current' (sum to 1) or 'max_peak' (max to 1).
        """
        if not self.compounds:
            return ChemicalProfile(self.sample_id, {}, self.metadata)
            
        if method == "total_ion_current":
            total = sum(self.compounds.values())
            if total == 0:
                return ChemicalProfile(self.sample_id, self.compounds.copy(), self.metadata)
            new_compounds = {k: v / total for k, v in self.compounds.items()}
            
        elif method == "max_peak":
            max_val = max(self.compounds.values())
            if max_val == 0:
                return ChemicalProfile(self.sample_id, self.compounds.copy(), self.metadata)
            new_compounds = {k: v / max_val for k, v in self.compounds.items()}
            
        else:
            raise ValidationError(f"Unknown normalization method: {method}")
            
        return ChemicalProfile(self.sample_id, new_compounds, self.metadata)
        
    def bray_curtis_distance(self, other: 'ChemicalProfile') -> float:
        """
        Calculate Bray-Curtis dissimilarity between two profiles.
        Assumes both are normalized similarly.
        """
        all_compounds = set(self.compounds.keys()) | set(other.compounds.keys())
        
        diff_sum = 0.0
        total_sum = 0.0
        
        for cmp in all_compounds:
            val1 = self.compounds.get(cmp, 0.0)
            val2 = other.compounds.get(cmp, 0.0)
            
            diff_sum += abs(val1 - val2)
            total_sum += (val1 + val2)
            
        if total_sum == 0:
            return 0.0
            
        return diff_sum / total_sum
