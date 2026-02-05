"""Model interpretability and explainability subpackage.

This subpackage provides tools for understanding and explaining machine
learning model predictions, including permutation importance, SHAP-like
explanations, LIME, partial dependence, and advanced feature selection.
"""

from __future__ import annotations

from . import explainers, feature_selection

__all__ = ["explainers", "feature_selection"]
