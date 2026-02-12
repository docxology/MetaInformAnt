"""Survival and time-to-event analysis subpackage.

Provides Kaplan-Meier estimation, Cox proportional hazards modelling,
competing risks analysis, recurrent event modelling, and time-varying
covariate handling."""
from __future__ import annotations

from . import time_to_event

__all__ = ['time_to_event']
