"""Signal Detection Theory (SDT) metrics.

This module implements key metrics for analyzing decision making under uncertainty
using Signal Detection Theory (Sensitivity d', Bias).
"""

from __future__ import annotations

import math
import numpy as np
from typing import Tuple, Union, Optional

def d_prime(hit_rate: float, false_alarm_rate: float, 
            correction: float = 1e-6) -> float:
    """Calculate sensitivity index d' (d-prime).
    
    d' = Z(Hit Rate) - Z(False Alarm Rate)
    
    Args:
        hit_rate: Proportion of signal trials correctly identified (0 to 1)
        false_alarm_rate: Proportion of noise trials incorrectly identified as signal
        correction: Small epsilon to avoid infinite Z-scores at rates of 0 or 1
        
    Returns:
        d' value. Higher values indicate better sensitivity.
        0 indicates inability to distinguish signal from noise.
    """
    from scipy.stats import norm
    
    # Clip rates to avoid infinity
    h = np.clip(hit_rate, correction, 1 - correction)
    fa = np.clip(false_alarm_rate, correction, 1 - correction)
    
    return norm.ppf(h) - norm.ppf(fa)

def criterion_c(hit_rate: float, false_alarm_rate: float, 
               correction: float = 1e-6) -> float:
    """Calculate decision criterion c (bias).
    
    c = -(Z(Hit Rate) + Z(False Alarm Rate)) / 2
    
    Args:
        hit_rate: Proportion of signal trials correctly identified
        false_alarm_rate: Proportion of noise trials incorrectly identified
        correction: Small epsilon
        
    Returns:
        c value. 
        c = 0: Unbiased
        c > 0: Conservative (tendency to say "No")
        c < 0: Liberal (tendency to say "Yes")
    """
    from scipy.stats import norm
    
    h = np.clip(hit_rate, correction, 1 - correction)
    fa = np.clip(false_alarm_rate, correction, 1 - correction)
    
    return -(norm.ppf(h) + norm.ppf(fa)) / 2.0

def likelihood_ratio_beta(hit_rate: float, false_alarm_rate: float,
                         correction: float = 1e-6) -> float:
    """Calculate likelihood ratio β (beta).
    
    β = exp(d' * c)
    
    Args:
        hit_rate: Hit rate
        false_alarm_rate: False alarm rate
        correction: Small epsilon
        
    Returns:
        β value.
        β = 1: Unbiased
        β > 1: Conservative
        β < 1: Liberal
    """
    d = d_prime(hit_rate, false_alarm_rate, correction)
    c = criterion_c(hit_rate, false_alarm_rate, correction)
    
    return math.exp(d * c)

def sdt_metrics(hits: int, misses: int, false_alarms: int, correct_rejections: int) -> dict[str, float]:
    """Calculate standard SDT metrics from raw counts.
    
    Args:
        hits: Signal present, response "Yes"
        misses: Signal present, response "No"
        false_alarms: Signal absent, response "Yes"
        correct_rejections: Signal absent, response "No"
        
    Returns:
        Dictionary containing hit_rate, fa_rate, d_prime, c, beta
    """
    n_signal = hits + misses
    n_noise = false_alarms + correct_rejections
    
    # Handle cases with no trials
    if n_signal == 0: n_signal = 1
    if n_noise == 0: n_noise = 1
    
    hit_rate = hits / n_signal
    fa_rate = false_alarms / n_noise
    
    return {
        "hit_rate": hit_rate,
        "false_alarm_rate": fa_rate,
        "d_prime": d_prime(hit_rate, fa_rate),
        "c": criterion_c(hit_rate, fa_rate),
        "beta": likelihood_ratio_beta(hit_rate, fa_rate)
    }
