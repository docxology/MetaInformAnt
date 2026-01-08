# Perception Module Specification

## Overview
This module models the transformation of physical stimuli into internal psychological representations and decisions.

## Mathematical Models

### Psychophysics

1.  **Weber's Law**: The just-noticeable difference is proportional to stimulus magnitude.
    $$ \frac{\Delta I}{I} = k $$
2.  **Fechner's Law**: Sensation scales logarithmically with intensity.
    $$ S = k \ln \left( \frac{I}{I_0} \right) $$
3.  **Stevens' Power Law**: Sensation scales as a power function.
    $$ S = k I^a $$
    *   $a < 1$: Compressive (brightness, loudness)
    *   $a > 1$: Expansive (electric shock)

### Signal Detection Theory (SDT)

SDT assumes internal Gaussian distributions for Noise $N(\mu_n, \sigma)$ and Signal+Noise $N(\mu_s, \sigma)$. We assume equal variance $\sigma_n = \sigma_s = 1$.

1.  **Sensitivity ($d'$)**: Distance between the means in standard deviation units.
    $$ d' = Z(\text{Hit Rate}) - Z(\text{False Alarm Rate}) $$
2.  **Bias ($c$)**: Distance of the criterion from the intersection point.
    $$ c = -\frac{Z(H) + Z(FA)}{2} $$

## Implementation Details
- **Error Handling**: Hit/FA rates of 0 or 1 result in infinite Z-scores. The implementation applies a standard correction (default `1e-6` or `1/(2N)` adjustment logic) to clip values within `(epsilon, 1-epsilon)`.
