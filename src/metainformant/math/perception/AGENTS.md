# Agents Guide for Perception Module

## Context
This module serves cognitive modeling tasks where sensory processing is a component.

## Best Practices
- **SDT Inputs**: Always ensure hit rates and false alarm rates are probabilities [0.0, 1.0].
- **Clipping**: Be aware that `d_prime` automatically clips values to avoid infinity. If precise handling of 0/1 rates is needed (e.g., Log-Linear correction), implement it manually before calling or check the `correction` parameter.
- **Arrays**: Psychophysics functions accept numpy arrays for efficient batch processing of stimuli.

## Common Pitfalls
- **Weber Contrast**: Undefined for background intensity of 0.
- **d' Interpretation**: A d' of 0 means chance performance. Negative d' usually implies label swapping or significant confusion.
