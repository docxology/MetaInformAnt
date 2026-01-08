import numpy as np
import logging
import sys
from metainformant.math import ddm
from metainformant.math import coalescent

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("math_verification")

def verify_ddm_optimization():
    logger.info("Verifying DDM optimization...")
    
    # Generate synthetic data
    v_true = 0.8
    a_true = 1.5
    z_true = 0.5
    t0_true = 0.3
    
    # Simple simulation of RTs (Wiener process approximation)
    # This is just for testing the fitting machinery, not the model physics
    n_trials = 100
    rt_data = []
    choice_data = []
    
    for _ in range(n_trials):
        # Deterministic part + noise
        dt = a_true / v_true + t0_true + np.random.normal(0, 0.1)
        rt_data.append(abs(dt)) # Ensure positive
        choice_data.append(1) # Always upper boundary for simplicity in this synthetic test
        
    logger.info(f"Generated {n_trials} trials.")
    
    # Fit parameters
    params = ddm.fit_ddm_parameters(rt_data, choice_data)
    
    logger.info(f"Fitted parameters: {params}")
    
    if params['log_likelihood'] == float('-inf'):
        logger.error("DDM fitting failed (log likelihood is -inf)")
        return False
        
    logger.info("DDM optimization verified.")
    return True

def verify_coalescent_optimization():
    logger.info("Verifying Coalescent optimization...")
    
    # Test watterson_theta with small and large n
    # Note: s must be <= L
    theta_small = coalescent.watterson_theta(10, 5, sequence_length=10)
    logger.info(f"Theta (n=10, s=5, L=10): {theta_small}")
    
    theta_large = coalescent.watterson_theta(1000, 50, sequence_length=100)
    logger.info(f"Theta (n=1000, s=50, L=100): {theta_large}")
    
    if theta_small <= 0 or theta_large <= 0:
        logger.error("Watterson theta calculation failed")
        return False
        
    logger.info("Coalescent optimization verified.")
    return True

def verify_perception():
    logger.info("Verifying Perception submodule...")
    import math
    from metainformant.math.perception import d_prime, stevens_power_law, sdt_metrics
    
    # Test d_prime
    # H=0.84 (Z=+1), FA=0.16 (Z=-1) => d' = 1 - (-1) = 2
    dp = d_prime(0.8413, 0.1587)
    logger.info(f"d'(H=0.84, FA=0.16) = {dp:.4f} (Expected ~2.0)")
    if abs(dp - 2.0) > 0.01:
        logger.error("d_prime calculation incorrect")
        return False
        
    # Test Stevens Power Law (vectorized)
    intensities = np.array([10, 100])
    sensation = stevens_power_law(intensities, exponent=0.5) # sqrt
    logger.info(f"Stevens(I=[10,100], a=0.5) = {sensation}")
    if not np.allclose(sensation, [math.sqrt(10), 10.0]):
         logger.error("Stevens law calculation incorrect")
         return False
         
    return True

def main():
    success = True
    if not verify_ddm_optimization():
        success = False
    if not verify_coalescent_optimization():
        success = False
    if not verify_perception():
        success = False
        
    if success:
        logger.info("All verifications passed!")
        sys.exit(0)
    else:
        logger.error("Some verifications failed.")
        sys.exit(1)

if __name__ == "__main__":
    main()
