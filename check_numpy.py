import sys
try:
    import numpy
    print(f"Numpy available: {numpy.__version__}")
    import metainformant.gwas.analysis.structure as structure
    print("Structure module imported")
except ImportError as e:
    print(f"ImportError: {e}")
except Exception as e:
    print(f"Error: {e}")
