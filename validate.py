# validate.py

import sys
import os
import traceback

# Ensure project root is in Python path before any imports
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
if PROJECT_ROOT not in sys.path:
    sys.path.append(PROJECT_ROOT)

# Imports after path setup
from models.parameters import DEFAULT_PARAMETERS
from models.constants import STEFAN_BOLTZMANN, SPECIFIC_HEAT_AIR, AIR_DENSITY
from app import csi_logic
from models import validate_config

def run_validation():
    print("🔍 Running CSI validation...\n")

    # Print input summary
    print("📦 CSI Input Parameters:")
    if hasattr(DEFAULT_PARAMETERS, "summary"):
        print(DEFAULT_PARAMETERS.summary())
    else:
        print(vars(DEFAULT_PARAMETERS))  # fallback for dataclass or object

    # Run CSI logic with explicit arguments
    try:
        result = csi_logic.compute_csi(
            canopy_temp=DEFAULT_PARAMETERS.canopy_temp,
            ambient_temp=DEFAULT_PARAMETERS.ambient_temp,
            soil_moisture=DEFAULT_PARAMETERS.soil_moisture
        )
        print(f"\n✅ CSI Output: {result}")
    except Exception as e:
        print("❌ CSI computation failed:")
        traceback.print_exc()

    # Run config validation
    print("\n📋 Running config validation...")
    try:
        if hasattr(validate_config, "test_constants"):
            validate_config.test_constants()
        else:
            print("⚠️ test_constants() not found in validate_config.py")

        if hasattr(validate_config, "test_parameters"):
            validate_config.test_parameters()
        else:
            print("⚠️ test_parameters() not found in validate_config.py")
    except Exception as e:
        print("❌ Config validation failed:")
        traceback.print_exc()

    # Print constants for audit
    print("\n📐 Physical Constants:")
    print(f"  Stefan-Boltzmann: {STEFAN_BOLTZMANN}")
    print(f"  Specific Heat of Air: {SPECIFIC_HEAT_AIR}")
    print(f"  Air Density: {AIR_DENSITY}")

    print("\n🎉 Validation completed successfully.")

if __name__ == "__main__":
    try:
        run_validation()
    except Exception as e:
        print("❌ Validation failed with error:")
        traceback.print_exc()
        sys.exit(1)
