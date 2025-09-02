# validate_config.py

from models.constants import STEFAN_BOLTZMANN, SPECIFIC_HEAT_AIR, AIR_DENSITY
from models.parameters import DEFAULT_PARAMETERS

def test_constants():
    assert abs(STEFAN_BOLTZMANN - 5.67e-8) < 1e-10
    assert SPECIFIC_HEAT_AIR > 900
    assert AIR_DENSITY > 1.0
    print("✅ Constants validated.")

def test_parameters():
    required_fields = [
        "leaf_area_index",
        "canopy_height",
        "canopy_temp",
        "ambient_temp",
        "soil_moisture"
    ]
    for field in required_fields:
        assert hasattr(DEFAULT_PARAMETERS, field), f"Missing field: {field}"
        value = getattr(DEFAULT_PARAMETERS, field)
        assert value is not None, f"{field} is None"
        if isinstance(value, (int, float)):
            assert value > 0, f"{field} must be > 0"

    print("✅ Parameters validated.")

if __name__ == "__main__":
    test_constants()
    test_parameters()
