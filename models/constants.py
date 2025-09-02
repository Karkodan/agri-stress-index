# models/constants.py

from dataclasses import dataclass

@dataclass(frozen=True)
class PhysicalConstants:
    stefan_boltzmann: float = 5.67e-8
    specific_heat_air: float = 1005
    air_density: float = 1.225

# Optional: expose constants directly
STEFAN_BOLTZMANN = PhysicalConstants().stefan_boltzmann
SPECIFIC_HEAT_AIR = PhysicalConstants().specific_heat_air
AIR_DENSITY = PhysicalConstants().air_density
