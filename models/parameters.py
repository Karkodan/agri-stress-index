# models/parameters.py

from dataclasses import dataclass, field
from typing import Optional

@dataclass(frozen=True)
class ModelParameters:
    leaf_area_index: float = 3.5               # Unitless
    canopy_height: float = 2.0                 # meters
    solar_radiation: float = 800.0             # W/m²
    wind_speed: float = 2.5                    # m/s
    air_temperature: float = 298.15            # Kelvin
    relative_humidity: float = 0.60            # 60%
    soil_moisture: float = 0.25                # m³/m³ (must be > 0 for CSI)
    canopy_temp: float = 305.15                # Kelvin
    ambient_temp: float = 300.15               # Kelvin

    def summary(self) -> str:
        return (
            f"🌿 LAI: {self.leaf_area_index}, "
            f"Height: {self.canopy_height}m, "
            f"Radiation: {self.solar_radiation}W/m², "
            f"Wind: {self.wind_speed}m/s, "
            f"Air Temp: {self.air_temperature}K, "
            f"Humidity: {self.relative_humidity * 100:.0f}%, "
            f"Soil Moisture: {self.soil_moisture}, "
            f"Canopy Temp: {self.canopy_temp}K, "
            f"Ambient Temp: {self.ambient_temp}K"
        )

# Default instance for validation and CSI logic
DEFAULT_PARAMETERS = ModelParameters()
