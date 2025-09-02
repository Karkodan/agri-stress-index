from models.canopy_temp import calculate_virtual_canopy_temp
from models.leaf_temp import compute_leaf_temperature
from models.radiation import compute_canopy_radiation
from models.soil_moisture import normalize_soil_moisture

def compute_crop_stress_index(
    air_temp_c: float,
    relative_humidity_pct: float,
    solar_radiation_w_m2: float,
    wind_speed_m_s: float,
    soil_moisture_pct: float
) -> dict:
    canopy_temp = calculate_virtual_canopy_temp(
        air_temp_c,
        relative_humidity_pct,
        solar_radiation_w_m2,
        wind_speed_m_s,
        soil_moisture_pct
    )

    stress_score = (canopy_temp / 40) * (1 - soil_moisture_pct / 100) * (relative_humidity_pct / 100)

    return {
        "canopy_temp": round(canopy_temp, 2),
        "stress_score": round(stress_score, 2),
        "risk_level": "High" if stress_score > 0.7 else "Moderate"
    }
