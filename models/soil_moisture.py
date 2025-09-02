def normalize_soil_moisture(raw_value: float) -> float:
    """Clamp and normalize soil moisture to 0–100%"""
    return max(0.0, min(100.0, raw_value))
