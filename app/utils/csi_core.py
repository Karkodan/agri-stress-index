# csi_core.py

def compute_csi(canopy_temp: float, ambient_temp: float, soil_moisture: float) -> float | None:
    try:
        normalized_temp = canopy_temp - ambient_temp
        if soil_moisture == 0:
            return None  # Avoid division by zero
        csi_score = normalized_temp / soil_moisture
        return round(csi_score, 3)
    except Exception:
        return None
