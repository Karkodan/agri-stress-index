def compute_csi(canopy_temp, ambient_temp, soil_moisture):
    if soil_moisture == 0:
        return None
    delta_temp = canopy_temp - ambient_temp
    return round(delta_temp / soil_moisture, 2)
