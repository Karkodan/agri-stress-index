def calculate_virtual_canopy_temp(
    air_temp_c: float,
    relative_humidity_pct: float,
    solar_radiation_w_m2: float,
    wind_speed_m_s: float,
    soil_moisture_pct: float
) -> float:
    """
    Placeholder canopy temperature model.
    Replace with HESS2021 logic or MATLAB bridge later.
    """
    # Simple proxy logic for now
    canopy_temp = (
        air_temp_c +
        (solar_radiation_w_m2 / 1000) -
        (soil_moisture_pct / 10) -
        (wind_speed_m_s * 0.5)
    )
    return round(canopy_temp, 2)
