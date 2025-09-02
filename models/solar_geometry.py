# models/solar_geometry.py

import numpy as np

def solar_dec_azimuth(DOY, LAT, tm):
    """
    Compute solar zenith and azimuth angles.

    Args:
        DOY (int): Day of year
        LAT (float): Latitude (degrees)
        tm (float): Solar time (hours)

    Returns:
        tuple: (ZEN [radians], AZIM [radians], declination [radians])
    """
    # Convert latitude to radians
    lat_rad = np.radians(LAT)

    # Solar declination angle (Cooper's formula)
    decl = 23.45 * np.sin(np.radians(360 * (DOY - 81) / 365))
    decl_rad = np.radians(decl)

    # Hour angle
    hour_angle = np.radians(15 * (tm - 12))

    # Zenith angle
    cos_ZEN = (np.sin(lat_rad) * np.sin(decl_rad) +
               np.cos(lat_rad) * np.cos(decl_rad) * np.cos(hour_angle))
    cos_ZEN = np.clip(cos_ZEN, -1, 1)
    ZEN = np.arccos(cos_ZEN)

    # Azimuth angle (optional, placeholder)
    AZIM = np.arctan2(-np.sin(hour_angle),
                      np.tan(decl_rad) * np.cos(lat_rad) -
                      np.sin(lat_rad) * np.cos(hour_angle))

    return ZEN, AZIM, decl_rad
