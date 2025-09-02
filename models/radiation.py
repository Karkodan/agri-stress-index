# models/radiation.py

import numpy as np
from models.constants import PhysicalConstants
from models.parameters import ModelParameters
from models.solar_geometry import solar_dec_azimuth  # to be implemented

def compute_canopy_radiation(LAT, R_inc, Ta, Tl, DOY, tm, cloud_frac=0.5, constants=None, parameters=None):
    """
    Compute canopy radiation balance.

    Args:
        LAT (float): Latitude (degrees)
        R_inc (float): Incoming solar radiation (W/m²)
        Ta (float): Air temperature (°C)
        Tl (float): Leaf temperature (°C)
        DOY (int): Day of year
        tm (float): Solar time (hours)
        cloud_frac (float): Cloud fraction [0–1]
        constants (PhysicalConstants): Optional override
        parameters (ModelParameters): Optional override

    Returns:
        tuple: (Q_short [W/m²], Qp [μmol/m²/s], Q_long [W/m²])
    """
    c = constants or PhysicalConstants()
    p = parameters or ModelParameters()

    # Convert temperatures to Kelvin
    Ta_K = Ta + 273.15 if Ta < 270 else Ta
    Tl_K = Tl + 273.15 if Tl < 270 else Tl

    # Solar geometry
    ZEN, AZIM, dd = solar_dec_azimuth(DOY, LAT, tm)
    R = max(0, R_inc * np.cos(ZEN))

    # Spectral fractions
    PAR_frac = 0.5
    NIR_frac = 0.5
    PAR_top = PAR_frac * R
    NIR_top = NIR_frac * R

    # Leaf optical properties
    sigm_PAR = 0.15
    sigm_NIR = 0.65
    rho_ch_PAR = (1 - np.sqrt(1 - sigm_PAR)) / (1 + np.sqrt(1 - sigm_PAR))
    rho_ch_NIR = (1 - np.sqrt(1 - sigm_NIR)) / (1 + np.sqrt(1 - sigm_NIR))

    # Extinction and reflection
    kd = 0.5
    kb = 1 / (2 * np.cos(ZEN))
    rho_PAR = rho_ch_PAR * 2 * kb / (kb + kd)
    rho_NIR = rho_ch_NIR * 2 * kb / (kb + kd)
    kestinctPAR = kb * np.sqrt(1 - sigm_PAR)
    kestinctNIR = kb * np.sqrt(1 - sigm_NIR)

    LAI = p.leaf_area_index
    Q_PAR = PAR_top * (1 - rho_PAR) * (1 - np.exp(-kestinctPAR * LAI))
    Q_NIR = NIR_top * (1 - rho_NIR) * (1 - np.exp(-kestinctNIR * LAI))
    q = 4.6  # μmol/J conversion factor
    Qp = q * Q_PAR
    Q_short = Q_PAR + Q_NIR

    # Longwave radiation
    sigm = c.stefan_boltzmann
    epsil = 0.97
    epsil_clear = 9.2e-6 * (273.16 + Ta_K)**2
    epsil_a = epsil_clear * (1 - 0.84 * cloud_frac) + 0.84 * cloud_frac
    kestinctLW = kd
    Q_long = (epsil_a * sigm * Ta_K**4 - epsil * sigm * Tl_K**4) * (1 - np.exp(-kestinctLW * LAI))

    return Q_short, Qp, Q_long
