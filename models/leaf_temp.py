# models/leaf_temp.py

import numpy as np
from models.constants import PhysicalConstants
from models.parameters import ModelParameters

def compute_leaf_temperature(pa, gh, Tair, D, gv, es, Qshort, csi, cloud_frac=0.5, constants=None, parameters=None):
    """
    Compute leaf temperature (°C) from energy balance.

    Args:
        pa (float): Atmospheric pressure (kPa)
        gh (float): Canopy conductance to heat (mol/m²/s)
        Tair (float): Air temperature (°C)
        D (float): Vapor pressure deficit (kPa)
        gv (float): Canopy conductance to water vapor (mol/m²/s)
        es (float): Saturation vapor pressure (kPa)
        Qshort (float): Absorbed shortwave radiation (W/m²)
        csi (float or np.ndarray): Cumulated leaf area above layer
        cloud_frac (float): Cloud fraction [0–1]
        constants (PhysicalConstants): Optional override
        parameters (ModelParameters): Optional override

    Returns:
        float: Leaf temperature (°C)
    """
    c = constants or PhysicalConstants()
    p = parameters or ModelParameters()

    # Slope of vapor pressure curve (Pa/K)
    delta = 17.502 * 240.97 * es / (240.97 + Tair)**2
    ss = delta / pa

    # Clear sky emissivity (Swinbank)
    epsaclear = 9.2e-6 * (Tair + 273.15)**2
    epsa = epsaclear * (1 - 0.84 * cloud_frac) + 0.84 * cloud_frac

    # Radiation terms
    T_k = Tair + 273.15
    sigm = c.stefan_boltzmann
    Cp = c.specific_heat_air
    P = pa * 1000  # Convert kPa to Pa
    kd = 0.5       # Extinction coefficient (assumed)
    epsil = 0.97   # Leaf emissivity (assumed)
    lambd = 2.45e6 # Latent heat of vaporization (J/kg)

    if isinstance(csi, (list, np.ndarray)) and len(csi) > 2:
        numer = Qshort + (epsa * sigm * T_k**4 - epsil * sigm * T_k**4) * kd * np.exp(-kd * csi) - lambd * gv * D / P
        denomin = Cp * gh + lambd * gv * ss + 4 * sigm * epsil * T_k**3 * kd * np.exp(-kd * csi)
    else:
        LAI = p.leaf_area_index
        numer = Qshort + (epsa * sigm * T_k**4 - epsil * sigm * T_k**4) * (1 - np.exp(-kd * LAI)) - lambd * gv * D / P
        denomin = Cp * gh + lambd * gv * ss + 4 * sigm * epsil * T_k**3 * (1 - np.exp(-kd * LAI))

    Tlnew = Tair + numer / denomin
    return Tlnew
