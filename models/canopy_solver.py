# models/canopy_solver.py

from models.constants import PhysicalConstants
from models.parameters import ModelParameters
from models.radiation import compute_canopy_radiation
from models.leaf_temp import compute_leaf_temperature

def solve_canopy_temperature(env_inputs: dict, constants=None, parameters=None) -> dict:
    """
    Unified canopy temperature solver.

    Args:
        env_inputs (dict): Environmental inputs:
            - LAT, DOY, tm, R_inc, Tair, vpd, pressure, cloud_frac
        constants (PhysicalConstants): Optional override
        parameters (ModelParameters): Optional override

    Returns:
        dict: {
            "Q_short": ..., "Qp": ..., "Q_long": ...,
            "leaf_temperature": ...
        }
    """
    c = constants or PhysicalConstants()
    p = parameters or ModelParameters()

    try:
        # Step 1: Compute radiation
        Q_short, Qp, Q_long = compute_canopy_radiation(
            LAT=env_inputs["LAT"],
            R_inc=env_inputs["R_inc"],
            Ta=env_inputs["Tair"],
            Tl=env_inputs["Tair"],  # initial guess
            DOY=env_inputs["DOY"],
            tm=env_inputs["tm"],
            cloud_frac=env_inputs.get("cloud_frac", 0.5),
            constants=c,
            parameters=p
        )

        # Step 2: Compute leaf temperature
        leaf_temp = compute_leaf_temperature(
            pa=env_inputs.get("pressure", 101.3),
            gh=0.02,
            Tair=env_inputs["Tair"],
            D=env_inputs.get("vpd", 2.0),
            gv=0.03,
            es=4.2,
            Qshort=Q_short,
            csi=1.5,
            cloud_frac=env_inputs.get("cloud_frac", 0.5),
            constants=c,
            parameters=p
        )

        return {
            "Q_short": Q_short,
            "Qp": Qp,
            "Q_long": Q_long,
            "leaf_temperature": leaf_temp
        }

    except Exception as e:
        print("Canopy solver error:", e)
        return None
