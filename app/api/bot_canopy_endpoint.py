# api/bot_canopy_endpoint.py

from fastapi import APIRouter, Query, HTTPException
from services.weather_fetcher import fetch_weather
from models.canopy_solver import solve_canopy_temperature

router = APIRouter()

@router.get("/bot/canopy-temp-live")
def get_live_canopy_temp(lat: float = Query(...), lon: float = Query(...)):
    weather = fetch_weather(lat, lon)

    env = {
        "LAT": lat,
        "DOY": weather["doy"],
        "tm": weather["hour_of_day"],
        "R_inc": weather["solar_radiation_w_m2"],
        "Tair": weather["air_temp_c"],
        "vpd": weather["vpd_kpa"],
        "pressure": 101.3,
        "cloud_frac": 0.3
    }

    result = solve_canopy_temperature(env)

    if not result or "leaf_temperature" not in result:
        raise HTTPException(status_code=500, detail="Canopy solver failed")

    return {
        "location": f"{lat},{lon}",
        "canopy_temperature": round(result["leaf_temperature"], 2),
        "Q_short": round(result["Q_short"], 2),
        "Qp": round(result["Qp"], 2),
        "Q_long": round(result["Q_long"], 2)
    }
