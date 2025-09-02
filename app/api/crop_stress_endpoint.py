# api/crop_stress_endpoint.py

from fastapi import APIRouter, Request
from models.virtual_sensor import compute_crop_stress_index

router = APIRouter()

@router.post("/api/crop-stress")
async def crop_stress_endpoint(request: Request):
    payload = await request.json()
    result = compute_crop_stress_index(
        air_temp_c=payload["air_temp_c"],
        relative_humidity_pct=payload["relative_humidity_pct"],
        solar_radiation_w_m2=payload["solar_radiation_w_m2"],
        wind_speed_m_s=payload["wind_speed_m_s"],
        soil_moisture_pct=payload["soil_moisture_pct"]
    )
    return result
