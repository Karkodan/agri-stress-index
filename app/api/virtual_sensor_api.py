from fastapi import FastAPI
from pydantic import BaseModel
from virtual_sensor import calculate_virtual_canopy_temp

app = FastAPI(title="AgriGenAI Virtual Sensor")

class SensorInput(BaseModel):
    air_temp_c: float
    relative_humidity_pct: float
    solar_radiation_w_m2: float
    wind_speed_m_s: float
    soil_moisture_pct: float

@app.post("/canopy-temp")
def get_canopy_temp(data: SensorInput):
    result = calculate_virtual_canopy_temp(**data.dict())
    return {"canopy_temperature": round(result, 2)}
