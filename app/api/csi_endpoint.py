from fastapi import APIRouter, HTTPException, Header
from pydantic import BaseModel, root_validator
from typing import Optional, Literal, Dict
from datetime import datetime
from app.utils.csi_core import compute_csi
import logging

router = APIRouter(tags=["Crop Stress Index"])

logger = logging.getLogger("csi_logger")
logger.setLevel(logging.INFO)

STRESS_SCALE: Dict[str, str] = {
    "0–20": "No stress",
    "20–40": "Mild stress",
    "40–60": "Moderate stress",
    "60–100": "Severe stress"
}

class FlexibleCSIRequest(BaseModel):
    lat: Optional[float]
    lon: Optional[float]
    timestamp: Optional[str]
    leaf_area_index: Optional[float]
    canopy_height: Optional[float]
    solar_radiation: Optional[float]
    wind_speed: Optional[float]
    air_temperature: Optional[float]
    relative_humidity: Optional[float]
    soil_moisture: Optional[float]
    canopy_temperature: Optional[float]
    ambient_temperature: Optional[float]

    @root_validator
    def validate_minimum_fields(cls, values):
        if not values.get("soil_moisture"):
            raise ValueError("soil_moisture is required")
        if not (
            values.get("canopy_temperature") and values.get("ambient_temperature")
        ) and not (
            values.get("leaf_area_index") and values.get("solar_radiation")
        ):
            raise ValueError("Provide either canopy+ambient temps or full model inputs")
        if values.get("timestamp"):
            try:
                datetime.fromisoformat(values["timestamp"])
            except ValueError:
                raise ValueError("Timestamp must be ISO 8601 format")
        return values

class CSIInputUsed(BaseModel):
    lat: Optional[float]
    lon: Optional[float]
    timestamp: datetime
    canopy_temperature: Optional[float]
    ambient_temperature: Optional[float]
    soil_moisture: Optional[float]
    leaf_area_index: Optional[float]
    canopy_height: Optional[float]
    solar_radiation: Optional[float]
    wind_speed: Optional[float]
    air_temperature: Optional[float]
    relative_humidity: Optional[float]

class CSIResponse(BaseModel):
    CSI: float
    stress_level: Literal["No stress", "Mild stress", "Moderate stress", "Severe stress"]
    recommendation: str
    scale: Dict[str, str]
    input_used: CSIInputUsed
    timestamp: str

def interpret_stress(score: float):
    if score < 20:
        return "No stress", "Normal irrigation"
    elif score < 40:
        return "Mild stress", "Monitor closely"
    elif score < 60:
        return "Moderate stress", "Consider irrigation or shade intervention"
    else:
        return "Severe stress", "Immediate action needed"

@router.post("/", response_model=CSIResponse, summary="Compute CSI from flexible input")
def compute_adaptive_csi(payload: FlexibleCSIRequest, x_api_key: Optional[str] = Header(None)):
    if x_api_key != "your-secret-key":
        raise HTTPException(status_code=403, detail="Invalid API key")

    try:
        score = compute_csi(
            canopy_temp=payload.canopy_temperature,
            ambient_temp=payload.ambient_temperature,
            soil_moisture=payload.soil_moisture
        )
        score_rounded = round(float(score), 2)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"CSI computation failed: {str(e)}")

    parsed_timestamp = (
        datetime.fromisoformat(payload.timestamp)
        if payload.timestamp
        else datetime.now()
    )

    input_used = CSIInputUsed(
        lat=payload.lat,
        lon=payload.lon,
        timestamp=parsed_timestamp,
        canopy_temperature=payload.canopy_temperature,
        ambient_temperature=payload.ambient_temperature,
        soil_moisture=payload.soil_moisture,
        leaf_area_index=payload.leaf_area_index,
        canopy_height=payload.canopy_height,
        solar_radiation=payload.solar_radiation,
        wind_speed=payload.wind_speed,
        air_temperature=payload.air_temperature,
        relative_humidity=payload.relative_humidity,
    )

    stress_level, recommendation = interpret_stress(score_rounded)
    logger.info(f"CSI computed: {score_rounded} | Stress: {stress_level} | Lat: {payload.lat}, Lon: {payload.lon}")

    return CSIResponse(
        CSI=score_rounded,
        stress_level=stress_level,
        recommendation=recommendation,
        scale=STRESS_SCALE,
        input_used=input_used,
        timestamp=parsed_timestamp.isoformat()
    )
