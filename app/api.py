# app/api.py

from fastapi import APIRouter
from pydantic import BaseModel
from datetime import datetime
from app import csi_logic
import logging

router = APIRouter()

# Optional: configure logging for audit trails
logger = logging.getLogger("csi_logger")
logging.basicConfig(level=logging.INFO)

class CSIRequest(BaseModel):
    canopy_temp: float
    ambient_temp: float
    soil_moisture: float

@router.post("/compute-csi")
def compute_csi(data: CSIRequest):
    result = csi_logic.compute_csi(
        canopy_temp=data.canopy_temp,
        ambient_temp=data.ambient_temp,
        soil_moisture=data.soil_moisture
    )

    response = {
        "version": "1.0.0",
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "input": {
            "canopy_temp": data.canopy_temp,
            "ambient_temp": data.ambient_temp,
            "soil_moisture": data.soil_moisture
        },
        "output": {
            "csi": result
        },
        "status": "success"
    }

    logger.info(f"CSI Request: {data.dict()} → Result: {result}")
    return response
