# app/api.py

from fastapi import APIRouter
from pydantic import BaseModel, Field
from datetime import datetime
from app import csi_logic
import logging

router = APIRouter()

# Configure logging for audit trails
logger = logging.getLogger("csi_logger")
logging.basicConfig(level=logging.INFO)

class CSIRequest(BaseModel):
    canopy_temp: float = Field(..., description="Canopy temperature in °C")
    ambient_temp: float = Field(..., description="Ambient air temperature in °C")
    soil_moisture: float = Field(..., description="Volumetric soil moisture (%)")

@router.post(
    "/",
    summary="Compute Crop Stress Index (CSI)",
    description="""
Calculates the Crop Stress Index (CSI) using canopy temperature, ambient temperature, and soil moisture.

**Use Cases:**
- Field-level stress diagnostics
- Mobile dashboard integration
- Scientific modeling and investor demos

**Returns:**
- CSI value
- Input echo
- Timestamp
- API version
""",
    response_description="Structured CSI result with metadata",
    tags=["CSI"]
)
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
