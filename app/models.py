from pydantic import BaseModel
from typing import Optional
from datetime import datetime

# 🌾 Request model for CSI computation
class CSIRequest(BaseModel):
    lat: float
    lon: float
    timestamp: datetime  # Accepts ISO 8601 strings like "2025-09-01T18:32:00"
    canopy_temp: float
    ambient_temp: float
    soil_moisture: float

# 🌿 Echoed input block in response
class CSIInputUsed(BaseModel):
    lat: float
    lon: float
    timestamp: datetime
    canopy_temp: float
    ambient_temp: float
    soil_moisture: float

# 🌞 Response model for CSI API
class CSIResponse(BaseModel):
    CSI: Optional[float]
    input_used: CSIInputUsed
