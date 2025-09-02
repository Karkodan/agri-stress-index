from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api import router as csi_router  # ✅ Updated path

app = FastAPI(
    title="AgriGenAI CSI API",
    description="Flexible Crop Stress Index computation for field-scale diagnostics",
    version="1.0.0"
)

# 🌐 Enable CORS for external access (adjust origins in production)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # 🔒 Replace with specific domains in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 🔌 Register CSI router
app.include_router(csi_router, prefix="/csi")

# Optional: Health check endpoint
@app.get("/health")
def health_check():
    return {"status": "ok", "message": "CSI API is live"}
