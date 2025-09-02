# main.py

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from app.api import router as csi_router  # ✅ Modular CSI router

# 🧠 GPU detection
try:
    import tensorflow as tf
    gpu_check_enabled = True
except ImportError:
    gpu_check_enabled = False

# 🚀 FastAPI app metadata
app = FastAPI(
    title="AgriGenAI CSI API",
    description="""
Flexible Crop Stress Index (CSI) computation for field-scale diagnostics.

**Version:** 1.0.0  
**Maintainer:** Sheriff, Solo Founder of AgriGenAI  
**Docs:** [GitHub](https://github.com/agrigenai/csi-api)  
**Contact:** sheriff@agrigen.ai
""",
    version="1.0.0",
    contact={
        "name": "Sheriff",
        "email": "sheriff@agrigen.ai",
        "url": "https://agrigen.ai"
    },
    license_info={
        "name": "MIT",
        "url": "https://opensource.org/licenses/MIT"
    }
)

# 🌐 CORS setup — restrict origins in production
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # 🔒 Replace with specific domains like https://agrigenai.com in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 📁 Serve static files (for favicon)
app.mount("/static", StaticFiles(directory="app/static"), name="static")

# 🌿 Serve favicon
@app.get("/favicon.ico", include_in_schema=False)
def favicon():
    return FileResponse("app/static/favicon.ico")

# 🔌 Register CSI router under /csi
app.include_router(csi_router, prefix="/csi", tags=["CSI"])

# 🧪 Log GPU status on startup
@app.on_event("startup")
def log_gpu_status():
    if gpu_check_enabled:
        gpus = tf.config.list_physical_devices('GPU')
        print(f"[CSI] GPU devices detected: {gpus}")
    else:
        print("[CSI] TensorFlow not available — running in CPU-only mode")

# 🩺 Health check endpoint
@app.get("/health", tags=["System"])
def health_check():
    if gpu_check_enabled:
        gpus = tf.config.list_physical_devices('GPU')
        return {
            "gpu_available": len(gpus) > 0,
            "gpu_count": len(gpus),
            "device_names": [gpu.name for gpu in gpus]
        }
    else:
        return {"gpu_available": False, "gpu_count": 0, "device_names": []}

# 🏠 Root endpoint for sanity check
@app.get("/", tags=["System"])
def root():
    return {
        "status": "AgriGenAI CSI API is running",
        "version": "1.0.0",
        "docs": "/docs",
        "health": "/health",
        "csi": "/csi"
    }
