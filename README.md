# 🌾 Crop Stress Index (CSI) API

The CSI API provides field-scale diagnostics of crop stress using canopy temperature modeling and environmental parameters. Designed for agronomists, researchers, and agricultural platforms, it delivers interpretable outputs and audit-friendly reproducibility.

---

## 🚀 Features

- 🔍 Modular CSI computation logic (converted from MATLAB to Python)
- 🌐 FastAPI endpoint `/csi/compute` for real-time diagnostics
- 📦 Dockerized deployment for reproducibility and portability
- 📊 Validated outputs against MATLAB reference models
- 🧾 Audit trail support with environment snapshots and logging

---

## 📌 Endpoint

### `POST /csi/compute`

**Payload:**

```json
{
  "canopy_temp": 34.2,
  "air_temp": 30.0,
  "vpd": 2.1,
  "solar_radiation": 650,
  "humidity": 45.0
}

Response:

json
{
  "csi_score": 0.78,
  "stress_level": "Moderate",
  "interpretation": "Crop is experiencing moderate thermal stress. Consider irrigation or shade strategies."
}

🧪 Validation
✅ CSI outputs match MATLAB reference values within ±0.01 tolerance
✅ Tested with real payloads from Tamil Nadu field trials
✅ Docker container runs cleanly on Ubuntu 22.04 with Conda environment

🛠️ Setup
Clone and Run
bash
git clone https://github.com/agrigenai/csi-api.git
cd csi-api
docker build -t csi-api .
docker run -p 8000:8000 csi-api
Local Dev
bash
conda env create -f environment.yml
uvicorn main:app --reload
📚 Documentation
Swagger UI: http://localhost:8000/docs

CSI Model Reference: docs/model_reference.md

Validation Logs: logs/validation.log

📈 Roadmap
🔜 Azure App Service deployment
🔜 Public API access for agronomist feedback
🔜 Enhanced logging and audit trail integration
🔜 Multilingual support for regional advisory

🧑‍🔬 About
This API is part of AgriGenAI’s mission to deliver real-world agricultural intelligence through modular, audit-friendly AI systems. Built by Sheriff, solo founder and systems architect, with a focus on reproducibility, mobile-first design, and investor-grade polish.

📄 License
MIT © 2025 AgriGenAI

📬 Contact
For collaboration, field trials, or integration inquiries:
Sheriff Founder & Systems Architect, AgriGenAI 📧 sheriff@agricps.co 🌐 https://www.agricps.co 
Prefer encrypted communication or audit trail logging? Let me know—we support secure channels for sensitive data exchange.

