from fastapi.testclient import TestClient
from app.main import app  # Adjust if your FastAPI entry point differs

client = TestClient(app)

# 🎯 Base payload
base_payload = {
    "lat": 12.9716,
    "lon": 77.5946,
    "timestamp": "2025-09-01T08:00:00Z",
    "canopy_temp": 35.2,
    "ambient_temp": 30.0,
    "soil_moisture": 0.25
}

def test_csi_api_success():
    response = client.post("/csi/compute", json=base_payload)
    assert response.status_code == 200
    data = response.json()
    assert "CSI" in data
    assert isinstance(data["CSI"], float)
    assert data["stress_level"] in ["No stress", "Mild stress", "Moderate stress", "Severe stress"]
    assert "recommendation" in data
    assert isinstance(data["scale"], dict)

def test_invalid_timestamp():
    payload = base_payload.copy()
    payload["timestamp"] = "not-a-date"
    response = client.post("/csi/compute", json=payload)
    assert response.status_code == 422

def test_missing_field():
    payload = base_payload.copy()
    del payload["soil_moisture"]
    response = client.post("/csi/compute", json=payload)
    assert response.status_code == 422

def test_extreme_values_severe_stress():
    payload = base_payload.copy()
    payload["canopy_temp"] = 60.0
    payload["soil_moisture"] = 0.01
    response = client.post("/csi/compute", json=payload)
    assert response.status_code == 200
    assert response.json()["stress_level"] == "Severe stress"

def test_batch_payloads():
    batch = [
        {"canopy_temp": 30.0, "ambient_temp": 28.0, "soil_moisture": 0.3},
        {"canopy_temp": 35.0, "ambient_temp": 30.0, "soil_moisture": 0.2},
        {"canopy_temp": 45.0, "ambient_temp": 30.0, "soil_moisture": 0.05}
    ]
    for i, case in enumerate(batch):
        payload = base_payload.copy()
        payload.update(case)
        response = client.post("/csi/compute", json=payload)
        assert response.status_code == 200, f"Batch {i} failed"
        assert "CSI" in response.json()
