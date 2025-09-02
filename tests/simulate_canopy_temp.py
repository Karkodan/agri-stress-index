import requests

def simulate_canopy_temp(moisture):
    payload = {
        "air_temp_c": 32,
        "relative_humidity_pct": 55,
        "solar_radiation_w_m2": 850,
        "wind_speed_m_s": 1.5,
        "soil_moisture_pct": moisture
    }
    response = requests.post("http://localhost:8000/canopy-temp", json=payload)
    result = response.json()["canopy_temperature"]
    print(f"Moisture: {moisture}% → Canopy Temp: {result:.2f}°C")

if __name__ == "__main__":
    print("🧪 Simulating canopy temperature across moisture levels...")
    for moisture in range(10, 100, 10):
        simulate_canopy_temp(moisture)
