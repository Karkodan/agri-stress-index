# services/weather_fetcher.py

import os
from dotenv import load_dotenv
import requests
from datetime import datetime

load_dotenv("/home/karkodan/AgriGenAI/.env")  # Explicit path for clarity

TOMORROW_API_KEY = os.getenv("TOMORROW_API_KEY")

def fetch_weather(lat: float, lon: float) -> dict:
    if not TOMORROW_API_KEY:
        raise ValueError("Missing TOMORROW_API_KEY in environment")

    url = "https://api.tomorrow.io/v4/weather/realtime"
    params = {
        "location": f"{lat},{lon}",
        "apikey": TOMORROW_API_KEY,
        "units": "metric"
    }
    response = requests.get(url, params=params)
    data = response.json()["data"]["values"]

    return {
        "air_temp_c": data["temperature"],
        "vpd_kpa": data.get("vaporPressureDeficit", 1.2),
        "wind_speed_m_s": data["windSpeed"],
        "solar_radiation_w_m2": data.get("solarGHI", 500),
        "hour_of_day": datetime.utcnow().hour,
        "doy": datetime.utcnow().timetuple().tm_yday
    }
