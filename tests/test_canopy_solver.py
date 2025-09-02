from models.canopy_solver import solve_canopy_temperature

env = {
    "LAT": 13.08,
    "DOY": 241,
    "tm": 12.0,
    "R_inc": 800,
    "Tair": 30.0,
    "vpd": 2.0,
    "pressure": 101.3,
    "cloud_frac": 0.3
}

result = solve_canopy_temperature(env)
print(f"🌿 Leaf temperature: {result['leaf_temperature']:.2f} °C")
print(f"☀️ Q_short: {result['Q_short']:.2f} W/m²")
print(f"🌱 Qp: {result['Qp']:.2f} μmol/m²/s")
print(f"🌙 Q_long: {result['Q_long']:.2f} W/m²")
