from models.radiation import compute_canopy_radiation

result = compute_canopy_radiation(
    LAT=13.08, R_inc=800, Ta=30, Tl=32, DOY=241, tm=12.0, cloud_frac=0.3
)

print(f"☀️ Q_short: {result[0]:.2f} W/m²")
print(f"🌱 Qp: {result[1]:.2f} μmol/m²/s")
print(f"🌙 Q_long: {result[2]:.2f} W/m²")
