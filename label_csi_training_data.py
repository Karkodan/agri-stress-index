import pandas as pd
from models.virtual_sensor import compute_crop_stress_index
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import joblib

# Load your processed plant parameters
df = pd.read_csv("data/Plant_Parameters.processed.csv")

# Rename columns to match expected input names
df = df.rename(columns={
    "AirTemp": "air_temp_c",
    "Humidity": "relative_humidity_pct",
    "Radiation": "solar_radiation_w_m2",
    "WindSpeed": "wind_speed_m_s",
    "SoilMoisture": "soil_moisture_pct",
    "Moisture": "soil_moisture_pct",  # fallback
    "Temperature": "air_temp_c"       # fallback
})

# Required columns for CSI computation
required_cols = [
    "air_temp_c",
    "relative_humidity_pct",
    "solar_radiation_w_m2",
    "wind_speed_m_s",
    "soil_moisture_pct"
]

# Fill missing columns with None
for col in required_cols:
    if col not in df.columns:
        print(f"⚠️ Missing column: {col} — filling with None")
        df[col] = None

# Drop rows with missing values
df = df.dropna(subset=required_cols)

# Compute CSI stress score for each row
def compute_score(row):
    result = compute_crop_stress_index(
        row["air_temp_c"],
        row["relative_humidity_pct"],
        row["solar_radiation_w_m2"],
        row["wind_speed_m_s"],
        row["soil_moisture_pct"]
    )
    return result["stress_score"]

df["stress_score"] = df.apply(compute_score, axis=1)

# Save labeled data
df.to_csv("data/training_data.with_csi.csv", index=False)
print("✅ CSI-labeled training data saved to: data/training_data.with_csi.csv")

# -------------------------------
# 🧪 Train a model on CSI scores
# -------------------------------

# Select features and target
X = df[required_cols]
y = df["stress_score"]

# Split into train/test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train a Random Forest model
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Evaluate
y_pred = model.predict(X_test)
rmse = mean_squared_error(y_test, y_pred, squared=False)
print(f"📊 Model RMSE: {rmse:.2f}")

# Save model
joblib.dump(model, "models/csi_model.pkl")
print("✅ CSI model saved to: models/csi_model.pkl")
