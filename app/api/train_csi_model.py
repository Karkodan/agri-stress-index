import os
import pandas as pd
import numpy as np
import joblib
import ast
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score

# Paths
data_root = "/home/karkodan/agri-stress-index/data"
enriched_dir = os.path.join(data_root, "camels_india/enriched_with_csi")
legacy_files = [
    os.path.join(data_root, "training_data.with_csi.csv"),
    os.path.join(data_root, "camels_india_03001_with_csi.csv")
]
model_dir = "/home/karkodan/agri-stress-index/models"
os.makedirs(model_dir, exist_ok=True)

# Manifest log
manifest = []

# Helper to extract stress_score from stringified dict
def extract_score(val):
    try:
        parsed = ast.literal_eval(val)
        return parsed.get('stress_score', None)
    except Exception:
        return None

# Aggregate CSI data
all_dfs = []

def process_file(path, source):
    try:
        df = pd.read_csv(path)
        total_rows = len(df)
        if 'csi' not in df.columns:
            manifest.append({
                "file": os.path.basename(path),
                "source": source,
                "total_rows": total_rows,
                "valid_csi_rows": 0,
                "status": "❌ Missing 'csi'"
            })
            return

        df['csi_score'] = df['csi'].apply(extract_score)
        valid_rows = df['csi_score'].dropna()
        all_dfs.append(valid_rows.to_frame())

        manifest.append({
            "file": os.path.basename(path),
            "source": source,
            "total_rows": total_rows,
            "valid_csi_rows": len(valid_rows),
            "status": "✅ OK" if len(valid_rows) > 0 else "❌ No valid CSI"
        })

    except Exception as e:
        manifest.append({
            "file": os.path.basename(path),
            "source": source,
            "total_rows": 0,
            "valid_csi_rows": 0,
            "status": f"⚠️ Error: {str(e)[:60]}"
        })

# Process enriched files
for file in os.listdir(enriched_dir):
    if file.endswith(".csv"):
        process_file(os.path.join(enriched_dir, file), "enriched")

# Process legacy files
for file in legacy_files:
    if os.path.exists(file):
        process_file(file, "legacy")

# Save manifest
manifest_df = pd.DataFrame(manifest)
manifest_df.to_csv(os.path.join(model_dir, "csi_manifest.csv"), index=False)

# Combine all valid CSI scores
combined_df = pd.concat(all_dfs, ignore_index=True)
print(f"📊 Total CSI samples aggregated: {len(combined_df)}")

# Drop missing
combined_df = combined_df.dropna()

# Simulate environmental features (replace with real ones later)
np.random.seed(42)
X = pd.DataFrame({
    "temperature": np.random.uniform(25, 40, size=len(combined_df)),
    "humidity": np.random.uniform(30, 90, size=len(combined_df)),
    "wind_speed": np.random.uniform(0.5, 5.0, size=len(combined_df)),
    "solar_radiation": np.random.uniform(400, 800, size=len(combined_df)),
    "soil_moisture": np.random.uniform(0.1, 0.4, size=len(combined_df))
})
y = combined_df['csi_score']

# Normalize features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

# Train model
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Evaluate
y_pred = model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"✅ Model trained successfully")
print(f"📉 MSE: {mse:.4f}")
print(f"📈 R² Score: {r2:.4f}")

# Save model and scaler
joblib.dump(model, os.path.join(model_dir, "csi_model.pkl"))
joblib.dump(scaler, os.path.join(model_dir, "csi_scaler.pkl"))

# Save training log
with open(os.path.join(model_dir, "training_log.txt"), "w") as log:
    log.write(f"Samples: {len(combined_df)}\n")
    log.write(f"MSE: {mse:.4f}\n")
    log.write(f"R²: {r2:.4f}\n")
    log.write("Model: RandomForestRegressor(n_estimators=100)\n")
    log.write("Scaler: StandardScaler()\n")

print("💾 Model, scaler, manifest, and training log saved to /models/")
