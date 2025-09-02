import pandas as pd
import os
from tqdm import tqdm
from pathlib import Path

CSI_DIR = Path("/data/csi/")
PARAMS_FILE = Path("/home/karkodan/agri-stress-index/data/Plant_Parameters.processed.csv")
OUTPUT_FILE = Path("/data/aggregated/csi_training_data.csv")

def load_environmental_data():
    env_df = pd.read_csv(PARAMS_FILE, parse_dates=["timestamp"])
    env_df.dropna(subset=["catchment_id", "timestamp"], inplace=True)
    env_df.sort_values(["catchment_id", "timestamp"], inplace=True)
    return env_df

def aggregate_csi_with_features(env_df):
    all_rows = []
    skipped = 0

    for file in tqdm(sorted(CSI_DIR.glob("*.csv"))):
        csi_df = pd.read_csv(file, parse_dates=["timestamp"])
        csi_df.dropna(subset=["catchment_id", "timestamp", "csi"], inplace=True)

        # Merge with nearest timestamp per catchment
        merged = pd.merge_asof(
            csi_df.sort_values("timestamp"),
            env_df,
            on="timestamp",
            by="catchment_id",
            direction="nearest",
            tolerance=pd.Timedelta("1H")  # adjustable
        )

        if merged.empty:
            skipped += len(csi_df)
            continue

        all_rows.append(merged)

    combined_df = pd.concat(all_rows, ignore_index=True)
    print(f"✅ Aggregated samples: {len(combined_df)}")
    print(f"⚠️ Skipped samples due to missing joins: {skipped}")
    return combined_df

def main():
    env_df = load_environmental_data()
    combined_df = aggregate_csi_with_features(env_df)
    combined_df.to_csv(OUTPUT_FILE, index=False)
    print(f"💾 Saved aggregated training data to: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
