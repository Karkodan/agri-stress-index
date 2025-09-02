import os
import pandas as pd

# 📍 Paths
RAW_DATA_PATH = os.path.expanduser("~/agri-stress-index/data/Plant_Parameters.csv")
PROCESSED_DATA_PATH = os.path.expanduser("~/agri-stress-index/data/Plant_Parameters.processed.csv")

def load_raw_data(path):
    """Load raw canopy temperature training data."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"❌ Dataset not found at {path}")
    df = pd.read_csv(path)
    print(f"✅ Loaded raw data: {df.shape[0]} rows, {df.shape[1]} columns")
    return df

def encode_features(df):
    """Encode categorical features for ML model."""
    categorical_cols = [col for col in ["crop_type", "growth_stage"] if col in df.columns]
    if not categorical_cols:
        print("⚠️ No categorical columns found to encode.")
        return df
    df_encoded = pd.get_dummies(df, columns=categorical_cols)
    print(f"✅ Encoded features: {df_encoded.shape[1]} columns")
    return df_encoded

def save_processed_data(df, path):
    """Save processed dataset for training."""
    df.to_csv(path, index=False)
    print(f"✅ Saved processed data to {path}")

def main():
    print("🚀 Starting sample data load and preprocessing...")
    df_raw = load_raw_data(RAW_DATA_PATH)
    df_processed = encode_features(df_raw)
    save_processed_data(df_processed, PROCESSED_DATA_PATH)
    print("🏁 Data preprocessing complete.")

if __name__ == "__main__":
    main()
