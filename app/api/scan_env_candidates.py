import pandas as pd
from pathlib import Path

DATA_DIR = Path("/home/karkodan/agri-stress-index/data")
CANDIDATE_COLUMNS = ["timestamp", "date_time", "datetime", "recorded_at"]
CATCHMENT_COLUMNS = ["catchment_id", "location_id", "site_id"]

def scan_csv_headers():
    for file in DATA_DIR.rglob("*.csv"):
        try:
            header = pd.read_csv(file, nrows=0).columns.str.lower()
            has_time = any(col in header for col in CANDIDATE_COLUMNS)
            has_catchment = any(col in header for col in CATCHMENT_COLUMNS)

            if has_time and has_catchment:
                print(f"✅ Candidate: {file}")
                print(f"   Columns: {', '.join(header)}\n")

        except Exception as e:
            print(f"⚠️ Skipped: {file} — {e}")

if __name__ == "__main__":
    scan_csv_headers()
