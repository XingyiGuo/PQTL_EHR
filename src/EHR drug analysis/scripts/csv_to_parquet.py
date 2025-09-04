import pandas as pd
from pathlib import Path

BASE = Path("data")

# Simple helper: convert known CSVs to Parquet with parsed date columns
def to_parquet(csv_path, date_cols=None):
    df = pd.read_csv(csv_path, parse_dates=date_cols or [])
    df.to_parquet(csv_path.with_suffix(".parquet"))

if __name__ == "__main__":
    # Flat files
    to_parquet(BASE / "pat_cancer.csv", ["cancer_start_date"])
    to_parquet(BASE / "pat_record.csv", ["min_record_date", "max_record_date"])
    to_parquet(BASE / "pat_demo.csv", ["birth_datetime"])

    # Drug exposures
    for p in (BASE / "drug_exposure_by_drug").glob("*.csv"):
        to_parquet(p, ["drug_exposure_start_datetime"])

    # Confounder detail
    for p in (BASE / "confounder_by_phecode_detail").glob("*.csv"):
        to_parquet(p, ["concept_date"])

    print("Conversion complete. Parquet files are written alongside the CSVs.")
