from extractor import extractor
import pandas as pd

def test_status():
    status_df = extractor.read_status()
    assert len(status_df) >= 1
    assert status_df.at["PMH22", "year"] == 2022
    assert str(status_df.at["PMH22", "year"]) == "2022"
    assert pd.isna(status_df.at["VH??", "year"])