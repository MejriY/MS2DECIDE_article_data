from pathlib import Path
import matchms
import pandas as pd
from functools import cache

class Mgf:
    def __init__(self, file: Path):
        spectra = list(matchms.importing.load_from_mgf(str(file)))
        metadatas = [s.metadata for s in spectra]
        df = pd.DataFrame(metadatas).astype({"feature_id": int, "scans": int})
        df_2_orig_features = set(df["feature_id"])
        if("spectype" in df.columns):
            corr = df["spectype"] == "CORRELATED MS"
            l1 = (df["ms_level"] == "1")
            to_remove = df[corr & l1].index
            df = df.drop(to_remove)
            df_2_remaining_features = set(df["feature_id"])
            assert df_2_orig_features == df_2_remaining_features, (df_2_orig_features, df_2_remaining_features)
            assert df["spectype"].isna().all()
            del df["spectype"]
        if("file_name" in df.columns):
            assert df["file_name"].isna().all()
            del df["file_name"]
        if("num_peaks" in df.columns):
            del df["num_peaks"]
        nb_rep = df.groupby("feature_id").nunique()
        features_with_repeated_ms_level = nb_rep[nb_rep["ms_level"] > 1].index
        l2 = (df["ms_level"] == "2")
        ft_rep_lines = df["feature_id"].isin(set(features_with_repeated_ms_level))
        to_remove = df[l2 & ft_rep_lines].index
        df = df.drop(to_remove)
        max_repetitions = df.groupby("feature_id").nunique().max()
        assert max_repetitions.max() == 1
        df = df.drop_duplicates()
        df_2_remaining_features = set(df["feature_id"])
        assert df_2_orig_features == df_2_remaining_features
        assert (df["feature_id"] == df["scans"]).all()
        self.df = df.set_index("feature_id").drop(columns=["scans"]).rename(columns={"charge": "Charge", "ms_level": "MS level"})
