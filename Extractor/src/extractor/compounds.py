import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import json
from extractor.gnps import IteratedQueries
from collections import Counter
from ms2decide.Tool import Tool
from ms2decide.K import K
from ms2decide.Tanimotos import Tanimotos

def get_task_ids(task_ids_file):
    with open(task_ids_file) as task_ids_data:
        task_ids = set(json.load(task_ids_data))
    return task_ids


def get_iterated_queries(task_ids_file, dir_tasks_cached):
    task_ids = get_task_ids(task_ids_file)
    qs = IteratedQueries.from_task_ids(task_ids, dir_tasks_cached)
    return qs


def get_tanimotos_by_id(ids, inchis, score_gnps):
    tanimotos_by_id = {}
    for id in ids:
        inchi_g = inchis.loc[id, "InChI GNPS"]
        inchi_s = inchis.loc[id, "InChI Sirius"]
        inchi_i = inchis.loc[id, "InChI ISDB-LOTUS"]

        by_tool = {
            Tool.GNPS.name: inchi_g if pd.notna(inchi_g) else "*",
            Tool.SIRIUS.name: inchi_s if pd.notna(inchi_s) else "*",
            Tool.ISDB.name: inchi_i if pd.notna(inchi_i) else "*",
        }
        tanimotos = Tanimotos(by_tool)
        tanimotos.compute_tanimoto()

        if pd.isna(inchi_g):
            assert tanimotos.tgs == 0
            assert tanimotos.tgi == 0
            if pd.notna(score_gnps[id]) and not pd.isna(inchi_s):
                tanimotos.tgs = 0.7
            if pd.notna(score_gnps[id]) and not pd.isna(inchi_i):
                tanimotos.tgi = 0.7
        if pd.isna(inchi_s):
            assert tanimotos.tgs == 0
            assert tanimotos.tsi == 0
        if pd.isna(inchi_i):
            assert tanimotos.tgi == 0
            assert tanimotos.tsi == 0
        if pd.isna(inchi_s):
            tanimotos.tgs = 0.25
            tanimotos.tsi = 0.25
        tanimotos_by_id[id] = tanimotos
    return tanimotos_by_id

def get_k_series(ids, scores_df, tanimotos_by_id):
    k_by_id = {}
    for id in ids:
        similarities = {
            Tool.GNPS.name: scores_df.loc[id, Tool.GNPS.name],
            Tool.SIRIUS.name: scores_df.loc[id, Tool.SIRIUS.name],
            Tool.ISDB.name: scores_df.loc[id, Tool.ISDB.name],
        }
        tanimotos = tanimotos_by_id[id]
        k = K(similarities, tanimotos).k()
        k_by_id[id] = k
    k_df = pd.Series(k_by_id, name="K")
    return k_df

def add_ranks_columns(df, column_name, score_column):
    rank_min_column = f"Rank min {column_name}"
    rank_max_column = f"Rank max {column_name}"
    ranks_column = f"Ranks {column_name}"
    rank_column = f"Rank {column_name}"
    df[rank_min_column] = (
        df[score_column].fillna(0).rank(method="min").astype(int)
    )
    df[rank_max_column] = (
        df[score_column].fillna(0).rank(method="max").astype(int)
    )
    df[ranks_column] = df.apply(
        lambda x: (
            str(x[rank_min_column])
            if x[rank_min_column] == x[rank_max_column]
            else (str(x[rank_min_column]) + "–" + str(x[rank_max_column]))
        ),
        axis=1,
    )
    if (df[rank_min_column] == df[rank_max_column]).all():
        df.rename(columns={ranks_column: rank_column}, inplace=True)

class Compounds:
    def __init__(self, compounds: pd.DataFrame):
        self.df = compounds
    
    @classmethod
    def from_tsv(cls, compounds_file, name_id_df):
        df = pd.read_csv(compounds_file, sep="\t").join(name_id_df, on = "Chemical name").set_index("Id").sort_index()
        return cls(df)
    
    def add_precursors(self, precursors, suffix=""):
        suffix_spaced = " " + suffix if suffix else ""
        self.df["Precursor m/z" + suffix_spaced] = precursors
    
    def add_retention_times(self, retention_seconds, suffix=""):
        suffix_spaced = " " + suffix if suffix else ""
        self.df["Retention time (sec)" + suffix_spaced] = retention_seconds
    
    def add_relative_molecular_weights(self):
        self.df["Relative molecular weight"] = self.df["InChI"].apply(lambda i: Descriptors.MolWt(Chem.inchi.MolFromInchi(i)) if pd.notna(i) else None)

    def add_gnps_standard(self):
        self.df["GNPS standard binary novelty indicator"] = self.df["Analog Score GNPS; peaks ≥ 6; Δ mass ≤ 0.02"] < 0.7
    
    def add_diffs(self):
        self.df["Precursor m/z GNPS − relative molecular weight"] = (
        self.df["Precursor m/z GNPS"] - self.df["Relative molecular weight"]
    )
        
    def quantification_table_minutes(self, precursors = None, retention_seconds = None):
        precursors_series = precursors if precursors is not None else (self.df["Precursor m/z"] if "Precursor m/z" in self.df.columns else pd.Series())
        retention_seconds_series = retention_seconds if retention_seconds is not None else (self.df["Retention time (sec)"] if "Retention time (sec)" in self.df.columns else pd.Series())
        qt = pd.DataFrame(index=self.df.index).rename_axis("row ID")
        qt["row m/z"] = precursors_series
        qt["row retention time"] = retention_seconds_series / 60.0
        qt["1.mzXML Peak area"] = 0
        return qt
    
    def join_iterated_queries(self, task_ids_file, dir_tasks_cached):
        qs = get_iterated_queries(task_ids_file, dir_tasks_cached)
        all_df = qs.all_df()
        self.df = self.df.join(all_df)
    
    def join_sirius_data(self, sirius_tsv):
        sirius_df = pd.read_csv(sirius_tsv, sep="\t").set_index("mappingFeatureId")
        sirius_df["Score Sirius approximate"] = sirius_df["ConfidenceScoreApproximate"].replace({float("-inf"): 0})
        sirius_df["Adduct Sirius"] = sirius_df["adduct"].map(lambda s: s.replace(" ", ""))
        sirius_df = sirius_df.rename(columns={"InChI": "InChI Sirius"}).loc[
            :, ["InChI Sirius", "Score Sirius approximate", "Adduct Sirius"]
        ]
        self.df = self.df.join(sirius_df)

    def add_adduct_summary(self):
        self.df["Adduct GNPS and Sirius"] = self.df.apply(
            lambda r: str(dict(Counter(r.filter(like = "Adduct ")))), axis=1
        )

    def join_isdb_data(self, isdb_tsv):
        isdb_df = (
            pd.read_csv(isdb_tsv, sep="\t")
            .set_index("Id")
            .rename(columns={"InChI": "InChI ISDB-LOTUS", "Score": "Score ISDB-LOTUS"})
        )
        self.df = self.df.join(isdb_df)

    def join_fermo(self, fermo_csv):
        more_df = (
            pd.read_csv(fermo_csv)
            .set_index("feature_id")
            .rename(columns={"novelty_score": "Score FERMO"})
        )
        self.df = self.df.join(more_df)
        self.df["Score compl FERMO"] = 1 - self.df["Score FERMO"]

    def join_k_tuples(self):
        self._join_k_tuples("Score GNPS iterated discounted", "Standard InChI GNPS iterated", suffix="")
        self._join_k_tuples("Analog Score GNPS; peaks ≥ 6; Δ mass ≤ 0.02", "Analog Standard InChI GNPS; peaks ≥ 6; Δ mass ≤ 0.02", suffix="without iterated")
    
    def _join_k_tuples(self, score_gnps_column_name, inchi_gnps_column_name, suffix):
        suffix_spaced = " " + suffix if suffix else ""
        self.df["cg" + suffix_spaced] = self.df[score_gnps_column_name].fillna(0)
        self.df["cs" + suffix_spaced] = self.df["Score Sirius approximate"].fillna(0.5)
        assert self.df["Score ISDB-LOTUS"].notna().all()
        self.df["ci" + suffix_spaced] = self.df["Score ISDB-LOTUS"]

        inchis_df = pd.DataFrame(self.df.loc[:, [inchi_gnps_column_name, "InChI Sirius", "InChI ISDB-LOTUS"]].rename(columns={inchi_gnps_column_name: "InChI GNPS"}))
        tanimotos_by_id = get_tanimotos_by_id(self.df.index, inchis_df, self.df[score_gnps_column_name])
        tanimoto_values_by_id = {i: {"tgs": t.tgs, "tgi": t.tgi, "tsi": t.tsi} for (i, t) in tanimotos_by_id.items()}
        tanimotos_df = pd.DataFrame.from_dict(tanimoto_values_by_id, orient="index")
        self.df = self.df.join(tanimotos_df, rsuffix=suffix_spaced)

        k_df = get_k_series(self.df.index, self.df.loc[:, ["cg" + suffix_spaced, "cs" + suffix_spaced, "ci" + suffix_spaced]].rename({"cg" + suffix_spaced: Tool.GNPS.name, "cs" + suffix_spaced: Tool.SIRIUS.name, "ci" + suffix_spaced: Tool.ISDB.name}, axis=1), tanimotos_by_id)
        self.df = self.df.join(k_df, rsuffix=suffix_spaced)

    def add_ranks(self):
        add_ranks_columns(
            self.df,
            "GNPS original",
            "Analog Score GNPS; peaks ≥ 6; Δ mass ≤ 0.02",
        )
        add_ranks_columns(
            self.df,
            "GNPS standard ranking",
            "GNPS standard binary novelty indicator",
        )
        add_ranks_columns(
            self.df,
            "GNPS iterated",
            "Score GNPS iterated discounted",
        )
        add_ranks_columns(self.df, "Sirius", "Score Sirius approximate")
        add_ranks_columns(
            self.df, "ISDB-LOTUS", "Score ISDB-LOTUS"
        )
        add_ranks_columns(
            self.df, "FERMO", "Score compl FERMO"
        )
        add_ranks_columns(self.df, "K", "K")
        add_ranks_columns(self.df, "K without iterated", "K without iterated")

    def get_counts(self):
        # selected_columns = self.df.columns.map(lambda s: (s.startswith("Standard InChI GNPS; ")) or (s == "Id"))
        standards = self.df.filter(like="Standard InChI GNPS")
        # columns = standards.columns
        # es = columns.map(lambda s: s.endswith(" exact"))
        # ces = columns[es]
        # nes = columns.map(lambda s: not s.endswith(" exact"))
        # cnes = columns[nes]
        # list(ces) + list(cnes)
        counts = standards.agg("count")
        counts.index.name = "Column"
        counts.name = "Count"
        return counts
    
    def prefix_semantic_ids(self):
        duplicated_precursors_indices = self.df.set_index("Precursor m/z").index.duplicated(keep=False)
        duplicated_precursors = self.df.set_index("Precursor m/z").index[duplicated_precursors_indices]
        sids = self.df.apply(lambda r: (str(r["Precursor m/z"]) + ";" + str(r["Retention time (seconds)"])) if r["Precursor m/z"] in duplicated_precursors else str(r["Precursor m/z"]), axis=1)
        assert sids.duplicated().sum() == 0
        assert (sids.dtype.name == "object")
        assert (sids.map(type).eq(str).all())
        self.df.insert(0, "Semantic id", sids)
