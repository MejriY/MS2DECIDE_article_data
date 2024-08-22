import os
import pandas as pd
import json
from extractor.gnps import GnpsCacher
from extractor.gnps import GnpsAnnotations
from extractor.gnps import GnpsIteratedNp
from extractor.gnps import GnpsParametersFile
from extractor.gnps import GnpsInchiScore
from extractor.gnps import GnpsTasks
from extractor.mgfs import MgfFiles
from extractor.support import add_ranks_columns
from rdkit import Chem
from rdkit.Chem import Descriptors
import matchms
from ms2decide.MgfInstance import MgfInstance
from ms2decide.Tool import Tool
from ms2decide.IsdbAnnotation import get_cfm_annotation
from ms2decide.K import K
from ms2decide.Tanimotos import Tanimotos
from shutil import rmtree
from extractor.pleiocarpa.datadirs import *
from collections import Counter
import extractor.support as support

def compute_isdb():
    GENERATED_DIR_ISDB.mkdir(parents=True, exist_ok=True)

    mgf = MgfInstance(INPUT_DIR / "2 - MZmine" / "Pleiocarpa.mgf")
    l = get_cfm_annotation(mgf)
    inchis = pd.Series({k: m.inchi for (k, m) in l.items()})
    scores = pd.Series({k: m.score for (k, m) in l.items()})
    df = pd.DataFrame({"InChI": inchis, "Score": scores})
    df.index.name = "Id"
    df.to_csv(GENERATED_DIR_ISDB / "ISDB-LOTUS annotations.tsv", sep="\t")

def transform_isdb():
    df = pd.read_csv(GENERATED_DIR_ISDB / "ISDB-LOTUS annotations.tsv", sep="\t").set_index("Id")
    df.replace({"#": None}).to_csv(GENERATED_DIR_ISDB / "ISDB-LOTUS annotations.tsv", sep="\t")

def generate_summary():
    GENERATED_DIR_TABLES.mkdir(parents=True, exist_ok=True)

    mgf_file_1 = INPUT_DIR / "2 - MZmine" / "Pleiocarpa.mgf"
    spectra_from_mgf_file_1 = list(matchms.importing.load_from_mgf(str(mgf_file_1)))
    ms_1 = [s.metadata for s in spectra_from_mgf_file_1]
    df_1 = pd.DataFrame(ms_1).set_index("feature_id")
    assert (df_1.index == df_1["scans"]).all()
    df_1 = df_1.drop(columns=["scans"])
    df_1 = df_1.rename(lambda x: x + " raw", axis=1).rename(columns={"charge raw": "Charge raw", "ms_level raw": "MS level raw"})

    mgf_file_2 = INPUT_DIR / "2 - MZmine" / "Pleiocarpa_sirius.mgf"
    spectra_from_mgf_file_2 = list(matchms.importing.load_from_mgf(str(mgf_file_2)))
    ms_2 = [s.metadata for s in spectra_from_mgf_file_2]
    df_2 = pd.DataFrame(ms_2)
    df_2_orig_features = set(df_2["feature_id"].astype(int))
    corr = df_2["spectype"] == "CORRELATED MS"
    l1 = (df_2["ms_level"] == "1")
    to_remove = df_2[corr & l1].index
    df_2 = df_2.drop(to_remove)
    df_2_remaining_features = set(df_2["feature_id"].astype(int))
    assert df_2_orig_features == df_2_remaining_features
    assert df_2["spectype"].isna().all()
    del df_2["spectype"]
    assert df_2["file_name"].isna().all()
    del df_2["file_name"]
    del df_2["num_peaks"]
    nb_rep = df_2.groupby("feature_id").nunique()
    features_with_repeated_ms_level = nb_rep[nb_rep["ms_level"] > 1].index
    l2 = (df_2["ms_level"] == "2")
    ft_rep_lines = df_2["feature_id"].isin(set(features_with_repeated_ms_level))
    to_remove = df_2[l2 & ft_rep_lines].index
    df_2 = df_2.drop(to_remove)
    max_repetitions = df_2.groupby("feature_id").nunique().max()
    assert max_repetitions.max() == 1
    df_2 = df_2.drop_duplicates()
    df_2_remaining_features = set(df_2["feature_id"].astype(int))
    assert df_2_orig_features == df_2_remaining_features
    assert (df_2["feature_id"] == df_2["scans"]).all()
    df_2 = df_2.drop(columns=["scans"])
    df_2 = df_2.set_index("feature_id").rename(lambda x: x + " Sirius", axis=1).rename(columns={"charge Sirius": "Charge Sirius", "ms_level Sirius": "MS level Sirius"})

    assert len(df_1) == len(df_2)
    compounds_joined = pd.concat([df_1, df_2], axis=1)
    compounds_joined.index = compounds_joined.index.astype(int)
    assert (abs(compounds_joined["retention_time raw"] - compounds_joined["retention_time Sirius"]) <= 0.005).all()
    compounds_joined = compounds_joined.drop(columns=["retention_time Sirius"]).rename(columns={"retention_time raw": "Retention time"})
    assert (abs(compounds_joined["precursor_mz raw"] - compounds_joined["precursor_mz Sirius"]) <= 0.0005).all()
    compounds_joined = compounds_joined.drop(columns=["precursor_mz Sirius"]).rename(columns={"precursor_mz raw": "Precursor m/z"})
    assert len(compounds_joined.columns) == 6
    compounds_joined = compounds_joined.reindex(columns = ["Precursor m/z", "Retention time", "Charge raw", "Charge Sirius", "MS level raw", "MS level Sirius"])
    compounds_joined.index.name = "Id"

    task_ids_file = INPUT_DIR_GNPS_TASKS / "Gnps task ids.json"
    with open(task_ids_file) as task_ids_data:
        task_ids = set(json.load(task_ids_data))

    ts = GnpsTasks(GENERATED_DIR_TABLES / "Fetched/", task_ids)
    ts.load()
    compounds_joined = compounds_joined.join(ts.all_matches())

    sirius_df = pd.read_csv(INPUT_DIR_SIRIUS / "structure_identifications.tsv", sep="\t").set_index("mappingFeatureId")
    sirius_df["Score Sirius"] = sirius_df["ConfidenceScoreExact"].replace({float("-inf"): 0})
    sirius_df["Adduct Sirius"] = sirius_df["adduct"].map(lambda s: s.replace(" ", ""))
    sirius_df = sirius_df.rename(columns={"InChI": "InChI Sirius"}).loc[:, ["InChI Sirius", "Score Sirius", "Adduct Sirius"]]
    compounds_joined = compounds_joined.join(sirius_df)

    compounds_joined["Adduct GNPS and Sirius"] = compounds_joined.apply(lambda r: str(dict(Counter(r[r.index.map(lambda c: c.startswith("Adduct "))]))), axis=1)

    isdb_df = (
        pd.read_csv(GENERATED_DIR_ISDB / "ISDB-LOTUS annotations.tsv", sep="\t")
        .set_index("Id")
        .rename(columns={"InChI": "InChI ISDB-LOTUS", "Score": "Score ISDB-LOTUS"})
    )
    compounds_joined = compounds_joined.join(isdb_df)
    
    ids = compounds_joined.index

    compounds_joined["cg"] = compounds_joined["Score GNPS iterated discounted"].fillna(0)
    compounds_joined["cs"] = compounds_joined["Score Sirius"].fillna(0.5)
    compounds_joined["ci"] = compounds_joined["Score ISDB-LOTUS"].fillna(0)

    tanimotos_by_id = {}
    tanimoto_values_by_id = {}
    for id in ids:
        inchi_g = compounds_joined.loc[id, "Standard InChI GNPS iterated"]
        inchi_s = compounds_joined.loc[id, "InChI Sirius"]
        inchi_i = compounds_joined.loc[id, "InChI ISDB-LOTUS"]

        by_tool = {
            Tool.GNPS.name: inchi_g if not pd.isna(inchi_g) else "*",
            Tool.SIRIUS.name: inchi_s if not pd.isna(inchi_s) else "*",
            Tool.ISDB.name: inchi_i if not pd.isna(inchi_i) else "*",
        }
        tanimotos = Tanimotos(by_tool)
        (tgs, tgi, tsi) = tanimotos.compute_tanimoto()

        if pd.isna(inchi_g):
            assert tgs == 0
            assert tgi == 0
        if pd.isna(inchi_s):
            assert tgs == 0
            assert tsi == 0
        if pd.isna(inchi_i):
            assert tgi == 0
            assert tsi == 0
        if pd.isna(inchi_s):
            tgs = 0.25
            tsi = 0.25
        tanimotos_by_id[id] = tanimotos
        tanimoto_values_by_id[id] = {"tgs": tgs, "tgi": tgi, "tsi": tsi}

    tanimotos_df = pd.DataFrame.from_dict(tanimoto_values_by_id, orient="index")
    compounds_joined = compounds_joined.join(tanimotos_df)

    k_by_id = {}
    for id in ids:
        similarities = {
            Tool.GNPS.name: compounds_joined.loc[id, "cg"],
            Tool.SIRIUS.name: compounds_joined.loc[id, "cs"],
            Tool.ISDB.name: compounds_joined.loc[id, "ci"],
        }
        tanimotos = tanimotos_by_id[id]
        k = K(similarities, tanimotos).k()
        k_by_id[id] = k
    k_df = pd.Series(k_by_id, name="K")
    compounds_joined = compounds_joined.join(k_df)

    add_ranks_columns(
        compounds_joined,
        "Rank min GNPS original",
        "Rank max GNPS original",
        "Ranks GNPS original",
        "Score GNPS; peaks ≥ 6; Δ mass ≤ 0.02",
    )
    add_ranks_columns(
        compounds_joined,
        "Rank min GNPS iterated",
        "Rank max GNPS iterated",
        "Ranks GNPS iterated",
        "Score GNPS iterated discounted",
    )
    add_ranks_columns(compounds_joined, "Rank min Sirius", "Rank max Sirius", "Ranks Sirius", "Score Sirius")
    add_ranks_columns(
        compounds_joined, "Rank min ISDB-LOTUS", "Rank max ISDB-LOTUS", "Ranks ISDB-LOTUS", "Score ISDB-LOTUS"
    )
    add_ranks_columns(compounds_joined, "Rank min K", "Rank max K", "Ranks K", "K")

    duplicated_precursors_indices = compounds_joined.set_index("Precursor m/z").index.duplicated(keep=False)
    duplicated_precursors = compounds_joined.set_index("Precursor m/z").index[duplicated_precursors_indices]
    sids = compounds_joined.apply(lambda r: (str(r["Precursor m/z"]) + ";" + str(r["Retention time"])) if r["Precursor m/z"] in duplicated_precursors else r["Precursor m/z"], axis=1)
    compounds_joined.insert(0, "Semantic id", sids)
    assert compounds_joined["Semantic id"].duplicated().sum() == 0

    compounds_joined.to_csv(GENERATED_DIR_TABLES / "Compounds joined.tsv", sep="\t")

    counts = compounds_joined.filter(like="Standard InChI GNPS").agg("count")
    counts.index.name = "Column"
    counts.name = "Count"
    counts.to_csv(GENERATED_DIR_TABLES / "Counts.tsv", sep="\t")

    by_k = compounds_joined.sort_values("Rank min K").loc[:, ["cg", "cs", "ci", "tgs", "tgi", "tsi", "K", "Ranks K"]]
    by_k.to_csv(GENERATED_DIR_TABLES / "Compounds by K.tsv", sep="\t")

    # y_df = pd.read_csv(INPUT_DIR / "Pleiocarpa_annotation.tsv", sep="\t").rename(columns = {"ID": "Id"}).set_index("Id").rename(columns= lambda x: x + " Yassine")
    # compounds_joined = compounds_joined.join(y_df)
    # compounds_joined["K diff"] = (compounds_joined["K"] - compounds_joined["K Yassine"]).round(4)
    # compounds_joined.to_csv(GENERATED_DIR_TABLES / "Compounds with Yassine.tsv", sep="\t")

def generate_article_data():
    GENERATED_DIR_ARTICLE.mkdir(parents=True, exist_ok=True)

    # os.chdir("Data/Extractor/")
    # unreported_df = support.unreported_ones()
    # rounded_unreported_masses = unreported_df["Precursor m/z"].map(lambda p: round(p))

    compounds = pd.read_csv(GENERATED_DIR_TABLES / "Compounds joined.tsv", sep="\t").set_index("Id")
    # compounds["Unreported"] = compounds["Precursor m/z"].map(lambda p: round(p) in rounded_unreported_masses)
    by_k = (
        compounds
        .sort_values("Rank min K").loc[:, ["Semantic id", "Adduct GNPS and Sirius", "cg", "cs", "ci", "tgs", "tgi", "tsi", "K", "Ranks K"]]
    )
    compounds.rename(columns=lambda x: x.replace(" ", "")).to_csv(GENERATED_DIR_ARTICLE / "Compounds.csv")
    by_k.rename(columns=lambda x: x.replace(" ", "")).replace({",": ";"}, regex=True).replace({"{": ""}, regex=True).replace({"}": ""}, regex=True).replace({"'": ""}, regex=True).to_csv(GENERATED_DIR_ARTICLE / "K.csv")
    
# compute_isdb()
# transform_isdb()
generate_summary()
generate_article_data()
