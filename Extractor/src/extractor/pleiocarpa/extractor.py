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
    GENERATED_DIR_SUMMARY.mkdir(parents=True, exist_ok=True)

    task_ids_file = INPUT_DIR_GNPS_TASKS / "Gnps task ids.json"
    with open(task_ids_file) as task_ids_data:
        task_ids = set(json.load(task_ids_data))

    ts = GnpsTasks(GENERATED_DIR_SUMMARY / "Fetched/", task_ids)
    ts.load()
    compounds_joined = ts.all_matches()

    sirius_df = pd.read_csv(INPUT_DIR_SIRIUS / "structure_identifications.tsv", sep="\t").set_index("mappingFeatureId")
    sirius_df["Score Sirius"] = sirius_df["ConfidenceScoreExact"].replace({float("-inf"): 0})
    sirius_df = sirius_df.rename(columns={"InChI": "InChI Sirius"}).loc[:, ["InChI Sirius", "Score Sirius"]]
    compounds_joined = compounds_joined.join(sirius_df)

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
        inchi_g = compounds_joined.loc[id, "InChI GNPS iterated"]
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
            tgs = 0.25
            tsi = 0.25
        assert not pd.isna(inchi_i)
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

    compounds_joined.to_csv(GENERATED_DIR_SUMMARY / "Compounds joined.tsv", sep="\t")

    counts = compounds_joined.filter(like="Standard InChI GNPS; ").agg("count")
    counts.index.name = "Column"
    counts.name = "Count"
    counts.to_csv(GENERATED_DIR_SUMMARY / "Counts.tsv", sep="\t")

    by_k = compounds_joined.loc[:, ["cg", "cs", "ci", "tgs", "tgi", "tsi", "K", "Ranks K"]].sort_values("Ranks K")
    by_k.to_csv(GENERATED_DIR_SUMMARY / "Compounds by K.tsv", sep="\t")

# compute_isdb()
# transform_isdb()
generate_summary()
