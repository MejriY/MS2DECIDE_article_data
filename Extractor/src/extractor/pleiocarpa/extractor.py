import os
import pandas as pd
import json
from extractor.gnps import GnpsAnnotations
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
from extractor.compounds import Compounds
from rdkit import RDLogger
from extractor.mgf import Mgf

def clean():
    assert rmtree.avoids_symlink_attacks
    if GENERATED_DIR.exists():
        rmtree(GENERATED_DIR)
    if GENERATED_DIR_ARTICLE.exists():
        rmtree(GENERATED_DIR_ARTICLE)

def compute_isdb():
    support.compute_isdb(INPUT_DIR / "2 - MZmine" / "Pleiocarpa.mgf", GENERATED_DIR_ISDB / "ISDB-LOTUS annotations.tsv")

def generate_summary():
    RDLogger.DisableLog("rdApp.*")

    GENERATED_DIR_TABLES.mkdir(parents=True, exist_ok=True)

    mgf_file_1 = INPUT_DIR / "2 - MZmine" / "Pleiocarpa.mgf"
    df_1 = Mgf(mgf_file_1).df.rename(lambda x: x + " raw", axis=1)

    mgf_file_2 = INPUT_DIR / "2 - MZmine" / "Pleiocarpa_sirius.mgf"
    df_2 = Mgf(mgf_file_2).df.rename(lambda x: x + " Sirius", axis=1)

    assert len(df_1) == len(df_2)
    compounds_joined = pd.concat([df_1, df_2], axis=1)
    compounds_joined.index = compounds_joined.index.astype(int)
    assert (abs(compounds_joined["retention_time raw"] - compounds_joined["retention_time Sirius"]) <= 0.005).all()
    compounds_joined = compounds_joined.drop(columns=["retention_time Sirius"]).rename(columns={"retention_time raw": "Retention time (seconds)"})
    assert (abs(compounds_joined["precursor_mz raw"] - compounds_joined["precursor_mz Sirius"]) <= 0.0005).all()
    compounds_joined = compounds_joined.drop(columns=["precursor_mz Sirius"]).rename(columns={"precursor_mz raw": "Precursor m/z"})
    assert len(compounds_joined.columns) == 6
    compounds_joined = compounds_joined.reindex(columns = ["Precursor m/z", "Retention time (seconds)", "Charge raw", "Charge Sirius", "MS level raw", "MS level Sirius"])
    compounds_joined.index.name = "Id"

    compounds = Compounds(compounds_joined)

    compounds.join_iterated_queries(INPUT_DIR_GNPS_TASKS / "Gnps task ids.json", GENERATED_DIR_TABLES / "Fetched/")
    compounds.join_sirius_data(INPUT_DIR_SIRIUS / "structure_identifications.tsv")
    compounds.add_adduct_summary()
    compounds.join_isdb_data(GENERATED_DIR_ISDB / "ISDB-LOTUS annotations.tsv")
    compounds.join_k_tuples()
    compounds.add_ranks()
    compounds.prefix_semantic_ids()
    compounds_joined = compounds.df.reset_index().set_index("Semantic id")

    descriptors_df = pd.read_csv(INPUT_DIR / "Descriptors.tsv", sep="\t").set_index("Semantic id").rename({"Adduct": "Adduct manual"}, axis=1)
    compounds_joined = compounds_joined.join(descriptors_df)
    
    compounds_joined.to_csv(GENERATED_DIR_TABLES / "Compounds joined.tsv", sep="\t")

    counts = compounds_joined.filter(like="Standard InChI GNPS").agg("count")
    counts.index.name = "Column"
    counts.name = "Count"
    counts.to_csv(GENERATED_DIR_TABLES / "Counts.tsv", sep="\t")

    # y_df = pd.read_csv(INPUT_DIR / "Pleiocarpa_annotation.tsv", sep="\t").rename(columns = {"ID": "Id"}).set_index("Id").rename(columns= lambda x: x + " Yassine")
    # compounds_joined = compounds_joined.join(y_df)
    # compounds_joined["K diff"] = (compounds_joined["K"] - compounds_joined["K Yassine"]).round(4)
    # compounds_joined.to_csv(GENERATED_DIR_TABLES / "Compounds with Yassine.tsv", sep="\t")

def generate_network_data():
# merely importing this creates a logs directory
# from extractor.network import Network
    # Unfinished business
    compounds = Compounds.from_tsv(GENERATED_DIR_TABLES / "Compounds joined.tsv").df
#     net = Network(INPUT_DIR_CYTOSCAPE / "network_k.cys", compounds)
#     net.export()

def generate_article_data():
    GENERATED_DIR_ARTICLE.mkdir(parents=True, exist_ok=True)

    # unreported_df = support.unreported_ones()
    # rounded_unreported_masses = unreported_df["Precursor m/z"].map(lambda p: round(p))

    compounds = Compounds.from_tsv(GENERATED_DIR_TABLES / "Compounds joined.tsv").df
    # compounds["Unreported"] = compounds["Precursor m/z"].map(lambda p: round(p) in rounded_unreported_masses)
    adduct_series = compounds.apply(lambda row: (row["Adduct manual"] if pd.notna(row["Adduct manual"]) else row["Adduct GNPS and Sirius"]), axis=1)
    by_k = (
        compounds.sort_values(by = ["Rank min K", "Precursor m/z", "Retention time (seconds)"]))
    by_k.insert(0, "Adduct", adduct_series)
    by_k = by_k.loc[:, ["Semantic id", "Adduct", "cg", "cs", "ci", "tgs", "tgi", "tsi", "K", "Ranks K", "Type of node", "Comments"]]
    by_k.rename(columns=lambda x: x.replace(" ", "")).replace({"{": "", "}": "", "'": ""}, regex=True).to_csv(GENERATED_DIR_ARTICLE / "K.tsv", sep="\t")
    
    compounds.rename(columns=lambda x: x.replace(" ", "")).to_csv(GENERATED_DIR_ARTICLE / "Compounds.tsv", sep="\t")
    