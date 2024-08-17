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

    compounds_joined.to_csv(GENERATED_DIR_SUMMARY / "All matches.tsv", sep="\t")

# compute_isdb()
# transform_isdb()
generate_summary()
