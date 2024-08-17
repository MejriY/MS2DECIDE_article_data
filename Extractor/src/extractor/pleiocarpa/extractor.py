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

def generate_summary():
    GENERATED_DIR_SUMMARY.mkdir(parents=True, exist_ok=True)

    task_ids_file = INPUT_DIR_GNPS_TASKS / "Gnps task ids.json"
    with open(task_ids_file) as task_ids_data:
        task_ids = set(json.load(task_ids_data))

    ts = GnpsTasks(GENERATED_DIR_SUMMARY / "Fetched/", task_ids)
    ts.load()
    ts.all_matches().to_csv(GENERATED_DIR_SUMMARY / "All matches.tsv", sep="\t")


generate_summary()
