import os
import pandas as pd
import json
from extractor.gnps import GnpsAnnotationsFile
from extractor.gnps import GnpsCacher
from extractor.gnps import GnpsParametersFile
from extractor.gnps import GnpsInchiScore

compounds_file = "../Manufactured case/Compounds.tsv"
compounds = pd.read_csv(compounds_file, sep="\t").set_index("New id")
task_ids_file = "../Manufactured case/Gnps task ids.json"
with open(task_ids_file) as task_ids_data:
    task_ids = json.load(task_ids_data)

for task_id in task_ids: 
    p = GnpsCacher.cache_retrieve(task_id)
    x = GnpsCacher.cache_retrieve_parameters(task_id)
    isc = GnpsInchiScore(p, x)
    new_cols = {f"inchi_gnps_{isc.min_peaks}_{isc.max_delta_mass}": isc.inchis, f"score_gnps_{isc.min_peaks}_{isc.max_delta_mass}": isc.scores}
    assert isc.inchis.index == compounds.index
    assert isc.scores.index == compounds.index
    compounds.assign(**new_cols)

# def main():