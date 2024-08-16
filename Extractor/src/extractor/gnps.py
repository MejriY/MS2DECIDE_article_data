import pandas as pd
import os
import requests
import json
from pathlib import Path
import xml.etree.ElementTree as ET
from collections import Counter

class GnpsAnnotationsFile:
    def __init__(self, json_file: str | os.PathLike):
        self.json_file = json_file

    def df(self):
        with open(self.json_file) as json_data:
            js = json.load(json_data)
            assert len(js) == 1
            (k, v), = js.items()
            df = pd.DataFrame(v)
            return df
    
    def summary(self):
        summary = self.df().loc[:, ["#Scan#", "Adduct", "ExactMass", "Precursor_MZ", "SpecMZ", "Charge", "SpecCharge", "MQScore", "INCHI", "INCHI_AUX", "InChIKey"]].rename(columns = {"#Scan#": "Id"})
        summary["Id"] = summary["Id"].astype(int)
        return summary.set_index("Id").sort_index()

class GnpsParametersFile:
    def __init__(self, xml_file: str | os.PathLike):
        self.xml_file = xml_file

    def params(self):
        with open(self.xml_file) as f:
            xml = f.read()
            tree = ET.parse(self.xml_file)
            root = tree.getroot()
            subtags = set([p.tag for p in root])
            assert subtags == {"parameter"}, subtags
            names_mult = Counter([p.attrib["name"] for p in root])
            dupls = set([names for names, count in names_mult.items() if count > 1])
            assert dupls == {} or dupls == {"upload_file_mapping"}, dupls

            d = {p.attrib["name"]: p.text for p in root if p.attrib["name"] not in dupls}
            return d
        
class GnpsFetcher:
    def fetch(task_id: str):
        url = f"https://gnps.ucsd.edu/ProteoSAFe/result_json.jsp?task={task_id}&view=view_all_annotations_DB"
        with requests.get(url) as r:
            r.raise_for_status()
            return r.content

    def fetch_and_save(task_id: str, path: str | os.PathLike):
        with open(path, "wb") as f:
            f.write(GnpsFetcher.fetch(task_id))

    def fetch_parameters(task_id: str):
        url = f"https://gnps.ucsd.edu/ProteoSAFe/ManageParameters?task={task_id}"
        with requests.get(url) as r:
            r.raise_for_status()
            return r.content

    def fetch_parameters_and_save(task_id: str, path: str | os.PathLike):
        with open(path, "wb") as f:
            f.write(GnpsFetcher.fetch_parameters(task_id))
    
class GnpsCacher:
    def __init__(self, cache_dir: str | os.PathLike):
        self.cache_dir = cache_dir
        
    def cache_retrieve(self, task_id: str):
        os.makedirs(self.cache_dir, exist_ok=True)
        path = self.cache_dir / f"{task_id}.json"
        if not path.exists():
            GnpsFetcher.fetch_and_save(task_id, path)
        return path
    
    def cache_retrieve_parameters(self, task_id: str):
        os.makedirs(self.cache_dir, exist_ok=True)
        path = self.cache_dir / f"{task_id}.xml"
        if not path.exists():
            GnpsFetcher.fetch_parameters_and_save(task_id, path)
        return path
    
class GnpsInchiScore:
    def __init__(self, all_annotations: Path, parameters: Path):
        self.summary = GnpsAnnotationsFile(all_annotations).summary()
        ps = GnpsParametersFile(parameters).params()
        self.min_peaks = int(ps["MIN_MATCHED_PEAKS_SEARCH"])
        self.max_delta_mass = float(ps["MAX_SHIFT_MASS"])
    
    @property
    def inchis(self):
        return self.summary.loc[:, "INCHI"]
    
    @property
    def scores(self):
        return self.summary.loc[:, "MQScore"]
    
    @property
    def inchis_scores_df(self):
        new_cols = {
            f"InChI GNPS; peaks ≥ {self.min_peaks}; Δ mass ≤ {self.max_delta_mass}": self.inchis,
            f"Score GNPS; peaks ≥ {self.min_peaks}; Δ mass ≤ {self.max_delta_mass}": self.scores,
        }
        return pd.DataFrame(new_cols)
