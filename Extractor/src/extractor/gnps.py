import pandas as pd
import os
import requests
import json
from pathlib import Path
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass
from rdkit.Chem.rdchem import Mol
import rdkit.Chem as Chem
from ms2decide.ClosestGNPS import _get_iterative_parameters

DISCOUNTS = _get_iterative_parameters()

class GnpsAnnotations:
    def __init__(self, json_raw_data: str):
        self.json_raw_data = json_raw_data

    def from_file(json_file: str | os.PathLike):
        with open(json_file) as json_data:
            return GnpsAnnotations(json_data.read())
    
    def df(self):
        js = json.loads(self.json_raw_data)
        assert len(js) == 1
        (k, v), = js.items()
        if(len(v) == 0):
            return pd.DataFrame()
        df = pd.DataFrame(v)
        return df
    
    def summary(self):
        summary = self.df().loc[:, ["#Scan#", "Adduct", "ExactMass", "Precursor_MZ", "SpecMZ", "Charge", "SpecCharge", "MQScore", "INCHI", "INCHI_AUX", "InChIKey", "Smiles"]].rename(columns = {"#Scan#": "Id"})
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
        json_data = GnpsFetcher.fetch(task_id)
        if(GnpsAnnotations(json_data).df().empty):
            print(f"Warning: task {task_id} has no annotations; not saving it.")
        else:
            with open(path, "wb") as f:
                f.write(json_data)
        return json_data

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
        
    def cache_retrieve_annotations(self, task_id: str):
        os.makedirs(self.cache_dir, exist_ok=True)
        path = self.cache_dir / f"{task_id}.json"
        if path.exists():
            with open(path) as f:
                json_data = f.read()
        else:
            json_data = GnpsFetcher.fetch_and_save(task_id, path)
        return GnpsAnnotations(json_data)
    
    def cache_retrieve_parameters(self, task_id: str):
        os.makedirs(self.cache_dir, exist_ok=True)
        path = self.cache_dir / f"{task_id}.xml"
        if not path.exists():
            GnpsFetcher.fetch_parameters_and_save(task_id, path)
        return path
    
class GnpsInchiScore:
    def __init__(self, all_annotations: GnpsAnnotations, parameters: GnpsParametersFile):
        if(all_annotations.df().empty):
            self.summary = pd.DataFrame()
        else:
            self.summary = all_annotations.summary()
        ps = parameters.params()
        self.min_peaks = int(ps["MIN_MATCHED_PEAKS_SEARCH"])
        self.max_delta_mass = float(ps["MAX_SHIFT_MASS"])
    
    @property
    def ids(self):
        return self.summary.index
    
    @property
    def inchis(self):
        if(self.summary.empty):
            return pd.Series()
        return self.summary.loc[:, "INCHI"]
    
    @property
    def smiles(self):
        if(self.summary.empty):
            return pd.Series()
        return self.summary.loc[:, "Smiles"]
    
    @property
    def scores(self):
        if(self.summary.empty):
            return pd.Series()
        return self.summary.loc[:, "MQScore"]
    
    @property
    def inchis_scores_df(self):
        new_cols = {
            f"InChI GNPS; peaks ≥ {self.attempt.min_peaks}; Δ mass ≤ {self.attempt.max_delta_mass}": self.inchis,
            f"Smiles GNPS; peaks ≥ {self.attempt.min_peaks}; Δ mass ≤ {self.attempt.max_delta_mass}": self.smiles,
            f"Score GNPS; peaks ≥ {self.attempt.min_peaks}; Δ mass ≤ {self.attempt.max_delta_mass}": self.scores,
        }
        return pd.DataFrame(new_cols)

    @property
    def attempt(self):
        return GnpsIterativeAttempt(self.min_peaks, self.max_delta_mass if self.max_delta_mass != 0 else float('inf'))
    
    def match(self, id: int):
        if(id not in self.ids):
            return None
        return GnpsMatch(self.summary.loc[id, "INCHI"], self.summary.loc[id, "Smiles"], float(self.scores[id]))

@dataclass(frozen=True, order=True)
class GnpsIterativeAttempt:
    min_peaks: int
    max_delta_mass: float
    
    @property
    def discount(self):
        return DISCOUNTS[self.min_peaks][self.max_delta_mass]
    
@dataclass(frozen=True)
class GnpsMatch:
    inchi: str
    smiles: str
    score: float
    
    def sanitized_inchi(self):
        removed = self.inchi.replace('"', "").strip()
        if(not removed.startswith("InChI=")):
            return "InChI=" + removed
        return removed
    
    def to_readable(self):
        i = Chem.inchi.MolFromInchi(self.sanitized_inchi())
        if(i is None):
            mol = Chem.MolFromSmiles(self.smiles)
        else:
            mol = i
        return GnpsReadableMatch(mol, self.score)
    
@dataclass(frozen=True)
class GnpsReadableMatch:
    inchi: Mol | None
    score: float
    
class GnpsIteratedNp:
    def __init__(self, match_by_attempt: dict[GnpsIterativeAttempt, GnpsReadableMatch]):
        self.match_by_attempt = match_by_attempt
    
    def best_match_discounted(self):
        for attempt in self.match_by_attempt.keys():
            match = self.match_by_attempt[attempt]
            if(match is not None and match.inchi is not None):
                # print(f"Best match for {attempt} with discount {attempt.discount}: {match.inchi} with score {match.score} of type {type(match.score)}")
                return GnpsReadableMatch(match.inchi, match.score * attempt.discount)