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
from functools import singledispatchmethod
from io import BytesIO
from functools import cached_property
from functools import cache
from functools import total_ordering
from collections import OrderedDict
import numpy as np

def discounts():
    discs = _get_iterative_parameters()
    for d in discs.keys():
        discs[d][float("inf")] = discs[d][0]
    return discs


DISCOUNTS = discounts()


class GnpsAnnotations:
    _json_str: str
    _INDEX = "#Scan#"
    _COLS = [
        "Adduct",
        "ExactMass",
        "Precursor_MZ",
        "SpecMZ",
        "Charge",
        "SpecCharge",
        "MQScore",
        "INCHI",
        "INCHI_AUX",
        "InChIKey",
        "Smiles",
    ]

    def __init__(self, json_data: bytes | str | None = None):
        if isinstance(json_data, bytes):
            self._json_str = json_data.decode(errors="replace")
        else:
            self._json_str = json_data
        # assert self.inchis.index.isin(self.scores.index).all()
        # assert self.smiles.index.isin(self.scores.index).all()

    @classmethod
    def from_file(cls, json_file: bytes | str | os.PathLike):
        with open(json_file, "rb") as json_stream:
            return GnpsAnnotations(json_stream.read())

    @cache
    def raw_df(self):
        if self._json_str is None:
            v = []
        else:
            js = json.loads(self._json_str)
            assert len(js) == 1
            ((k, v),) = js.items()
        if len(v) == 0:
            df = pd.DataFrame(columns=[self._INDEX] + self._COLS).set_index(self._INDEX)
        else:
            df = pd.DataFrame(v).set_index(self._INDEX)
        int_index = df.index.astype(int)
        return df.set_index(int_index).sort_index()

    def _summary(self):
        summary = self.raw_df().loc[
            :,
            self._COLS,
        ]
        summary.index.rename("Id", inplace=True)
        return summary

    def ids(self):
        return self._summary().index

    def inchis(self):
        return self._summary()["INCHI"].rename("InChI")

    def smiles(self):
        return self._summary().loc[:, "Smiles"]

    def adducts(self):
        return self._summary().loc[:, "Adduct"]

    def scores(self):
        return self._summary().loc[:, "MQScore"].astype(float).rename("Score")
 
    @cache
    def inchis_smiles_series(self):
        return self._summary().apply(lambda x: GnpsInchiSmiles(x["INCHI"], x["Smiles"]), axis=1).rename("InchiSmiles")
    
    @cache
    def standard_inchis_series(self):
        return self.inchis_smiles_series().apply(lambda x: x.to_standard_inchi).rename("Standard InChI")

    def summary_df(self):
        return self._summary().join(self.standard_inchis_series()).rename(columns={"ExactMass": "Exact mass", "Precursor_MZ": "Precursor m/z", "MQScore": "Score", "INCHI": "InChI", "INCHI_AUX": "InChI aux"}).astype({"Exact mass": float, "Precursor m/z": float, "Score": float})
   
    def matches_series(self):
        return self.summary_df().join(self.inchis_smiles_series()).apply(lambda x: GnpsMatch(x["InchiSmiles"].to_mol, x["Standard InChI"], x["Score"]), axis=1)

class GnpsParameters:
    _xml_data: str | bytes

    def __init__(self, xml_data: bytes | str):
        self._xml_data = xml_data

    @classmethod
    def from_file(cls, xml_file: bytes | str | os.PathLike):
        with open(xml_file) as xml_stream:
            return GnpsParameters(xml_stream.read())

    def all(self) -> dict[str, str]:
        if hasattr(self, "_all"):
            return self._all

        if isinstance(self._xml_data, bytes):
            root = ET.parse(BytesIO(self._xml_data)).getroot()
        else:
            root = ET.fromstring(self._xml_data)
        subtags = set([p.tag for p in root])
        assert subtags == {"parameter"}, subtags
        names_mult = Counter([p.attrib["name"] for p in root])
        dupls = set([names for names, count in names_mult.items() if count > 1])
        assert dupls == {} or dupls == {"upload_file_mapping"}, dupls

        self._all = {p.attrib["name"]: p.text for p in root if p.attrib["name"] not in dupls}
        return self._all

    @property
    def min_peaks(self):
        return int(self.all()["MIN_MATCHED_PEAKS_SEARCH"])

    @property
    def max_delta_mass(self):
        raw = float(self.all()["MAX_SHIFT_MASS"])
        return raw if raw != 0 else float("inf")

    def to_query(self):
        return GnpsQuery(self.min_peaks, self.max_delta_mass)

class GnpsTaskFetcher:
    _task_id: str

    @classmethod
    def url_status(cls, task_id: str):
        return f"https://gnps.ucsd.edu/ProteoSAFe/status_json.jsp?task={task_id}"

    @classmethod
    def url_exact(cls, task_id: str):
        return f"https://gnps.ucsd.edu/ProteoSAFe/result_json.jsp?task={task_id}&view=view_all_annotations_DB"

    @classmethod
    def url_analog(cls, task_id: str):
        return f"https://gnps.ucsd.edu/ProteoSAFe/result_json.jsp?task={task_id}&view=view_all_analog_annotations_DB"

    @classmethod
    def url_parameters(cls, task_id: str):
        return f"https://gnps.ucsd.edu/ProteoSAFe/ManageParameters?task={task_id}"

    def __init__(self, task_id: str):
        self._task_id = task_id

    def done(self):
        if hasattr(self, "_done"):
            return True

        url = self.url_status(self._task_id)
        with requests.get(url) as r:
            r.raise_for_status()
            data = r.json()
            done: bool = data["status"] == "DONE"
            if done:
                self._done = True
            return done

    def fetch_exact(self):
        assert self.done()
        url = self.url_exact(self._task_id)
        with requests.get(url) as r:
            r.raise_for_status()
            answer = r.content
            assert isinstance(answer, bytes)
            return answer

    def fetch_analog(self):
        assert self.done()
        url = self.url_analog(self._task_id)
        with requests.get(url) as r:
            r.raise_for_status()
            answer = r.content
            assert isinstance(answer, bytes)
            return answer

    def fetch_parameters(self):
        url = self.url_parameters(self._task_id)
        with requests.get(url) as r:
            r.raise_for_status()
            answer = r.content
            assert isinstance(answer, bytes)
            return answer


class GnpsCachingTaskFetcher:
    def __init__(self, cache_dir: os.PathLike):
        self.cache_dir = cache_dir
        self.cache_dir_analog = cache_dir / "analog"
        self.cache_dir_exact = cache_dir / "exact"

    def exact_annotations(self, task_id: str):
        os.makedirs(self.cache_dir_exact, exist_ok=True)

        path = self.cache_dir_exact / f"{task_id}.json"
        if path.exists():
            with open(path, "rb") as f:
                json_data = f.read()
        else:
            fetcher = GnpsTaskFetcher(task_id)
            if not fetcher.done():
                return GnpsAnnotations()

            json_data = fetcher.fetch_exact()
            with open(path, "wb") as f:
                f.write(json_data)
        return GnpsAnnotations(json_data)

    def analog_annotations(self, task_id: str):
        os.makedirs(self.cache_dir_analog, exist_ok=True)

        path = self.cache_dir_analog / f"{task_id}.json"
        if path.exists():
            with open(path, "rb") as f:
                json_data = f.read()
        else:
            fetcher = GnpsTaskFetcher(task_id)
            if not fetcher.done():
                return GnpsAnnotations()

            json_data = fetcher.fetch_analog()
            with open(path, "wb") as f:
                f.write(json_data)
        return GnpsAnnotations(json_data)

    def parameters(self, task_id: str):
        os.makedirs(self.cache_dir, exist_ok=True)

        path = self.cache_dir / f"{task_id}.xml"
        if path.exists():
            with open(path, "rb") as f:
                data = f.read()
        else:
            fetcher = GnpsTaskFetcher(task_id)
            data = fetcher.fetch_parameters()
            with open(path, "wb") as f:
                f.write(data)
        return GnpsParameters(data)

    def queried(self, task_id: str):
        return GnpsQueried(self.parameters(task_id).to_query(), self.exact_annotations(task_id), self.analog_annotations(task_id))

@dataclass(frozen=True)
class GnpsInchiSmiles:
    raw_inchi: str
    smiles: str

    @cached_property
    def sanitized_inchi(self):
        removed = self.raw_inchi.replace('"', "").strip()
        if not removed.startswith("InChI="):
            return "InChI=" + removed
        return removed

    @cached_property
    def to_mol(self) -> Mol | None:
        i = Chem.inchi.MolFromInchi(self.sanitized_inchi)
        if i is None:
            mol = Chem.MolFromSmiles(self.smiles)
        else:
            mol = i
        return mol

    @cached_property
    def to_standard_inchi(self) -> str | None:
        mol = self.to_mol
        if mol is None:
            return None
        return Chem.inchi.MolToInchi(mol)


@dataclass(frozen=True)
@total_ordering
class GnpsQuery:
    min_peaks: int
    max_delta_mass: float

    def __lt__(self, other):
        return (-self.min_peaks, self.max_delta_mass) < (-other.min_peaks, other.max_delta_mass)

    @property
    def discount(self):
        return DISCOUNTS[self.min_peaks][self.max_delta_mass]


@dataclass(frozen=True)
class GnpsQueried:
    query: GnpsQuery
    exact_annotations: GnpsAnnotations
    analog_annotations: GnpsAnnotations

    @property
    def min_peaks(self):
        return self.query.min_peaks

    @property
    def max_delta_mass(self):
        return self.query.max_delta_mass

    @property
    def ids(self):
        exact_ids = self.exact_annotations.ids()
        analog_ids = self.analog_annotations.ids()
        return exact_ids.union(analog_ids).sort_values()

    def summary_df(self):
        cols = ["InChI", "Smiles", "Standard InChI", "Score", "Adduct"]
        query_descr = f"peaks ≥ {self.min_peaks}; Δ mass ≤ {self.max_delta_mass}"
        analogs = self.analog_annotations.summary_df().loc[:, cols].rename(lambda x: f"Analog {x} GNPS; {query_descr}", axis=1)
        exacts = self.exact_annotations.summary_df().loc[:, cols].rename(lambda x: f"Exact {x} GNPS; {query_descr}", axis=1)
        return analogs.join(exacts)


@dataclass(frozen=True)
class GnpsMatch:
    mol: Mol | None
    standard_inchi: str | None
    score: float

    def __post_init__(self):
        assert (self.mol is None) == (self.standard_inchi is None)
        assert 0 <= self.score <= 1


class IteratedQueries:
    def __init__(self, querieds: list[GnpsQueried]):
        self._all = OrderedDict(sorted({q.query: q for q in querieds}.items()))
        self._ids = set()
        for q in self._all.values():
            self._ids.update(q.ids.to_list())
    
    @classmethod
    def from_task_ids(cls, task_ids: list[str], cache_dir: os.PathLike):
        fetcher = GnpsCachingTaskFetcher(cache_dir)
        queries = [fetcher.queried(task_id) for task_id in task_ids]
        return cls(queries)
    
    def ids(self):
        return self._ids
    
    def _best_matches_discounted_series(self):
        dict = {id: self._best_match_discounted(id) for id in self.ids()}
        return pd.Series(dict, name="Best matches GNPS iterated").rename_axis("Id")
        # return pd.Index(self.ids()).map(lambda x: self._best_match_discounted(x)).rename("Best matches GNPS iterated")
    
    def _best_match_discounted(self, id):
        for q in self._all.values():
            match = q.analog_annotations.matches_series()[id] if id in q.analog_annotations.ids() else None
            if getattr(match, "mol", None) is not None:
                return GnpsMatch(match.mol, match.standard_inchi, match.score * q.query.discount)
            
    def all_df(self):
        best_matches_discounted = self._best_matches_discounted_series()
        stds = best_matches_discounted.map(lambda x: x.standard_inchi if x is not None else pd.NA).rename("Standard InChI GNPS iterated")
        scores = best_matches_discounted.map(lambda x: x.score if x is not None else np.NAN).rename("Score GNPS iterated discounted")
        assert scores.dtype == float
        all_series = [v.summary_df() for v in self._all.values()] + [stds, scores]
        return pd.concat(all_series, axis=1)
