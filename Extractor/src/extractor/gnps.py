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
        ((k, v),) = js.items()
        if len(v) == 0:
            return pd.DataFrame()
        df = pd.DataFrame(v)
        return df

    def summary(self):
        summary = (
            self.df()
            .loc[
                :,
                [
                    "#Scan#",
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
                ],
            ]
            .rename(columns={"#Scan#": "Id"})
        )
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
    def fetch_exact(task_id: str):
        url = f"https://gnps.ucsd.edu/ProteoSAFe/result_json.jsp?task={task_id}&view=view_all_annotations_DB"
        with requests.get(url) as r:
            r.raise_for_status()
            return r.content

    def fetch_analog(task_id: str):
        url = f"https://gnps.ucsd.edu/ProteoSAFe/result_json.jsp?task={task_id}&view=view_all_analog_annotations_DB"
        with requests.get(url) as r:
            r.raise_for_status()
            return r.content

    def fetch_exact_and_save(task_id: str, path: str | os.PathLike):
        json_data = GnpsFetcher.fetch_exact(task_id)
        if GnpsAnnotations(json_data).df().empty:
            print(f"Warning: task {task_id} has no annotations; not saving it.")
        else:
            with open(path, "wb") as f:
                f.write(json_data)
        return json_data

    def fetch_analog_and_save(task_id: str, path: str | os.PathLike):
        json_data = GnpsFetcher.fetch_analog(task_id)
        if GnpsAnnotations(json_data).df().empty:
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
        self.cache_dir_analog = cache_dir / "analog"
        self.cache_dir_exact = cache_dir / "exact"

    def cache_retrieve_analog_annotations(self, task_id: str):
        os.makedirs(self.cache_dir_analog, exist_ok=True)
        path = self.cache_dir_analog / f"{task_id}.json"
        if path.exists():
            with open(path) as f:
                json_data = f.read()
        else:
            json_data = GnpsFetcher.fetch_analog_and_save(task_id, path)
        return GnpsAnnotations(json_data)

    def cache_retrieve_exact_annotations(self, task_id: str):
        os.makedirs(self.cache_dir_exact, exist_ok=True)
        path = self.cache_dir_exact / f"{task_id}.json"
        if path.exists():
            with open(path) as f:
                json_data = f.read()
        else:
            json_data = GnpsFetcher.fetch_exact_and_save(task_id, path)
        return GnpsAnnotations(json_data)

    def cache_retrieve_parameters(self, task_id: str):
        os.makedirs(self.cache_dir, exist_ok=True)
        path = self.cache_dir / f"{task_id}.xml"
        if not path.exists():
            GnpsFetcher.fetch_parameters_and_save(task_id, path)
        return path


class GnpsInchiScore:
    def __init__(self, all_annotations: GnpsAnnotations, parameters: GnpsParametersFile):
        if all_annotations.df().empty:
            self.summary = pd.DataFrame()
        else:
            self.summary = all_annotations.summary()
        ps = parameters.params()
        self.min_peaks = int(ps["MIN_MATCHED_PEAKS_SEARCH"])
        self.max_delta_mass = float(ps["MAX_SHIFT_MASS"])
        assert self.inchis.index.isin(self.scores.index).all()
        assert self.smiles.index.isin(self.scores.index).all()

    @property
    def ids(self):
        return self.summary.index

    @property
    def inchis(self):
        if self.summary.empty:
            return pd.Series()
        return self.summary.loc[:, "INCHI"]

    @property
    def smiles(self):
        if self.summary.empty:
            return pd.Series()
        return self.summary.loc[:, "Smiles"]

    @property
    def inchis_smiles_df(self):
        new_cols = {"INCHI": self.inchis, "Smiles": self.smiles}
        return pd.DataFrame(new_cols)

    def standard_inchis(self):
        if not hasattr(self, "_standard_inchis"):
            self._standard_inchis = self.inchis_smiles_df.apply(
                lambda x: GnpsInchiSmiles(x["INCHI"], x["Smiles"]).to_standard_inchi(), axis=1
            )
        return self._standard_inchis

    @property
    def scores(self):
        if self.summary.empty:
            return pd.Series().astype(float)
        return self.summary.loc[:, "MQScore"].astype(float)

    @property
    def inchis_scores_df(self):
        new_cols = {
            f"InChI GNPS; peaks ≥ {self.attempt.min_peaks}; Δ mass ≤ {self.attempt.max_delta_mass}": self.inchis,
            f"Smiles GNPS; peaks ≥ {self.attempt.min_peaks}; Δ mass ≤ {self.attempt.max_delta_mass}": self.smiles,
            f"Standard InChI GNPS; peaks ≥ {self.attempt.min_peaks}; Δ mass ≤ {self.attempt.max_delta_mass}": self.standard_inchis(),
            f"Score GNPS; peaks ≥ {self.attempt.min_peaks}; Δ mass ≤ {self.attempt.max_delta_mass}": self.scores,
        }
        return pd.DataFrame(new_cols)

    @property
    def attempt(self):
        return GnpsIterativeAttempt(self.min_peaks, self.max_delta_mass if self.max_delta_mass != 0 else float("inf"))

    def short_matches(self):
        if not hasattr(self, "_short_matches_dict"):
            self._short_matches_dict = {id: GnpsReadableMatch(self.standard_inchis()[id], self.scores[id]) for id in self.ids}
        return self._short_matches_dict


# Still need to find out how to order descending by min peaks.
@dataclass(frozen=True, order=True)
class GnpsIterativeAttempt:
    min_peaks: int
    max_delta_mass: float

    @property
    def discount(self):
        return DISCOUNTS[self.min_peaks][self.max_delta_mass]


@dataclass(frozen=True)
class GnpsInchiSmiles:
    inchi: str
    smiles: str

    def sanitized_inchi(self):
        removed = self.inchi.replace('"', "").strip()
        if not removed.startswith("InChI="):
            return "InChI=" + removed
        return removed

    def to_mol(self):
        i = Chem.inchi.MolFromInchi(self.sanitized_inchi())
        if i is None:
            mol = Chem.MolFromSmiles(self.smiles)
        else:
            mol = i
        return mol

    def to_standard_inchi(self):
        mol = self.to_mol()
        if mol is None:
            return None
        return Chem.inchi.MolToInchi(mol)


@dataclass(frozen=True)
class GnpsReadableMatch:
    inchi: Mol | None
    score: float


class GnpsIteratedNp:
    def __init__(self, match_by_attempt: dict[GnpsIterativeAttempt, GnpsReadableMatch]):
        self.match_by_attempt = match_by_attempt

    def best_match_discounted(self):
        ordered_attempts = sorted(self.match_by_attempt.keys(), key=lambda x: (-x.min_peaks, x.max_delta_mass))
        for attempt in ordered_attempts:
            match = self.match_by_attempt[attempt]
            if match is not None and match.inchi is not None:
                # print(f"Best match for {attempt} with discount {attempt.discount}: {match.inchi} with score {match.score} of type {type(match.score)}")
                return GnpsReadableMatch(match.inchi, match.score * attempt.discount)


@dataclass
class GnpsTask:
    cache_dir: Path
    task_id: str

    def load(self):
        analog_annotations = GnpsCacher(self.cache_dir).cache_retrieve_analog_annotations(self.task_id)
        exact_annotations = GnpsCacher(self.cache_dir).cache_retrieve_exact_annotations(self.task_id)
        parameters_file = GnpsCacher(self.cache_dir).cache_retrieve_parameters(self.task_id)
        self.analog_isc = GnpsInchiScore(analog_annotations, GnpsParametersFile(parameters_file))
        self.exact_isc = GnpsInchiScore(exact_annotations, GnpsParametersFile(parameters_file))

    def inchis_scores_both_df(self):
        return pd.concat([self.analog_isc.inchis_scores_df, self.exact_isc.inchis_scores_df.rename(lambda x: x + " exact", axis=1)], axis=1)

    def inchis_scores_analog_df(self):
        return self.analog_isc.inchis_scores_df

    def inchis_scores_exact_df(self):
        return self.exact_isc.inchis_scores_df

    @property
    def attempt(self):
        assert self.analog_isc.attempt == self.exact_isc.attempt
        return self.analog_isc.attempt

    def __eq__(self, other):
        return self.task_id == other.task_id
    def __hash__(self):
        return hash(self.task_id)
        
@dataclass
class GnpsTasks:
    cache_dir: Path
    task_ids: set[str]

    def load(self):
        self.tasks = {task_id: GnpsTask(self.cache_dir, task_id) for task_id in self.task_ids}
        for task in self.tasks.values():
            task.load()

    def task_id_to_attempt_df(self):
        to_peaks = {task.task_id: task.attempt.min_peaks for task in self.tasks.values()}
        to_delta_mass = {task.task_id: task.attempt.max_delta_mass for task in self.tasks.values()}
        return pd.DataFrame([to_peaks, to_delta_mass]).T.rename(columns={0: "Min peaks", 1: "Max Δ mass"}).astype({"Min peaks": int}).sort_values(by=["Min peaks", "Max Δ mass"], ascending=[False, True]).rename_axis("Task id")
    
    def _ordered_attempts(self):
        if hasattr(self, "__ordered_attempts"):
            return self.__ordered_attempts
        all_attempts = [task.attempt for task in self.tasks.values()]
        self.__ordered_attempts = sorted(all_attempts, key=lambda x: (-x.min_peaks, x.max_delta_mass))
        return self.__ordered_attempts
    
    def inchis_scores_df(self):
        by_attempt = {task.attempt: task for task in self.tasks.values()}
        ordered_attempts = self._ordered_attempts()
        return pd.concat([by_attempt[attempt].inchis_scores_both_df() for attempt in ordered_attempts], axis=1)

    def _match_by_attempt_dict(self, id):
        match_by_attempt = {}
        for task in self.tasks.values():
            short_match_results = task.analog_isc.short_matches()
            if(id not in short_match_results.keys()):
                match = None
            else:
                match = short_match_results[id]
            if match is not None:
                match_by_attempt[task.analog_isc.attempt] = match
        return match_by_attempt

    def _best_match_discounted(self, id):
        match_by_attempt = self._match_by_attempt_dict(id)
        return GnpsIteratedNp(match_by_attempt).best_match_discounted()

    def ids(self):
        return self.inchis_scores_df().index
    
    def best_matches_discounted(self):
        return {id: self._best_match_discounted(id) for id in self.ids()}
    
    def _best_matches_inchi_series(self):
        return pd.Series(
        {i: v.inchi if v is not None else None for (i, v) in self.best_matches_discounted().items()}, name="Standard InChI GNPS iterated"
    )
    def _best_matches_scores_series(self):
        return pd.Series(
        {i: v.score if v is not None else None for (i, v) in self.best_matches_discounted().items()}, name="Score GNPS iterated discounted"
    )
    def best_matches_df(self):
        return pd.concat([self._best_matches_inchi_series(), self._best_matches_scores_series()], axis=1)
    
    def all_matches(self):
        df = pd.concat([self.inchis_scores_df(), self.best_matches_df()], axis=1).sort_index()
        df.index.name = "Id"
        return df