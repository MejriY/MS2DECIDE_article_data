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
from extractor import support
from rdkit import Chem
from rdkit.Chem import Descriptors
import matchms
from ms2decide.MgfInstance import MgfInstance
from ms2decide.Tool import Tool
from ms2decide.IsdbAnnotation import get_cfm_annotation
from ms2decide.K import K
from ms2decide.Tanimotos import Tanimotos
from shutil import rmtree
from extractor.manufactured.datadirs import *


def clean():
    assert rmtree.avoids_symlink_attacks
    rmtree(GENERATED_DIR)


def generate_gnps_input():
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    GENERATED_DIR_INPUTS.mkdir(parents=True, exist_ok=True)

    compounds_file = INPUT_DIR / "Compounds.tsv"
    compounds = pd.read_csv(compounds_file, sep="\t").set_index("Chemical name")
    names = set(compounds.index.to_list())
    assert len(names) == 96

    all_annotations = INPUT_DIR / "Mgf files/"
    mgfs = MgfFiles(all_annotations)
    assert mgfs.d.keys() == names, set(names) - set(mgfs.d.keys())

    inchis = compounds.loc[compounds["InChI"].notna(), "InChI"]
    compounds["Relative molecular weight"] = inchis.apply(lambda i: Descriptors.MolWt(Chem.inchi.MolFromInchi(i)))
    compounds["Precursor m/z"] = mgfs.precursors
    compounds["Retention time"] = mgfs.retentions
    compounds["Precursor m/z − relative molecular weight"] = (
        compounds["Precursor m/z"] - compounds["Relative molecular weight"]
    )

    d = mgfs.d
    for name in d.keys():
        spectrum = d[name]
        id = compounds.loc[name, "Id"]
        spectrum.set("scans", id)
        # No apparent effect when exported; seems that we need to build the spectrum using Spectrum(mz=sp.peaks.mz,intensities=sp.peaks.intensities,metadata=m).
        # spectrum.set("MSLEVEL", 2)

    all_spectra = list()
    for id in compounds["Id"]:
        name = compounds[compounds["Id"] == id].index[0]
        spectrum = d[name]
        all_spectra.append(spectrum)

    # This export PRECURSOR_MZ instead of PEPMASS, which apparently GNPS does not like.
    matchms.exporting.save_as_mgf(all_spectra, str(GENERATED_DIR_INPUTS / "All Matchms.mgf"))
    matchms.exporting.save_as_mgf(all_spectra, str(GENERATED_DIR_INPUTS / "All GNPS.mgf"), export_style="gnps")

    qt = pd.DataFrame()
    qt["row ID"] = compounds["Id"]
    qt["row m/z"] = compounds["Precursor m/z"]
    qt["row retention time"] = compounds["Retention time"] / 60.0
    qt["1.mzXML Peak area"] = 0
    qt.set_index("row ID", inplace=True)
    qt.to_csv(GENERATED_DIR_INPUTS / "Quantification table.csv")


def compute_isdb():
    GENERATED_DIR_ISDB.mkdir(parents=True, exist_ok=True)

    mgf = MgfInstance(GENERATED_DIR_INPUTS / "All GNPS.mgf")
    l = get_cfm_annotation(mgf)
    inchis = pd.Series({k: m.inchi for (k, m) in l.items()})
    scores = pd.Series({k: m.score for (k, m) in l.items()})
    df = pd.DataFrame({"InChI": inchis, "Score": scores})
    df.index.name = "Id"
    df.to_csv(GENERATED_DIR_ISDB / "ISDB-LOTUS annotations.tsv", sep="\t")


def generate_summary():
    GENERATED_DIR_TABLES.mkdir(parents=True, exist_ok=True)

    task_ids_file = GENERATED_DIR_GNPS_TASKS / "Gnps task ids.json"
    with open(task_ids_file) as task_ids_data:
        task_ids = set(json.load(task_ids_data))

    compounds_file = INPUT_DIR / "Compounds.tsv"
    compounds_by_id = pd.read_csv(compounds_file, sep="\t").set_index("Id")
    ids = compounds_by_id.index
    assert len(set(ids.to_list())) == 96

    ts = GnpsTasks(GENERATED_DIR_TABLES / "Fetched/", task_ids)
    ts.load()
    compounds_joined = compounds_by_id.join(ts.all_matches())

    sirius_df = pd.read_csv(SIRIUS_DIR / "structure_identifications.tsv", sep="\t").set_index("mappingFeatureId")
    sirius_df["Score Sirius"] = sirius_df["ConfidenceScoreExact"].replace({float("-inf"): 0})
    sirius_df = sirius_df.rename(columns={"InChI": "InChI Sirius"}).loc[:, ["InChI Sirius", "Score Sirius"]]
    compounds_joined = compounds_joined.join(sirius_df)

    isdb_df = (
        pd.read_csv(GENERATED_DIR_ISDB / "ISDB-LOTUS annotations.tsv", sep="\t")
        .set_index("Id")
        .rename(columns={"InChI": "InChI ISDB-LOTUS", "Score": "Score ISDB-LOTUS"})
    )
    compounds_joined = compounds_joined.join(isdb_df)

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

    compounds_joined.to_csv(GENERATED_DIR_TABLES / "Compounds joined.tsv", sep="\t")

    # y_df = support.y_manufactured_df().rename(columns=lambda x: x + " Yassine")
    # compounds_joined = compounds_joined.join(y_df)
    # compounds_joined["K diff"] = (compounds_joined["K"] - compounds_joined["K Yassine"]).round(4)
    # compounds_joined.to_csv(GENERATED_DIR_TABLES / "Compounds joined with Y.tsv", sep="\t")
    
    # selected_columns = compounds_joined.columns.map(lambda s: (s.startswith("Standard InChI GNPS; ")) or (s == "Id"))
    standards = compounds_joined.filter(like="Standard InChI GNPS")
    # columns = standards.columns
    # es = columns.map(lambda s: s.endswith(" exact"))
    # ces = columns[es]
    # nes = columns.map(lambda s: not s.endswith(" exact"))
    # cnes = columns[nes]
    # list(ces) + list(cnes)
    counts = standards.agg("count")
    counts.index.name = "Column"
    counts.name = "Count"
    counts.to_csv(GENERATED_DIR_TABLES / "Counts.tsv", sep="\t")

    by_k = compounds_joined.sort_values("Rank min K").loc[:, ["cg", "cs", "ci", "tgs", "tgi", "tsi", "K", "Ranks K"]]
    by_k.to_csv(GENERATED_DIR_TABLES / "Compounds by K.tsv", sep="\t")



def generate_article_data():
    GENERATED_DIR_ARTICLE.mkdir(parents=True, exist_ok=True)
    compounds = pd.read_csv(GENERATED_DIR_TABLES / "Compounds joined.tsv", sep="\t").set_index("Id")
    assert (compounds["Rank min K"] == compounds["Rank max K"]).all()
    by_k = (
        compounds
        .sort_values("Rank min K").loc[:, ["cg", "cs", "ci", "tgs", "tgi", "tsi", "K", "Ranks K"]]
        .rename({"Ranks K": "Rank K"}, axis=1)
    )
    compounds.rename(columns=lambda x: x.replace(" ", "_")).to_csv(GENERATED_DIR_ARTICLE / "Compounds.csv")
    by_k.rename(columns=lambda x: x.replace(" ", "_")).to_csv(GENERATED_DIR_ARTICLE / "K.csv")
