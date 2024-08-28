import pandas as pd
import json
from extractor.compounds import Compounds
from extractor.gnps import IteratedQueries
from extractor.manufactured.mgfs import MgfFiles
from extractor.support import add_ranks_columns
import matchms
from ms2decide.MgfInstance import MgfInstance
from ms2decide.Tool import Tool
from ms2decide.IsdbAnnotation import get_cfm_annotation
from ms2decide.K import K
from ms2decide.Tanimotos import Tanimotos
from shutil import rmtree
from extractor.manufactured.datadirs import *
from collections import Counter


TASK_IDS_FILE = GENERATED_DIR_GNPS_TASKS / "Gnps task ids.json"


def get_task_ids():
    with open(TASK_IDS_FILE) as task_ids_data:
        task_ids = set(json.load(task_ids_data))
    return task_ids


def get_iterated_queries():
    task_ids = get_task_ids()
    qs = IteratedQueries.from_task_ids(task_ids, GENERATED_DIR_GNPS_TASKS_CACHED)
    return qs


def clean():
    assert rmtree.avoids_symlink_attacks
    rmtree(GENERATED_DIR)


def generate_gnps_input():
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    GENERATED_DIR_INPUTS.mkdir(parents=True, exist_ok=True)

    compounds = Compounds.from_tsv(INPUT_DIR / "Compounds.tsv")
    mgfs = MgfFiles(INPUT_DIR / "Mgf files/", compounds.df["Chemical name"])

    all_spectra = mgfs.all_spectra()
    # This exports PRECURSOR_MZ instead of PEPMASS, which apparently GNPS does not like.
    mgfs.export_all_spectra(GENERATED_DIR_INPUTS / "All Matchms.mgf")
    mgfs.export_all_spectra(GENERATED_DIR_INPUTS / "All GNPS.mgf", export_style="gnps")

    compounds.quantification_table_minutes(mgfs.precursors_series(), mgfs.retentions_seconds_series()).to_csv(GENERATED_DIR_INPUTS / "Quantification table.csv")


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

    compounds_file = INPUT_DIR / "Compounds.tsv"
    compounds_obj = Compounds.from_tsv(compounds_file)
    ids = compounds_obj.df.index
    assert len(ids.to_list()) == 96
    assert len(set(ids.to_list())) == 96
    names = compounds_obj.df["Chemical name"].to_list()
    assert len(names) == 96
    assert len(set(names)) == 96

    all_annotations = INPUT_DIR / "Mgf files/"
    mgfs = MgfFiles(all_annotations)
    assert mgfs.d.keys() == names, set(names) - set(mgfs.d.keys())

    compounds_obj.add_relative_molecular_weights()
    compounds_obj.add_precursors(mgfs.precursors)
    compounds_obj.add_retention_times(mgfs.retentions_sec)
    compounds_obj.add_diffs()

    compounds_joined = compounds_obj.df.join(get_iterated_queries().all_df())

    sirius_df = pd.read_csv(SIRIUS_DIR / "structure_identifications.tsv", sep="\t").set_index("mappingFeatureId")
    sirius_df["Score Sirius"] = sirius_df["ConfidenceScoreExact"].replace({float("-inf"): 0})
    sirius_df["Adduct Sirius"] = sirius_df["adduct"].map(lambda s: s.replace(" ", ""))
    sirius_df = sirius_df.rename(columns={"InChI": "InChI Sirius"}).loc[
        :, ["InChI Sirius", "Score Sirius", "Adduct Sirius"]
    ]
    compounds_joined = compounds_joined.join(sirius_df)

    compounds_joined["Adduct GNPS and Sirius"] = compounds_joined.apply(
        lambda r: str(dict(Counter(r[r.index.map(lambda c: c.startswith("Adduct "))]))), axis=1
    )

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
    to_emph = {k: "\\emph{" + str(k) + "}" for k in range(91, 97)}
    by_k = (
        compounds.sort_values("Rank min K")
        .loc[:, ["cg", "cs", "ci", "tgs", "tgi", "tsi", "K", "Ranks K"]]
        .rename(index=to_emph)
        .rename({"Ranks K": "Rank K"}, axis=1)
    )
    inchis = compounds.columns[compounds.columns.map(lambda s: "InChI" in s)]
    # Percents in smiles may cause problems with csvsimple
    compounds.drop(columns=inchis).rename(columns=lambda x: x.replace(" ", "")).replace(
        {",": ";", "N\\-demethyl": r"N\\Hyphdash{}demethyl"}, regex=True
    ).to_csv(GENERATED_DIR_ARTICLE / "Compounds.csv")
    compounds.rename(columns=lambda x: x.replace(" ", "").replace("-", "").replace("/", "")).replace(
        {"N-demethyl": r"N\\Hyphdash{}demethyl", "N-methyl": r"N\\Hyphdash{}methyl"}, regex=True
    ).to_csv(GENERATED_DIR_ARTICLE / "Compounds.tsv", sep="\t")
    by_k.rename(columns=lambda x: x.replace(" ", "")).to_csv(GENERATED_DIR_ARTICLE / "K.csv")
