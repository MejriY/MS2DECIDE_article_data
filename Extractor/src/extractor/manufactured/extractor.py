import pandas as pd
import json
from extractor.compounds import Compounds
from extractor.gnps import IteratedQueries
from extractor.manufactured.mgfs import MgfFiles
import matchms
from ms2decide.MgfInstance import MgfInstance
from ms2decide.IsdbAnnotation import get_cfm_annotation
from shutil import rmtree
from extractor.manufactured.datadirs import *
from rdkit import RDLogger


TASK_IDS_FILE = GENERATED_DIR_GNPS_TASKS / "Gnps task ids.json"


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
    RDLogger.DisableLog("rdApp.*")

    GENERATED_DIR_TABLES.mkdir(parents=True, exist_ok=True)

    compounds = Compounds.from_tsv(INPUT_DIR / "Compounds.tsv")
    mgfs = MgfFiles(INPUT_DIR / "Mgf files/", compounds.df["Chemical name"])
    ids = compounds.df.index

    col1 = compounds.df.pop("Reported")
    compounds.df.insert(1, "Reported", col1)
    compounds.add_relative_molecular_weights()
    compounds.add_precursors(mgfs.precursors_series())
    compounds.add_retention_times(mgfs.retentions_seconds_series())
    compounds.add_diffs()

    compounds.join_iterated_queries(TASK_IDS_FILE, GENERATED_DIR_GNPS_TASKS_CACHED)
    compounds.join_sirius_data(SIRIUS_DIR / "structure_identifications.tsv")
    compounds.add_adduct_summary()
    compounds.join_isdb_data(GENERATED_DIR_ISDB / "ISDB-LOTUS annotations.tsv")
    compounds.join_k_tuples()
    compounds.add_ranks()
    compounds_joined = compounds.df


    compounds_joined.to_csv(GENERATED_DIR_TABLES / "Compounds joined.tsv", sep="\t")

    # y_df = support.y_manufactured_df().rename(columns=lambda x: x + " Yassine")
    # compounds_joined = compounds_joined.join(y_df)
    # compounds_joined["K diff"] = (compounds_joined["K"] - compounds_joined["K Yassine"]).round(4)
    # compounds_joined.to_csv(GENERATED_DIR_TABLES / "Compounds joined with Y.tsv", sep="\t")

    compounds.get_counts().to_csv(GENERATED_DIR_TABLES / "Counts.tsv", sep="\t")

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
