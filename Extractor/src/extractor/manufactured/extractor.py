import pandas as pd
import json
from extractor.compounds import Compounds
from extractor.manufactured.mgfs import MgfFiles
from extractor.manufactured.mgfs import name_id_df
from shutil import rmtree
from extractor.manufactured.datadirs import *
from rdkit import RDLogger
import extractor.support as support
import numpy
import matchms

def clean():
    assert rmtree.avoids_symlink_attacks
    if GENERATED_DIR.exists():
        rmtree(GENERATED_DIR)
    if GENERATED_DIR_ARTICLE.exists():
        rmtree(GENERATED_DIR_ARTICLE)


def generate_input():
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    GENERATED_DIR_INPUTS.mkdir(parents=True, exist_ok=True)

    mgfs_sirus = MgfFiles(INPUT_DIR / "Mgf_Sirius/")
    compounds = Compounds.from_tsv(INPUT_DIR / "Compounds.tsv", name_id_df(mgfs_sirus.names))
    assert set(compounds.df.index) == set(range(1, 96)), compounds.df.index
    names = compounds.df["Chemical name"]
    assert mgfs_sirus.names == set(names), set(names) - set(mgfs_sirus.names)

    mgfs_sirus.export_all_sirius(GENERATED_DIR_INPUTS / "All Sirius.mgf")

    mgfs_gnps = MgfFiles(INPUT_DIR / "FBMN/")
    compounds = Compounds.from_tsv(INPUT_DIR / "Compounds.tsv", name_id_df(mgfs_gnps.names))
    assert set(compounds.df.index) == set(range(1, 96)), compounds.df.index
    names = compounds.df["Chemical name"]
    assert mgfs_gnps.names == set(names), set(names) - set(mgfs_gnps.names)

    # This exports PRECURSOR_MZ instead of PEPMASS, which apparently GNPS does not like.
    # mgfs_gnps.export_all_level2(GENERATED_DIR_INPUTS / "All Matchms.mgf")
    mgfs_gnps.export_all_level2(GENERATED_DIR_INPUTS / "All GNPS.mgf", export_style="gnps")
    compounds.quantification_table_minutes(mgfs_gnps.precursors_series(), mgfs_gnps.retentions_seconds_series()).to_csv(GENERATED_DIR_INPUTS / "Quantification table.csv")

    # subsets = [list(range(83, 84)), list(range(84, 85))]
    # for subset in subsets:
    #     subset_str = f"{min(subset):02d}-{max(subset):02d}"
    #     mgfs.export_level2(subset, GENERATED_DIR_INPUTS / f"{subset_str} GNPS.mgf", export_style="gnps")
    #     compounds.quantification_table_minutes(mgfs.precursors_series(), mgfs.retentions_seconds_series()).loc[subset, :].to_csv(GENERATED_DIR_INPUTS / f"{subset_str} Quantification table.csv")

    compounds.df.to_csv(GENERATED_DIR_INPUTS / "Compounds.tsv", sep="\t")

def compute_isdb():
    support.compute_isdb(GENERATED_DIR_INPUTS / "All GNPS.mgf", GENERATED_DIR_ISDB / "ISDB-LOTUS annotations.tsv", tol=0.002)

def generate_summary():
    RDLogger.DisableLog("rdApp.*")

    GENERATED_DIR_TABLES.mkdir(parents=True, exist_ok=True)

    mgfs_sirius = MgfFiles(INPUT_DIR / "Mgf_Sirius/")
    mgfs_gnps = MgfFiles(INPUT_DIR / "FBMN/")
    compounds = Compounds.from_tsv(INPUT_DIR / "Compounds.tsv", name_id_df(mgfs_gnps.names))
    ids = compounds.df.index
    precursors = mgfs_gnps.precursors_series()
    retentions = mgfs_gnps.retentions_seconds_series()
    
    # col1 = compounds.df.pop("Reported")
    # compounds.df.insert(1, "Reported", col1)
    compounds.add_relative_molecular_weights()
    compounds.add_precursors(mgfs_gnps.precursors_series(), "GNPS")
    compounds.add_precursors(mgfs_gnps.precursors_series(), "Sirius")
    compounds.add_retention_times(retentions)
    # TODO
    # compounds.add_diffs()

    compounds.join_iterated_queries(INPUT_DIR_GNPS_TASKS / "Gnps task ids.json", GENERATED_DIR_GNPS_TASKS_CACHED)
    compounds.join_sirius_data(SIRIUS_DIR / "structure_identifications.tsv")
    compounds.add_adduct_summary()
    compounds.join_isdb_data(GENERATED_DIR_ISDB / "ISDB-LOTUS annotations.tsv")
    compounds.join_k_tuples()
    compounds.add_ranks()

    compounds.df.to_csv(GENERATED_DIR_TABLES / "Compounds joined.tsv", sep="\t")

    # y_df = support.y_manufactured_df().rename(columns=lambda x: x + " Yassine")
    # compounds_joined = compounds_joined.join(y_df)
    # compounds_joined["K diff"] = (compounds_joined["K"] - compounds_joined["K Yassine"]).round(4)
    # compounds_joined.to_csv(GENERATED_DIR_TABLES / "Compounds joined with Y.tsv", sep="\t")

    compounds.get_counts().to_csv(GENERATED_DIR_TABLES / "Counts.tsv", sep="\t")

    # by_k = compounds.df.sort_values("Rank min K").loc[:, ["cg", "cs", "ci", "tgs", "tgi", "tsi", "K", "Rank K"]]
    # by_k.to_csv(GENERATED_DIR_TABLES / "Compounds by K.tsv", sep="\t")


def generate_article_data():
    GENERATED_DIR_ARTICLE.mkdir(parents=True, exist_ok=True)
    compounds = Compounds.from_tsv(GENERATED_DIR_TABLES / "Compounds joined.tsv").df.replace(
        {"N-demethyl": r"N\\Hyphdash{}demethyl", "N-methyl": r"N\\Hyphdash{}methyl"}, regex=True
    )
    to_emph = {k: "\\emph{" + str(k) + "}" for k in range(91, 97)}
    new_names = compounds.apply(lambda x: x["Chemical name"] if(x["Reported"]) else "", axis=1).rename("Chemical name reported")
    by_k = (
        compounds.sort_values(by = ["Rank min K", "Id"])
        .loc[:, ["Chemical name", "cg", "cs", "ci", "tgs", "tgi", "tsi", "K", "Rank K"]].assign(**{"Chemical name reported": new_names})
        .rename(index=to_emph)
    )
    compounds.rename(columns=lambda x: x.replace(" ", "").replace("/", "").replace("-", "")).to_csv(GENERATED_DIR_ARTICLE / "Compounds.tsv", sep="\t")
    by_k.rename(columns=lambda x: x.replace(" ", "")).to_csv(GENERATED_DIR_ARTICLE / "K.tsv", sep="\t")
