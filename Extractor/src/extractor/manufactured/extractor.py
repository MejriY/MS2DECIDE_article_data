import pandas as pd
import json
from extractor.compounds import Compounds
from extractor.manufactured.mgfs import MgfFiles
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


def generate_gnps_input():
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    GENERATED_DIR_INPUTS.mkdir(parents=True, exist_ok=True)

    compounds = Compounds.from_tsv(INPUT_DIR / "Compounds.tsv")
    mgfs = MgfFiles(INPUT_DIR / "Mgf files/", compounds.df["Chemical name"])

    # This exports PRECURSOR_MZ instead of PEPMASS, which apparently GNPS does not like.
    mgfs.export_all_spectra(GENERATED_DIR_INPUTS / "All Matchms.mgf")
    mgfs.export_all_spectra(GENERATED_DIR_INPUTS / "All GNPS.mgf", export_style="gnps")
    mgfs.export_all_spectra_cut(GENERATED_DIR_INPUTS / "All MZmine cut.mgf")

    compounds.quantification_table_minutes(mgfs.precursors_series(), mgfs.retentions_seconds_series()).to_csv(GENERATED_DIR_INPUTS / "Quantification table.csv")

def compute_isdb():
    support.compute_isdb(GENERATED_DIR_INPUTS / "All GNPS.mgf", GENERATED_DIR_ISDB / "ISDB-LOTUS annotations.tsv")

def generate_sirius_input():
    GENERATED_DIR_SIRIUS.mkdir(parents=True, exist_ok=True)

    mzmine_path = MZMINE_DIR / "All MZmine cut Sirius.mgf"
    mzmine_str = mzmine_path.read_text()
    mzmine_patched_str = mzmine_str.replace("CHARGE=1?", "CHARGE=1+")
    patched_path = GENERATED_DIR_SIRIUS / "MZmine Sirius patched.mgf"
    patched_path.write_text(mzmine_patched_str)
    matchmzmine = matchms.importing.load_from_mgf(str(patched_path), metadata_harmonization=False)
    input_spectra = list(matchmzmine)

    compounds = Compounds.from_tsv(INPUT_DIR / "Compounds.tsv")
    ids = compounds.df.index
    input_ids = [s.get("scans") for s in input_spectra]
    input_ids_floats = [float(i) for i in input_ids]
    assert set(input_ids_floats) == set(ids.tolist()), f"{set(input_ids_floats) - set(ids.tolist())} {set(ids.tolist()) - set(input_ids)}"
    assert len(input_spectra) == len(ids) * 2
    
    filtered_spectra = list()
    for spectrum in input_spectra:
        mslevel = float(spectrum.metadata_dict()["ms_level"])
        if(mslevel == 1):
            mzs = spectrum.mz
            pepmass = spectrum.metadata_dict()["pepmass"][0]
            kept = mzs >= pepmass
            metadata_copy = spectrum.metadata_dict().copy()
            metadata_copy["num_peaks"] = sum(kept)
            spectrum_copy = matchms.Spectrum(spectrum.mz[kept], spectrum.intensities[kept], metadata_copy, metadata_harmonization=False)
            assert spectrum_copy.mz.size <= spectrum.mz.size
        else:
            spectrum_copy = spectrum
        filtered_spectra.append(spectrum_copy)
    
    out = GENERATED_DIR_SIRIUS / "Sirius.mgf"
    MgfFiles.export_sirius(filtered_spectra, out)
    patched_path.unlink()

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
    compounds = Compounds.from_tsv(GENERATED_DIR_TABLES / "Compounds joined.tsv").df
    to_emph = {k: "\\emph{" + str(k) + "}" for k in range(91, 97)}
    new_names = compounds.apply(lambda x: x["Chemical name"] if(x["Reported"]) else "", axis=1).rename("Chemical name reported")
    by_k = (
        compounds.sort_values(by = ["Rank min K", "Id"])
        .loc[:, ["Chemical name", "cg", "cs", "ci", "tgs", "tgi", "tsi", "K", "Rank K"]].assign(**{"Chemical name reported": new_names})
        .rename(index=to_emph)
    )
    compounds.rename(columns=lambda x: x.replace(" ", "").replace("-", "").replace("/", "")).replace(
        {"N-demethyl": r"N\\Hyphdash{}demethyl", "N-methyl": r"N\\Hyphdash{}methyl"}, regex=True
    ).to_csv(GENERATED_DIR_ARTICLE / "Compounds.tsv", sep="\t")
    by_k.rename(columns=lambda x: x.replace(" ", "")).to_csv(GENERATED_DIR_ARTICLE / "K.tsv", sep="\t")
