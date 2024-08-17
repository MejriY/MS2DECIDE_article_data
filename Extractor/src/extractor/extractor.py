import os
import pandas as pd
import json
from extractor.gnps import GnpsCacher
from extractor.gnps import GnpsAnnotations
from extractor.gnps import GnpsIteratedNp
from extractor.gnps import GnpsParametersFile
from extractor.gnps import GnpsInchiScore
from extractor.mgfs import MgfFiles
from rdkit import Chem
from rdkit.Chem import Descriptors
import matchms
from ms2decide.MgfInstance import MgfInstance
from ms2decide.Tool import Tool
from ms2decide.IsdbAnnotation import get_cfm_annotation
from ms2decide.K import K
from ms2decide.Tanimotos import Tanimotos
from shutil import rmtree
from extractor.datadirs import *


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
    GENERATED_DIR_SUMMARY.mkdir(parents=True, exist_ok=True)

    task_ids_file = GENERATED_DIR_GNPS_TASKS / "Gnps task ids.json"
    with open(task_ids_file) as task_ids_data:
        task_ids = json.load(task_ids_data)

    compounds_file = INPUT_DIR / "Compounds.tsv"
    compounds_by_id = pd.read_csv(compounds_file, sep="\t").set_index("Id")
    ids = compounds_by_id.index
    assert len(set(ids.to_list())) == 96

    compounds_joined = compounds_by_id
    iterated_np_by_id_dict = {id: {} for id in ids}
    for task_id in task_ids:
        all_annotations = GnpsCacher(GENERATED_DIR_SUMMARY / "Fetched/").cache_retrieve_annotations(task_id)
        parameters_file = GnpsCacher(GENERATED_DIR_SUMMARY / "Fetched/").cache_retrieve_parameters(task_id)
        isc = GnpsInchiScore(all_annotations, GnpsParametersFile(parameters_file))
        assert isc.inchis.index.isin(ids).all()
        assert isc.scores.index.isin(ids).all()
        for id in ids:
            match = isc.match(id)
            iterated_np_by_id_dict[id][isc.attempt] = match.to_readable() if match is not None else None
        print(f"Joining task {task_id}, min peaks {isc.min_peaks}, max delta mass {isc.max_delta_mass}")
        task_df = isc.inchis_scores_df
        compounds_joined = compounds_joined.join(task_df)

    iterated_np_by_id = {id: GnpsIteratedNp(iterated_np_by_id_dict[id]) for id in ids}
    best_match_discounted_by_id = {
        i: v.best_match_discounted() if v is not None else None for (i, v) in iterated_np_by_id.items()
    }
    compounds_joined["InChI GNPS iterated"] = pd.Series(
        {
            i: (Chem.inchi.MolToInchi(v.inchi) if v is not None else None)
            for (i, v) in best_match_discounted_by_id.items()
        }
    )
    compounds_joined["Score GNPS iterated discounted"] = pd.Series(
        {i: v.score if v is not None else None for (i, v) in best_match_discounted_by_id.items()}
    )

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
        inchi_g = compounds_joined.loc[id, "InChI GNPS iterated"]
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

    add_ranks_columns(compounds_joined, "Rank min GNPS original", "Rank max GNPS original", "Ranks GNPS original", "Score GNPS; peaks ≥ 6; Δ mass ≤ 0.02")
    add_ranks_columns(compounds_joined, "Rank min GNPS iterated", "Rank max GNPS iterated", "Ranks GNPS iterated", "Score GNPS iterated discounted")
    add_ranks_columns(compounds_joined, "Rank min Sirius", "Rank max Sirius", "Ranks Sirius", "Score Sirius")
    add_ranks_columns(compounds_joined, "Rank min ISDB-LOTUS", "Rank max ISDB-LOTUS", "Ranks ISDB-LOTUS", "Score ISDB-LOTUS")
    add_ranks_columns(compounds_joined, "Rank min K", "Rank max K", "Ranks K", "K")

    compounds_joined.to_csv(GENERATED_DIR_SUMMARY / "Compounds joined.tsv", sep="\t")

    yassine_values = pd.read_csv(GENERATED_DIR / "output gnps with new file" / "Manufactured annotation.tsv", sep="\t").set_index("ID")
    yassine_values.columns = yassine_values.columns.to_series().apply(lambda x: x + " Yassine")
    compounds_and_yassine = compounds_joined.join(yassine_values)
    compounds_and_yassine["K - K Yassine"] = (compounds_and_yassine["K"] - compounds_and_yassine["K Yassine"]).apply(lambda x: round(x, 5))
    compounds_and_yassine["Rank max K - ranking by k Yassine"] = compounds_and_yassine["Rank max K"] - compounds_and_yassine["ranking by k Yassine"]

    compounds_and_yassine.to_csv(GENERATED_DIR_SUMMARY / "Compounds and Yassine.tsv", sep="\t")


def add_ranks_columns(compounds_joined, rank_min_column, rank_max_column, ranks_column, score_column):
    compounds_joined[rank_min_column] = (
        compounds_joined[score_column].astype(float).fillna(0).rank(method="min").astype(int)
    )
    compounds_joined[rank_max_column] = (
        compounds_joined[score_column].astype(float).fillna(0).rank(method="max").astype(int)
    )
    compounds_joined[ranks_column] = compounds_joined.apply(
        lambda x: (
            str(x[rank_min_column])
            if x[rank_min_column] == x[rank_max_column]
            else (str(x[rank_min_column]) + "–" + str(x[rank_max_column]))
        ),
        axis=1,
    )


def compute():
    dict_inchi = {
        "GNPS": "InChI=1S/C43H54N4O6/c1-7-23-17-25-20-43(22-53-52-6)39-28(15-16-47(40(23)43)41(25)48)27-13-14-34(50-4)36(38(27)45-39)31-18-29-24(8-2)21-46(3)33(35(29)42(49)51-5)19-30-26-11-9-10-12-32(26)44-37(30)31/h9-14,23-25,29,31,33,35,40-41,44-45,48H,7-8,15-21H2,1-6H3/t23-,24+,25-,29-,31-,33-,35?,40-,41?,43?/m0/s1",
        "ISDB": "InChI=1S/C23H27N2O3/c1-25-8-7-23-16-4-2-3-5-17(16)24-20(26)11-18-21(22(23)24)15(10-19(23)28-13-25)14(12-25)6-9-27-18/h2-6,15,18-19,21-22H,7-13H2,1H3/q+1",
        "SIRIUS": "InChI=1S/C43H52N4O5/c1-7-24-15-23-20-43(42(49)52-6)39-27(13-14-47(21-23)40(24)43)29-17-30(36(50-4)19-34(29)45-39)31-16-28-25(8-2)22-46(3)35(37(28)41(48)51-5)18-32-26-11-9-10-12-33(26)44-38(31)32/h8-12,17,19,23-24,28,31,35,37,40,44-45H,7,13-16,18,20-22H2,1-6H3",
    }

    tanimotos = Tanimotos(dict_inchi)
