import os
import pandas as pd
import json
from extractor.gnps import GnpsCacher
from extractor.gnps import GnpsAnnotations
from extractor.gnps import GnpsInchiScore
from extractor.mgfs import MgfFiles
from rdkit import Chem
from rdkit.Chem import Descriptors
import matchms
from ms2decide.K import K
from ms2decide.Tanimotos import Tanimotos
from shutil import rmtree
from extractor.datadirs import *

def clean():
    assert rmtree.avoids_symlink_attacks
    rmtree(GENERATED_DIR)

def generate_gnps_input():
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    GENERATED_DIR_INPUT_GNPS.mkdir(parents=True, exist_ok=True)

    compounds_file = INPUT_DIR/"Compounds.tsv"
    compounds = pd.read_csv(compounds_file, sep="\t").set_index("Chemical name")
    names = set(compounds.index.to_list())
    assert len(names) == 96

    all_annotations = INPUT_DIR/"Mgf files/"
    mgfs = MgfFiles(all_annotations)
    assert mgfs.d.keys() == names, set(names) - set(mgfs.d.keys())

    inchis = compounds.loc[compounds["InChI"].notna(), "InChI"]
    compounds["Relative molecular weight"] = inchis.apply(lambda i: Descriptors.MolWt(Chem.inchi.MolFromInchi(i)))
    compounds["Precursor m/z"] = mgfs.precursors
    compounds["Retention time"] = mgfs.retentions
    compounds["Precursor m/z − relative molecular weight"] = compounds["Precursor m/z"] - compounds["Relative molecular weight"]

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
    # matchms.exporting.save_as_mgf(all_spectra, str(generated_dir / "All Matchms.mgf"))
    matchms.exporting.save_as_mgf(all_spectra, str(GENERATED_DIR_INPUT_GNPS / "All GNPS.mgf"), export_style="gnps")

    qt = pd.DataFrame()
    qt["row ID"] = compounds["Id"]
    qt["row m/z"] = compounds["Precursor m/z"]
    qt["row retention time"] = compounds["Retention time"] / 60.0
    qt["1.mzXML Peak area"] = 0
    qt.set_index("row ID", inplace=True)
    qt.to_csv(GENERATED_DIR_INPUT_GNPS / "Quantification table.csv")

def generate_summary():
    GENERATED_DIR_SUMMARY.mkdir(parents=True, exist_ok=True)
    
    task_ids_file = GENERATED_DIR_GNPS_TASKS / "Gnps task ids.json"
    with open(task_ids_file) as task_ids_data:
        task_ids = json.load(task_ids_data)

    compounds_file = INPUT_DIR/"Compounds.tsv"
    compounds_by_id = pd.read_csv(compounds_file, sep="\t").set_index("Id")
    ids = compounds_by_id.index
    assert len(set(ids.to_list())) == 96

    compounds_joined = compounds_by_id
    for task_id in task_ids:
        all_annotations_data = GnpsCacher(GENERATED_DIR_SUMMARY / "Fetched/").cache_retrieve_annotations_data(task_id)
        all_annotations = GnpsAnnotations(all_annotations_data)
        parameters = GnpsCacher(GENERATED_DIR_SUMMARY / "Fetched/").cache_retrieve_parameters(task_id)
        isc = GnpsInchiScore(all_annotations, parameters)
        assert isc.inchis.index.isin(ids).all()
        assert isc.scores.index.isin(ids.index).all()
        print(f"Joining task {task_id}, min peaks {isc.min_peaks}, max delta mass {isc.max_delta_mass}")
        task_df = isc.inchis_scores_df
        compounds_joined = compounds_joined.join(task_df)

    compounds_joined.to_csv(GENERATED_DIR_SUMMARY / "Compounds joined.csv")

def main():
    # todo
    pass

def compute():
    dict_inchi = {'GNPS': 'InChI=1S/C43H54N4O6/c1-7-23-17-25-20-43(22-53-52-6)39-28(15-16-47(40(23)43)41(25)48)27-13-14-34(50-4)36(38(27)45-39)31-18-29-24(8-2)21-46(3)33(35(29)42(49)51-5)19-30-26-11-9-10-12-32(26)44-37(30)31/h9-14,23-25,29,31,33,35,40-41,44-45,48H,7-8,15-21H2,1-6H3/t23-,24+,25-,29-,31-,33-,35?,40-,41?,43?/m0/s1',
                        'ISDB': 'InChI=1S/C23H27N2O3/c1-25-8-7-23-16-4-2-3-5-17(16)24-20(26)11-18-21(22(23)24)15(10-19(23)28-13-25)14(12-25)6-9-27-18/h2-6,15,18-19,21-22H,7-13H2,1H3/q+1',
                        'SIRIUS': 'InChI=1S/C43H52N4O5/c1-7-24-15-23-20-43(42(49)52-6)39-27(13-14-47(21-23)40(24)43)29-17-30(36(50-4)19-34(29)45-39)31-16-28-25(8-2)22-46(3)35(37(28)41(48)51-5)18-32-26-11-9-10-12-33(26)44-38(31)32/h8-12,17,19,23-24,28,31,35,37,40,44-45H,7,13-16,18,20-22H2,1-6H3'}
    
    tanimotos = Tanimotos(dict_inchi)
