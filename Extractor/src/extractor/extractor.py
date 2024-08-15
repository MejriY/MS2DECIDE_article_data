import os
import pandas as pd
import json
from extractor.gnps import GnpsAnnotationsFile
from extractor.gnps import GnpsCacher
from extractor.gnps import GnpsParametersFile
from extractor.gnps import GnpsInchiScore
from extractor.mgfs import MgfFiles
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem.Draw import rdDepictor
import matchms
from pathlib import Path
from ms2decide.K import K
from ms2decide.Tanimotos import Tanimotos
from shutil import rmtree

output_dir = Path("../Generated/") / "Manufactured case/"
assert rmtree.avoids_symlink_attacks
rmtree(output_dir, ignore_errors=True)

output_dir.mkdir(parents=True)

compounds_file = "../Manufactured case/Compounds.tsv"
compounds = pd.read_csv(compounds_file, sep="\t").set_index("Chemical name")
names = set(compounds.index.to_list())
assert len(names) == 96

all_annotations = Path("../Manufactured case/Mgf files/")
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

all_spectra = list()
for id in compounds["Id"]:
    name = compounds[compounds["Id"] == id].index[0]
    spectrum = d[name]
    all_spectra.append(spectrum)

matchms.exporting.save_as_mgf(all_spectra, str(output_dir / "All.mgf"))
# export_style="gnps" does not export the feature_ids. It does export PEPMASS, but apparently having the charge and the precursor_mz is ok to GNPS as well.
# matchms.exporting.save_as_mgf(all_spectra, str(output_dir / "All Gnps.mgf"), export_style="gnps")

qt = pd.DataFrame()
qt["row ID"] = compounds["Id"]
qt["row m/z"] = compounds["Precursor m/z"]
qt["row retention time"] = compounds["Retention time"] / 60.0
qt["1.mzXML Peak area"] = 0
qt.set_index("row ID", inplace=True)
# TODO find out how we found this value
qt.at[91, "row retention time"] = 721.696
qt.to_csv(output_dir / "Quantification table.csv")

dict_inchi = {'GNPS': 'InChI=1S/C43H54N4O6/c1-7-23-17-25-20-43(22-53-52-6)39-28(15-16-47(40(23)43)41(25)48)27-13-14-34(50-4)36(38(27)45-39)31-18-29-24(8-2)21-46(3)33(35(29)42(49)51-5)19-30-26-11-9-10-12-32(26)44-37(30)31/h9-14,23-25,29,31,33,35,40-41,44-45,48H,7-8,15-21H2,1-6H3/t23-,24+,25-,29-,31-,33-,35?,40-,41?,43?/m0/s1',
                      'ISDB': 'InChI=1S/C23H27N2O3/c1-25-8-7-23-16-4-2-3-5-17(16)24-20(26)11-18-21(22(23)24)15(10-19(23)28-13-25)14(12-25)6-9-27-18/h2-6,15,18-19,21-22H,7-13H2,1H3/q+1',
                      'SIRIUS': 'InChI=1S/C43H52N4O5/c1-7-24-15-23-20-43(42(49)52-6)39-27(13-14-47(21-23)40(24)43)29-17-30(36(50-4)19-34(29)45-39)31-16-28-25(8-2)22-46(3)35(37(28)41(48)51-5)18-32-26-11-9-10-12-33(26)44-38(31)32/h8-12,17,19,23-24,28,31,35,37,40,44-45H,7,13-16,18,20-22H2,1-6H3'}
 
tanimotos = Tanimotos(dict_inchi)

task_ids_file = "../Manufactured case/Gnps task ids.json"
with open(task_ids_file) as task_ids_data:
    task_ids = json.load(task_ids_data)

compounds_by_id = compounds.reset_index().set_index("Id")
compounds_joined = compounds_by_id

for task_id in task_ids[:20]:
    all_annotations = GnpsCacher(output_dir / "Fetched/").cache_retrieve(task_id)
    parameters = GnpsCacher(output_dir / "Fetched/").cache_retrieve_parameters(task_id)
    isc = GnpsInchiScore(all_annotations, parameters)
    assert isc.inchis.index.isin(compounds_by_id.index).all()
    assert isc.scores.index.isin(compounds_by_id.index).all()
    print(f"Joining task {task_id}, min peaks {isc.min_peaks}, max delta mass {isc.max_delta_mass}")
    task_df = isc.inchis_scores_df
    compounds_joined = compounds_joined.join(task_df)

compounds_joined.to_csv(output_dir / "Compounds joined.csv")


def to_image():
    hl_from_inchi = Chem.inchi.MolFromInchi(
        "InChI=1S/C15H24N2O2/c18-12-4-5-16-8-10-6-11(14(16)7-12)9-17-13(10)2-1-3-15(17)19/h10-14,18H,1-9H2/t10-,11-,12-,13+,14-/m0/s1"
    )

    # Draw.MolToImage(hl_from_inchi).save("hl_from_inchi.png")
    canon_smiles = pyAvalonTools.GetCanonSmiles(hl_from_inchi)
    hl_from_canon_smiles = Chem.MolFromSmiles(canon_smiles)
    # Draw.MolToImage(hl_from_canon_smiles).save("hl_from_canon_smiles.png")
    hl_from_smiles = Chem.MolFromSmiles("[C@@H]12CN3CC[C@H](O)C[C@H]3[C@H](CN3C(=O)CCC[C@H]23)C1")
    # Draw.MolToImage(hl_from_smiles).save("hl_from_smiles.png")

    d = rdMolDraw2D.MolDraw2DSVG(500, 500)
    s = rdMolDraw2D.PrepareAndDrawMolecule(d, hl_from_inchi)
    d.FinishDrawing()
    t = d.GetDrawingText()
    # with open("hl_from_inchi.svg", "w") as f:
    # f.write(t)

    t = Draw.MolToACS1996SVG(hl_from_inchi)
    # with open("hl_from_inchi_acs1996.svg", "w") as f:
    # f.write(t)

    # Draw.MolToSVG(hl_from_inchi)
    # ValueError: Bad Conformer Id
    # Attempt as per https://github.com/rdkit/rdkit/issues/4991
    # rdDepictor.Compute2DCoords(hl_from_inchi)
    # Draw.MolToSVG(hl_from_inchi)
    # ValueError: Bad Conformer Id
