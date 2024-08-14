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

output_dir = Path("../Generated/") / "Manufactured case"
if output_dir.exists():
    for f in output_dir.iterdir():
        f.unlink()
    output_dir.rmdir()

output_dir.mkdir(parents=True)

compounds_file = "../Manufactured case/Compounds.tsv"
compounds = pd.read_csv(compounds_file, sep="\t").set_index("Chemical name")
names = set(compounds.index.to_list())
assert len(names) == 96

p = Path("../Manufactured case/Mgf files/")
mgfs = MgfFiles(p)
assert mgfs.d.keys() == names, set(names) - set(mgfs.d.keys())

inchis = compounds.loc[compounds["InChI"].notna(), "InChI"]
compounds["Relative molecular weight"] = inchis.apply(lambda i: Descriptors.MolWt(Chem.inchi.MolFromInchi(i)))
compounds["Precursor m/z"] = mgfs.precursors
compounds["Precursor m/z − relative molecular weight"] = compounds["Precursor m/z"] - compounds["Relative molecular weight"]

d = mgfs.d
for name in d.keys():
    spectrum = d[name]
    id = compounds.loc[name, "Id"]
    spectrum.set("feature_id", id)

all_spectra = list()
for id in compounds["Id"]:
    name = compounds[compounds["Id"] == id].index[0]
    spectrum = d[name]
    all_spectra.append(spectrum)

matchms.exporting.save_as_mgf(all_spectra, str(output_dir / "All Matchms.mgf"))
matchms.exporting.save_as_mgf(all_spectra, str(output_dir / "All Gnps.mgf"), export_style="gnps")

s = next(iter(d.values()))
s.get("precursor_mz")

# kf_from_inchi = Chem.inchi.MolFromInchi(
#     "InChI=1S/C15H10O6/c16-8-3-1-7(2-4-8)15-14(20)13(19)12-10(18)5-9(17)6-11(12)21-15/h1-6,16-18,20"
# )
hl_from_inchi = Chem.inchi.MolFromInchi(
    "InChI=1S/C15H24N2O2/c18-12-4-5-16-8-10-6-11(14(16)7-12)9-17-13(10)2-1-3-15(17)19/h10-14,18H,1-9H2/t10-,11-,12-,13+,14-/m0/s1"
)
Descriptors.MolWt(hl_from_inchi)

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


task_ids_file = "../Manufactured case/Gnps task ids.json"
with open(task_ids_file) as task_ids_data:
    task_ids = json.load(task_ids_data)

for task_id in task_ids:
    p = GnpsCacher.cache_retrieve(task_id)
    x = GnpsCacher.cache_retrieve_parameters(task_id)
    isc = GnpsInchiScore(p, x)
    new_cols = {
        f"inchi_gnps_{isc.min_peaks}_{isc.max_delta_mass}": isc.inchis,
        f"score_gnps_{isc.min_peaks}_{isc.max_delta_mass}": isc.scores,
    }
    assert isc.inchis.index == compounds.index
    assert isc.scores.index == compounds.index
    compounds.assign(**new_cols)

# def main():
