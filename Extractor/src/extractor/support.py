from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem.Draw import rdDepictor

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
