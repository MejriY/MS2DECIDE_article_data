from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem.Draw import rdDepictor
from pathlib import PurePath
import json
import requests
from io import StringIO
import pandas as pd
from github import Github
from github import Auth
from github import ContentFile
import extractor.manufactured.datadirs as mdirs
from ms2decide.MgfInstance import MgfInstance
from ms2decide.IsdbAnnotation import get_cfm_annotation
import extractor.manufactured.extractor as manufextractor
import extractor.pleiocarpa.extractor as pleioextractor

def clean():
    manufextractor.clean()
    pleioextractor.clean()

def unreported_ones():
    compounds_manufactured = pd.read_csv(mdirs.GENERATED_DIR_TABLES / "Compounds joined.tsv", sep="\t").set_index("Id")
    return compounds_manufactured[compounds_manufactured["Chemical name"].map(lambda s: s.startswith("Unreported "))]

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

def creds(auth_path = PurePath(".")):
    file = PurePath(auth_path, "auth.json")
    with open(file) as f:
        auth = json.load(f)
        username = auth["username"]
        token = auth["token"]
    return (username, token)

def y_manufactured_df(auth_path = PurePath(".")):
    auth = Auth.Token(creds(PurePath("."))[1])
    with Github(auth=auth) as g:
        repo = g.get_repo("MejriY/MS2DECIDE_article_data")
        sha = "7ff8db2500f278770845be3fa70283242d1054ac"
        server_path = "/Generated/Manufactured case/output gnps with new file/Manufactured annotation.tsv"
        contents : ContentFile = repo.get_contents(server_path, ref=sha)
        assert contents.name == "Manufactured annotation.tsv", contents
        content_decoded = contents.decoded_content.decode("utf-8")
        assert content_decoded.startswith("ID\t"), content_decoded
    return pd.read_csv(StringIO(content_decoded), sep="\t").rename(columns={"ID": "Id"}).set_index("Id")

def compute_isdb(input_mgf, output_tsv):
    mgf = MgfInstance(input_mgf)
    l = get_cfm_annotation(mgf)
    inchis = pd.Series({k: m.inchi for (k, m) in l.items()})
    scores = pd.Series({k: m.score for (k, m) in l.items()})
    # Let’s stick to a slightly more reasonable number of significant digits, the file is not stable after export-import
    df = pd.DataFrame({"InChI": inchis, "Score": scores}).replace({"#": None}).astype({"Score": float}).assign(Score = lambda x: x["Score"].map(lambda f: round(f, 14), na_action="ignore"))
    df.index.name = "Id"
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_tsv, sep="\t")
