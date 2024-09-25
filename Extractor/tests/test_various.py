import pandas as pd
from extractor import support
from extractor.mgf import Mgf
from pathlib import PurePath
import json
import xml.etree.ElementTree as ET
import requests
from github import Github
from github import Auth
import urllib
import base64
from github import ContentFile
import pytest
from json import JSONDecodeError
from pathlib import Path
from rdkit import Chem
from rdkit import RDLogger
import re
import matchms
import tempfile

def test_read_df_manual():
    json_file = "tests/resources/26a5cbca3e844cc0b126f992c69df832.json"
    with open(json_file) as json_data:
        js = json.load(json_data)
        assert len(js) == 1
        (k, v), = js.items()
        df = pd.DataFrame(v)
    assert df.shape == (71, 48)

def test_parse_parameters_manual():
    xml_file = "tests/resources/26a5cbca3e844cc0b126f992c69df832.xml"
    with open(xml_file) as f:
        xml = f.read()
        tree = ET.parse(xml_file)
        root = tree.getroot()
        p0 = root[0]
        assert p0.tag == "parameter"
        assert p0.attrib["name"] == "ANALOG_SEARCH"

def test_rdkit_spits_out(capfd):
    RDLogger.DisableLog('rdApp.info')
    assert Chem.inchi.MolFromInchi("pure garbage") is None
    captured = capfd.readouterr()
    assert captured.out == ""
    assert re.compile(r"\[..:..:..\] ERROR: \n\n").match(captured.err)

def test_rdkit_silent(capfd):
    # https://github.com/rdkit/rdkit/issues/2320
    RDLogger.DisableLog("rdApp.*")
    assert Chem.inchi.MolFromInchi("pure garbage") is None
    captured = capfd.readouterr()
    assert captured.out == ""
    assert captured.err == ""

def test_load_invalid_utf8():
    input_path = Path("tests/resources/Invalid UTF-8 61f141ad8c3449fdb78a5086630b40f5.json")
    with open(input_path, encoding="utf-8") as f:
        with pytest.raises(UnicodeDecodeError):
            json_data = json.load(f)
    with open(input_path, encoding="utf-8") as f:
        with pytest.raises(UnicodeDecodeError):
            json_data = json.loads(f.read())
    with open(input_path, "rb") as f:
        bytes = f.read()
    with pytest.raises(UnicodeDecodeError):
        json.loads(bytes)
    with pytest.raises(UnicodeDecodeError):
        json.loads(bytes.decode("utf-8"))
    with pytest.raises(JSONDecodeError):
        json.loads(bytes.decode("utf-8", errors="backslashreplace"))
    json_data = json.loads(bytes.decode("utf-8", errors="replace"))
    assert len(json_data) == 1
    json_objects = json_data["blockData"]
    assert len(json_objects) >= 100, len(json_objects)
    for json_object in json_objects:
        if(json_object["#Scan#"]=="1158"):
            pbl_object = json_object
            break
    assert pbl_object["#Scan#"] == "1158", pbl_object
    assert pbl_object["Compound_Name"] == "(1S,5R,9S,13R)-14-formyl-5,9-dimethyltetracyclo[11.2.1.0�,�???.0???,???]hexadec-14-ene-5-carboxylic acid", pbl_object["Compound_Name"]

def test_direct_dl_fails():
    credsd = support.creds(PurePath("."))
    sha = "7ff8db2500f278770845be3fa70283242d1054ac"
    url = f"https://github.com/MejriY/MS2DECIDE_article_data/raw/{sha}/Generated/Manufactured%20case/output%20gnps%20with%20new%20file/Manufactured%20annotation.tsv"
    answer = requests.get(url, auth=credsd)
    assert answer.status_code == 404, answer.status_code

def test_list_root():
    auth = Auth.Token(support.creds(PurePath("."))[1])
    with Github(auth=auth) as g:
        repo = g.get_repo("MejriY/MS2DECIDE_article_data")
        sha = "7ff8db2500f278770845be3fa70283242d1054ac"
        server_path = "/"
        contents = repo.get_contents(urllib.parse.quote(server_path), ref=sha)
        assert len(contents) == 6, contents

def test_list_dir():
    auth = Auth.Token(support.creds(PurePath("."))[1])
    with Github(auth=auth) as g:
        repo = g.get_repo("MejriY/MS2DECIDE_article_data")
        sha = "7ff8db2500f278770845be3fa70283242d1054ac"
        server_path = "/Generated"
        contents = repo.get_contents(urllib.parse.quote(server_path), ref=sha)
        assert len(contents) == 1, contents

def test_list_subdir():
    auth = Auth.Token(support.creds(PurePath("."))[1])
    with Github(auth=auth) as g:
        repo = g.get_repo("MejriY/MS2DECIDE_article_data")
        sha = "7ff8db2500f278770845be3fa70283242d1054ac"
        server_path = "/Generated/Manufactured case/output gnps with new file"
        contents = repo.get_contents(server_path, ref=sha)
        assert len(contents) == 2, contents

def test_dl():
    auth = Auth.Token(support.creds(PurePath("."))[1])
    with Github(auth=auth) as g:
        repo = g.get_repo("MejriY/MS2DECIDE_article_data")
        sha = "7ff8db2500f278770845be3fa70283242d1054ac"
        server_path = "/Generated/Manufactured case/output gnps with new file/Manufactured annotation.tsv"
        contents : ContentFile = repo.get_contents(server_path, ref=sha)
        assert contents.name == "Manufactured annotation.tsv", contents
        content_encoded = contents.content
        content_decoded = contents.decoded_content.decode("utf-8")
        assert content_decoded.startswith("ID\t"), content_decoded

def test_pepmass_digits():
    input_path = Path("tests/resources/Pleiocarpa.mgf")
    spectras = matchms.importing.load_from_mgf(str(input_path))
    (s, ) = [s for s in spectras if s.metadata_dict().get("feature_id") == "39"]
    assert s.metadata_dict()["precursor_mz"] == 141.104
    mgf = Mgf(input_path)
    print(mgf.df.columns)
    assert mgf.df.loc[39, "precursor_mz"] == 141.104

def test_export_matchms():
    input_path = Path("tests/resources/Alchorneine.mgf")
    output_path = Path("out.mgf")
    spectras = matchms.importing.load_from_mgf(str(input_path))
    (s, ) = spectras
    assert s.metadata_dict()["precursor_mz"] == 222.1597604
    output_path.unlink(missing_ok=True)
    matchms.exporting.save_spectra([s], str(output_path), export_style="matchms", append=False)
    first_lines = output_path.read_text().splitlines()[:5]
    assert first_lines == [
        "BEGIN IONS",
        "CHARGE=1+",
        "RETENTION_TIME=92.382",
        "PRECURSOR_MZ=222.1597604",
        "100.251302248751 159.445 "
    ]

def test_export_gnps():
    input_path = Path("tests/resources/Alchorneine.mgf")
    output_path = Path("out.mgf")
    spectras = matchms.importing.load_from_mgf(str(input_path))
    (s, ) = spectras
    assert s.metadata_dict()["precursor_mz"] == 222.1597604
    output_path.unlink(missing_ok=True)
    matchms.exporting.save_spectra([s], str(output_path), export_style="gnps", append=False)
    first_lines = output_path.read_text().splitlines()[:4]
    assert first_lines == [
        "BEGIN IONS",
        "CHARGE=1",
        "PEPMASS=222.1597604",
        "100.251302248751 159.445 "
    ]

def test_export_ploum():
    input_path = Path("tests/resources/Alchorneine.mgf")
    output_path = Path("out.mgf")
    spectras = matchms.importing.load_from_mgf(str(input_path))
    (s, ) = spectras
    assert s.metadata_dict()["precursor_mz"] == 222.1597604
    s.set("ploum", "ploum42")
    output_path.unlink(missing_ok=True)
    matchms.exporting.save_spectra([s], str(output_path), export_style="matchms", append=False)
    first_lines = output_path.read_text().splitlines()[:6]
    assert first_lines == [
        "BEGIN IONS",
        "CHARGE=1+",
        "RETENTION_TIME=92.382",
        "PRECURSOR_MZ=222.1597604",
        "PLOUM=ploum42",
        "100.251302248751 159.445 "
    ]

def test_export_scans():
    input_path = Path("tests/resources/Alchorneine.mgf")
    output_path = Path("out.mgf")
    spectras = matchms.importing.load_from_mgf(str(input_path))
    (s, ) = spectras
    assert s.metadata_dict()["precursor_mz"] == 222.1597604
    s.set("scans", "ploum42")
    output_path.unlink(missing_ok=True)
    matchms.exporting.save_spectra([s], str(output_path), export_style="matchms", append=False)
    first_lines = output_path.read_text().splitlines()[:6]
    assert first_lines == [
        "BEGIN IONS",
        "CHARGE=1+",
        "RETENTION_TIME=92.382",
        "PRECURSOR_MZ=222.1597604",
        "SCANS=ploum42",
        "100.251302248751 159.445 "
    ]

def test_export_mslevel():
    input_path = Path("tests/resources/Alchorneine.mgf")
    output_path = Path("out.mgf")
    spectras = matchms.importing.load_from_mgf(str(input_path))
    (s, ) = spectras
    assert s.metadata_dict()["precursor_mz"] == 222.1597604
    s.set("mslevel", "ploum42")
    output_path.unlink(missing_ok=True)
    matchms.exporting.save_spectra([s], str(output_path), export_style="matchms", append=False)
    first_lines = output_path.read_text().splitlines()[:6]
    assert first_lines == [
        "BEGIN IONS",
        "CHARGE=1+",
        "RETENTION_TIME=92.382",
        "PRECURSOR_MZ=222.1597604",
        "MS_LEVEL=ploum42",
        "100.251302248751 159.445 "
    ]

def test_export_rtinseconds_fails():
    input_path = Path("tests/resources/Alchorneine.mgf")
    output_path = Path("out.mgf")
    spectras = matchms.importing.load_from_mgf(str(input_path))
    (s, ) = spectras
    assert s.metadata_dict()["precursor_mz"] == 222.1597604
    with pytest.raises(ValueError) as e_info:
        s.set("rtinseconds", "ploum42")
    with pytest.raises(ValueError) as e_info:
        s.set("RTINSECONDS", "ploum42")
    with pytest.raises(ValueError) as e_info:
        s.set("rtinseconds", "92")
    with pytest.raises(ValueError) as e_info:
        s.set("rtinseconds", 92)
    with pytest.raises(ValueError) as e_info:
        s.set("RTINSECONDS", "92")
    with pytest.raises(ValueError) as e_info:
        s.set("RTINSECONDS", 92)

def test_export_manual():
    input_path = Path("tests/resources/Alchorneine.mgf")
    output_path = Path("out.mgf")
    spectras = matchms.importing.load_from_mgf(str(input_path))
    (spectrum, ) = spectras
    assert spectrum.metadata_dict()["precursor_mz"] == 222.1597604
    metadata_copy = matchms.Metadata(spectrum.metadata_dict(), matchms_key_style=False)
    metadata_expanded = metadata_copy.set("ploum", "ploum42")
    assert metadata_expanded.data.keys() == {"charge", "retention_time", "precursor_mz", "ploum"}
    spectrum_copy = matchms.Spectrum(spectrum.mz, spectrum.intensities, metadata_expanded.data, metadata_harmonization=False)
    output_path.unlink(missing_ok=True)
    matchms.exporting.save_spectra([spectrum_copy], str(output_path), export_style="matchms", append=False)
    first_lines = output_path.read_text().splitlines()[:6]
    assert first_lines == [
        "BEGIN IONS",
        "CHARGE=1+",
        "RETENTION_TIME=92.382",
        "PRECURSOR_MZ=222.1597604",
        "PLOUM=ploum42",
        "100.251302248751 159.445 "
    ]
    output_path.unlink(missing_ok=True)

def test_export_rtinseconds_manual_fails():
    input_path = Path("tests/resources/Alchorneine.mgf")
    output_path = Path("out.mgf")
    spectras = matchms.importing.load_from_mgf(str(input_path))
    (spectrum, ) = spectras
    assert spectrum.metadata_dict()["precursor_mz"] == 222.1597604
    metadata_copy = matchms.Metadata(spectrum.metadata_dict(), matchms_key_style=False)
    # A warning occurs.
    metadata_expanded = metadata_copy.set("rtinseconds", 42)
    assert metadata_expanded.data.keys() == {"charge", "retention_time", "precursor_mz", "rtinseconds"}
    spectrum_copy = matchms.Spectrum(spectrum.mz, spectrum.intensities, metadata_expanded.data, metadata_harmonization=False)
    output_path.unlink(missing_ok=True)
    matchms.exporting.save_spectra([spectrum_copy], str(output_path), export_style="matchms", append=False)
    first_lines = output_path.read_text().splitlines()[:5]
    assert first_lines == [
        "BEGIN IONS",
        "CHARGE=1+",
        "RETENTION_TIME=92.382",
        "PRECURSOR_MZ=222.1597604",
        "100.251302248751 159.445 "
    ]
    output_path.unlink(missing_ok=True)
