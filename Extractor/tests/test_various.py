import pandas as pd
from zipfile import ZipFile 
import zipfile
from extractor.gnps import GnpsFile
import json

def test_read_tsv():
    path = "tests/resources/75b94088c0f24a0eb8b064b380300649.tsv"
    df = pd.read_csv(path, sep="\t")
    assert df.shape == (71, 46)

def test_extract():
    zip_file = "tests/resources/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-26a5cbca-view_all_analog_annotations_DB.zip"
    internal_path = "DB_result/75b94088c0f24a0eb8b064b380300649.tsv"
    with ZipFile(zip_file) as myzip:
        with myzip.open(internal_path) as myfile:
            df = pd.read_csv(myfile, sep="\t")
            assert df.shape == (71, 46)

def test_inside():
    zip_file = "tests/resources/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-26a5cbca-view_all_analog_annotations_DB.zip"
    df = GnpsFile(zip_file).df()
    assert df.shape == (71, 46)
    
def test_inside_manual():
    zip_file = "tests/resources/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-26a5cbca-view_all_analog_annotations_DB.zip"
    with ZipFile(zip_file) as myzip:
        p = zipfile.Path(myzip, at="DB_result/")
        l = set(p.iterdir())
        assert len(l) == 1
        (e, ) = l
        with e.open() as myfile:
            df = pd.read_csv(myfile, sep="\t")
            assert df.shape == (71, 46)

def test_read_df():
    json_file = "tests/resources/26a5cbca3e844cc0b126f992c69df832.json"
    with open(json_file) as json_data:
        js = json.load(json_data)
        assert len(js) == 1
        (k, v), = js.items()
        df = pd.DataFrame(v)
    assert df.shape == (71, 48)