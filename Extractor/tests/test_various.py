import pandas as pd
from extractor.gnps import GnpsAnnotations
from extractor.gnps import GnpsParametersFile
import json
import xml.etree.ElementTree as ET

def test_read_df():
    df = GnpsAnnotations.from_file("tests/resources/26a5cbca3e844cc0b126f992c69df832.json").df()
    assert df.shape == (71, 48)
    
def test_read_df_manual():
    json_file = "tests/resources/26a5cbca3e844cc0b126f992c69df832.json"
    with open(json_file) as json_data:
        js = json.load(json_data)
        assert len(js) == 1
        (k, v), = js.items()
        df = pd.DataFrame(v)
    assert df.shape == (71, 48)

def test_parse_parameters():
    xml_file = "tests/resources/26a5cbca3e844cc0b126f992c69df832.xml"
    ps = GnpsParametersFile(xml_file).params()
    assert len(ps) == 37

def test_parse_parameters_manual():
    xml_file = "tests/resources/26a5cbca3e844cc0b126f992c69df832.xml"
    with open(xml_file) as f:
        xml = f.read()
        tree = ET.parse(xml_file)
        root = tree.getroot()
        p0 = root[0]
        assert p0.tag == "parameter"
        assert p0.attrib["name"] == "ANALOG_SEARCH"

