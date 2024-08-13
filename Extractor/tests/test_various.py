import pandas as pd
from extractor.gnps import GnpsFile
import json

def test_read_df():
    df = GnpsFile("tests/resources/26a5cbca3e844cc0b126f992c69df832.json").df()
    assert df.shape == (71, 48)
    
def test_read_df_manual():
    json_file = "tests/resources/26a5cbca3e844cc0b126f992c69df832.json"
    with open(json_file) as json_data:
        js = json.load(json_data)
        assert len(js) == 1
        (k, v), = js.items()
        df = pd.DataFrame(v)
    assert df.shape == (71, 48)