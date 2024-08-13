import pandas as pd
from extractor.gnps import GnpsFile
import json
import extractor.gnps as gnps

json_file = "tests/resources/26a5cbca3e844cc0b126f992c69df832.json"
with open(json_file) as json_data:
    js = json.load(json_data)
    assert len(js) == 1
    (k, v), = js.items()
    dfj = pd.DataFrame(v)
assert dfj.shape == (71, 48)
zip_file = "tests/resources/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-26a5cbca-view_all_analog_annotations_DB.zip"
internal_path = "DB_result/75b94088c0f24a0eb8b064b380300649.tsv"
with ZipFile(zip_file) as myzip:
    with myzip.open(internal_path) as myfile:
        dfz = pd.read_csv(myfile, sep="\t")
        assert dfz.shape == (71, 46)
cj = dfj.columns
cz = dfz.columns
# difference?
print(set(cj) - set(cz))

def main():
    # Create dir if does not exist
    os.makedirs("fetched", exist_ok=True)
    print("Downloading.")
    gnps.GnpsFetcher.fetch_and_save("26a5cbca3e844cc0b126f992c69df832", "fetched/26a5cbca3e844cc0b126f992c69df832")
    print("Downloaded.")

