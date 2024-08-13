import os
import pandas as pd
from extractor.gnps import GnpsAnnotationsFile
from extractor.gnps import GnpsCacher
from extractor.gnps import GnpsParametersFile
from extractor.gnps import GnpsInchiScore

compounds_file = "../Manufactured case/Compounds.tsv"
compounds = pd.read_csv(compounds_file, sep="\t").set_index("Id")

p = GnpsCacher.cache_retrieve("26a5cbca3e844cc0b126f992c69df832")
x = GnpsCacher.cache_retrieve_parameters("26a5cbca3e844cc0b126f992c69df832")
isc = GnpsInchiScore(p, x)
new_cols = {f"inchi_gnps_{isc.min_peaks}_{isc.max_delta_mass}": isc.inchis, f"score_gnps_{isc.min_peaks}_{isc.max_delta_mass}": isc.scores}
# assert summary.index == compounds.index
compounds.assign(**new_cols)

summary = GnpsAnnotationsFile(p).summary()
ps = GnpsParametersFile(x).params()
print(ps)
min_peaks = int(ps["MIN_MATCHED_PEAKS"])
max_delta_mass = float(ps["tolerance.PM_tolerance"])
new_cols = {f"inchi_gnps_{min_peaks}_{max_delta_mass}": summary.loc[:, "INCHI"]}
# assert summary.index == compounds.index
compounds.assign(**new_cols)

def main():
    compounds_file = "../Manufactured case/Compounds.tsv"
    compounds = pd.read_csv(compounds_file, sep="\t").set_index("Id")

    p = GnpsCacher.cache_retrieve("26a5cbca3e844cc0b126f992c69df832")
    df = GnpsAnnotationsFile(p).summary()
    # print(df)
    x = GnpsCacher.cache_retrieve_parameters("26a5cbca3e844cc0b126f992c69df832")
    # x = GnpsFetcher.fetch_parameters_and_save("26a5cbca3e844cc0b126f992c69df832")
    ps = GnpsParametersFile(x).params()
    print(ps)
    min_peaks = int(ps["MIN_MATCHED_PEAKS"])
    max_delta_mass = float(ps["tolerance.PM_tolerance"])
    compounds[f"inchi_gnps_{min_peaks}_{max_delta_mass}"] = df.loc[:, "INCHI"]