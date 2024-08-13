import os
import pandas as pd
from extractor.gnps import GnpsFile
from extractor.gnps import GnpsFetcher

def main():
    # Create dir if does not exist
    os.makedirs("fetched", exist_ok=True)
    print("Downloading.")
    GnpsFetcher.fetch_and_save("26a5cbca3e844cc0b126f992c69df832", "fetched/26a5cbca3e844cc0b126f992c69df832.json")
    print("Downloaded.")

