import os
import pandas as pd
from extractor.gnps import GnpsFile
from extractor.gnps import GnpsCacher

def main():
    p = GnpsCacher.cache_retrieve("26a5cbca3e844cc0b126f992c69df832")
    df = GnpsFile(p).df()
    print(df)
