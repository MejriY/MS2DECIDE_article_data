from pathlib import Path
import os
import shutil
import pandas as pd
import re

REPO_DIR = Path("../")
INPUT_DIR = REPO_DIR / "Pleiokomenine C/"
GENERATED_DIR = REPO_DIR / "Generated/" / "Pleiokomenine C/"
GENERATED_DIR_CARTESIAN = GENERATED_DIR / "Cartesian coordinates/"
GENERATED_DIR_ARTICLE = REPO_DIR / "../" / "Article/" / "Generated/" / "PleiokomenineC/"

def transform_coordinates():
    GENERATED_DIR_CARTESIAN.mkdir(parents=True, exist_ok=True)

    input_dir = INPUT_DIR / "Cartesian coordinates/"
    inputs = os.listdir(input_dir)
    for i in inputs:
        assert i.endswith(".txt")
        nb = i[:-4]
        transform_coordinates_file(input_dir / f"{nb}.txt", GENERATED_DIR_CARTESIAN / f"{nb}.tsv")

def transform_coordinates_file(input_file, output_file):
    input = pd.read_fwf(input_file, skiprows = 4, names = ["Index", "Element", "x", "y", "z"], index_col = "Index")
    assert len(input) == 110
    input.to_csv(output_file, sep="\t")

def generate_article_data():
    GENERATED_DIR_ARTICLE.mkdir(parents=True, exist_ok=True)
    shutil.copytree(GENERATED_DIR_CARTESIAN, GENERATED_DIR_ARTICLE, dirs_exist_ok=True)
    conf_input = INPUT_DIR / "Boltzmann analysis/" / "Conformational population.tsv"
    pd.read_csv(conf_input, sep="\t").set_index("Conformation").rename(columns=lambda x: re.sub('[ ()/]', '', x)).to_csv(GENERATED_DIR_ARTICLE / "Conformationalpopulation.tsv", sep="\t")
    shutil.copy(INPUT_DIR / "Boltzmann analysis/" / "Eight conformations.tsv", GENERATED_DIR_ARTICLE / "Eightconformations.tsv")

