from pathlib import Path
import os
import pandas as pd

REPO_DIR = Path("../")
INPUT_DIR = REPO_DIR / "Pleiokomenine C/"
GENERATED_DIR = REPO_DIR / "Generated/" / "Pleiokomenine C/"

def transform_coordinates():
    output_dir = GENERATED_DIR / "Cartesian coordinates/"
    output_dir.mkdir(parents=True, exist_ok=True)

    input_dir = INPUT_DIR / "Cartesian coordinates/"
    inputs = os.listdir(input_dir)
    for i in inputs:
        assert i.endswith(".txt")
        nb = i[:-4]
        transform_coordinates_file(input_dir / f"{nb}.txt", output_dir / f"{nb}.tsv")

def transform_coordinates_file(input_file, output_file):
    input = pd.read_fwf(input_file, skiprows = 4, names = ["Index", "Element", "x", "y", "z"], index_col = "Index")
    assert len(input) == 110
    print(input)
    input.to_csv(output_file, sep="\t")
