import extractor.manufactured.datadirs as mdirs
import shutil
import pandas as pd
from ms2decide.Tool import Tool
from ms2decide.K import K
from ms2decide.Tanimotos import Tanimotos


def generate_article_data():
    gen_article_dir = mdirs.REPO_DIR / "../" / "Article/" / "Generated/" / "Elicitation/"
    gen_article_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(mdirs.REPO_DIR / "Elicitation/" / "AnswersDual.tsv", gen_article_dir / "AnswersDual.tsv")
    cases = pd.read_csv(mdirs.REPO_DIR / "Elicitation/" / "Cases.tsv", sep="\t").set_index("Id")
    cases["K"] = cases.apply(lambda r: round(get_k(r), 12), axis=1)
    cases.to_csv(gen_article_dir / "Cases.tsv", sep="\t")


def get_k(row_series):
    similarities = {
        Tool.GNPS.name: row_series["cg"],
        Tool.SIRIUS.name: row_series["cs"],
        Tool.ISDB.name: row_series["ci"],
    }
    tanimotos = Tanimotos({})
    tanimotos.compute_tanimoto()
    tanimotos.tgs = row_series["tgs"]
    tanimotos.tgi = row_series["tgi"]
    tanimotos.tsi = row_series["tsi"]
    tanimotos.data = similarities
    return K(similarities, tanimotos).k()
