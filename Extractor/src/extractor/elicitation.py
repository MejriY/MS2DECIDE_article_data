import extractor.manufactured.datadirs as mdirs
import shutil

def generate_article_data():
    gen_article_dir = mdirs.REPO_DIR / "../" / "Article/" / "Generated/" / "Elicitation/"
    gen_article_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(mdirs.REPO_DIR / "Elicitation/" / "AnswersDual.tsv", gen_article_dir / "AnswersDual.tsv")
    shutil.copyfile(mdirs.REPO_DIR / "Elicitation/" / "Cases.tsv", gen_article_dir / "Cases.tsv")
