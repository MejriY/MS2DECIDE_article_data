import pytest
from extractor.doier import Doier
from extractor.article import Article, ArticleDoi

article_AAA = Article(title="Automated grading systems for programming assignments: A literature review", authors=("H Aldriye", "A Alkhalaf", "M Alkhalaf"), authors_partial=False, year=2019, abstract="abstract", link="https://example.com/", PDF="https://pdfs.semanticscholar.org/3ec5/6c7bc0338717acca9d5e519c0a27ca5ce70d.pdf")
article_DMP = Article(title="Automatic question generation and answer assessment: a survey", authors=("B Das", "M Majumder", "S Phadikar"), authors_partial=False, year=2021, abstract="abstract", link="https://example.com/")
article_IAK = Article(title="Review of recent systems for automatic assessment of programming assignments", authors=("P Ihantola", "T Ahoniemi", "V Karavirta"), authors_partial=True, year=2010, abstract="abstract", link="https://example.com/", PDF="https://www.researchgate.net/profile/Petri-Ihantola/publication/216714976_Review_of_recent_systems_for_automatic_assessment_of_programming_assignments/links/0deec5182391eb2954000000/Review-of-recent-systems-for-automatic-assessment-of-programming-assignments.pdf")

@pytest.mark.slow
def test_AAA():
    d = Doier()
    art = d.query_and_merge(article_AAA)
    assert "10.14569/ijacsa.2019.0100328" == art.DOI

@pytest.mark.slow
def test_DMP():
    d = Doier()
    art = d.query_and_merge(article_DMP)
    assert "10.1186/s41039-021-00151-1" == art.DOI

@pytest.mark.slow
def test_IAK():
    d = Doier()
    art = d.query_and_merge(article_IAK)
    assert "10.1145/1930464.1930480" == art.DOI

@pytest.mark.slow
def test_none():
    d = Doier()
    art = d.query_article(Article(title="Does not exist", authors=("A", "B", "C"), authors_partial=False, year=2000, abstract="abstract", link=None))
    assert art is None
