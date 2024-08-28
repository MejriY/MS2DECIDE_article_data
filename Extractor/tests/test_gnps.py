from pathlib import Path
from extractor.gnps import GnpsAnnotations
from extractor.gnps import GnpsParameters
from extractor.gnps import GnpsMatch
from extractor.gnps import GnpsTaskFetcher
from extractor.gnps import GnpsCachingTaskFetcher
from extractor.gnps import GnpsInchiSmiles
from extractor.gnps import GnpsQuery
from extractor.gnps import GnpsQueried
from extractor.gnps import IteratedQueries
import responses
import json
from rdkit import Chem
from rdkit import RDLogger

def test_annot_empty():
    a = GnpsAnnotations()
    df = a.raw_df()
    assert df.shape == (0, 11)

def test_annot_empty_json():
    a = GnpsAnnotations("{ \"blockData\" : [] }")
    df = a.raw_df()
    assert df.shape == (0, 11)

def test_annot_summary():
    input_path = Path("tests/resources/1c3cbef069874546aca7cca96cb01a05.json")
    a = GnpsAnnotations.from_file(input_path)
    summary = a.summary_df()
    assert summary.shape == (95, 12)
    assert summary.iloc[0:5].index.to_list() == list(range(1, 6))

def test_annot():
    input_path = Path("tests/resources/26a5cbca3e844cc0b126f992c69df832.json")
    a = GnpsAnnotations.from_file(input_path)
    df = a.raw_df()
    assert df.shape == (71, 47)
    insms = a.inchis_smiles_series()
    assert insms.shape == (71, )
    insm = insms.iloc[0]
    assert type(insm) == GnpsInchiSmiles
    summary = a.summary_df()
    assert summary.shape == (71, 12)
    matches = a.matches_series()
    assert matches.shape == (71, )
    assert type(matches.iloc[0]) == GnpsMatch

def test_annot_encoding():
    input_path = Path("tests/resources/Invalid UTF-8 61f141ad8c3449fdb78a5086630b40f5.json")
    a = GnpsAnnotations.from_file(input_path)
    df = a.raw_df()
    assert df.shape == (521, 47)
    assert df.loc[1158, "Compound_Name"] == "(1S,5R,9S,13R)-14-formyl-5,9-dimethyltetracyclo[11.2.1.0�,�???.0???,???]hexadec-14-ene-5-carboxylic acid"

def test_parameters():
    input_path = Path("tests/resources/26a5cbca3e844cc0b126f992c69df832.xml")
    ps = GnpsParameters.from_file(input_path).all()
    assert len(ps) == 37

def test_fetch_live():
    task_id = "1c3cbef069874546aca7cca96cb01a05"
    assert GnpsTaskFetcher(task_id).done()
    # assert a.raw_df().shape == (95, 47)

@responses.activate
def test_caching_done(tmp_path):
    task_id = "1c3cbef069874546aca7cca96cb01a05"
    resp_status = responses.get(
        GnpsTaskFetcher.url_status(task_id),
        json=json.loads(Path("tests/resources/status 1c3cbef069874546aca7cca96cb01a05.json").read_text()),
    )
    resp_exact = responses.get(
        GnpsTaskFetcher.url_exact(task_id),
        json=json.loads(Path("tests/resources/1c3cbef069874546aca7cca96cb01a05.json").read_text()),
    )
    c = GnpsCachingTaskFetcher(tmp_path)
    a = c.exact_annotations(task_id)
    assert a.raw_df().shape == (95, 47)
    a2 = c.exact_annotations(task_id)
    assert a2.raw_df().shape == (95, 47)
    assert resp_status.call_count == 1
    assert resp_exact.call_count == 1

@responses.activate
def test_caching_not_done(tmp_path):
    task_id = "1c3cbef069874546aca7cca96cb01a05"
    resp_status = responses.get(
        GnpsTaskFetcher.url_status(task_id),
        json=json.loads(Path("tests/resources/status fake not done.json").read_text()),
    )
    c = GnpsCachingTaskFetcher(tmp_path)
    a = c.exact_annotations(task_id)
    assert a.raw_df().shape == (0, 11)
    a2 = c.exact_annotations(task_id)
    assert a2.raw_df().shape == (0, 11)
    assert resp_status.call_count == 2

def test_inchi_smiles():
    RDLogger.DisableLog("rdApp.*")
    correct_inchi = "InChI=1S/C9H6O3/c10-7-3-1-6-2-4-9(11)12-8(6)5-7/h1-5,10H"
    incorrect_inchi = "1S/C15H22O4/c1-9-5-6-13-10(2)7-11(16)14(3,4)15(13)8-12(17)18/h5,10,13H,6-8H2,1-4H3,(H,17,18)/t10-,13-/m1/s1"
    s = GnpsInchiSmiles(f'"{incorrect_inchi}"', "O=C1OC=2C=C(O)C=CC2C=C1")
    assert s.raw_inchi == f'"{incorrect_inchi}"'
    assert s.sanitized_inchi == f"InChI={incorrect_inchi}"
    assert Chem.inchi.MolFromInchi(s.sanitized_inchi) is None
    assert s.to_mol is not None
    assert s.to_mol.GetNumAtoms() == 12
    assert s.to_standard_inchi == correct_inchi

def test_order_attempts():
    a61 = GnpsQuery(6, 1)
    a62 = GnpsQuery(6, 2)
    a41 = GnpsQuery(4, 1)
    a42 = GnpsQuery(4, 2)
    assert a61 < a62
    assert a62 < a41
    assert a41 < a42
    assert sorted([a41, a62, a42, a61]) == [a61, a62, a41, a42]

def test_queried():
    p = GnpsParameters.from_file(Path("tests/resources/1c3cbef069874546aca7cca96cb01a05.xml"))
    a1 = GnpsAnnotations.from_file(Path("tests/resources/1c3cbef069874546aca7cca96cb01a05.json"))
    a2 = GnpsAnnotations.from_file(Path("tests/resources/1c3cbef069874546aca7cca96cb01a05.json"))

    q = GnpsQueried(p.to_query(), a1, a2)
    assert q.query == p.to_query()
    assert q.summary_df().shape == (95, 10)

def test_iterated():
    p1 = GnpsParameters.from_file(Path("tests/resources/1c3cbef069874546aca7cca96cb01a05.xml"))
    a1 = GnpsAnnotations.from_file(Path("tests/resources/1c3cbef069874546aca7cca96cb01a05.json"))
    q1 = GnpsQueried(p1.to_query(), a1, a1)
    p2 = GnpsParameters.from_file(Path("tests/resources/4dc38f3e2d0a4536b2ddb91d6526fca7.xml"))
    a2 = GnpsAnnotations.from_file(Path("tests/resources/4dc38f3e2d0a4536b2ddb91d6526fca7.json"))
    q2 = GnpsQueried(p2.to_query(), a2, a2)
    i = IteratedQueries([q1, q2])
    # b = i._best_matches_discounted_series()
    # print(b)
    # print(i.all_df())
    assert i.all_df().shape == (95, 22)

def test_gnps_non_monotonous():
    small_delta_a = GnpsAnnotations.from_file("tests/resources/1c3cbef069874546aca7cca96cb01a05.json")
    small_delta_params = GnpsParameters.from_file("tests/resources/1c3cbef069874546aca7cca96cb01a05.xml")
    big_delta_a = GnpsAnnotations.from_file("tests/resources/4dc38f3e2d0a4536b2ddb91d6526fca7.json")
    big_delta_params = GnpsParameters.from_file("tests/resources/4dc38f3e2d0a4536b2ddb91d6526fca7.xml")
    # small_delta = GnpsInchiScore(small_delta_a, small_delta_params)
    # big_delta = GnpsInchiScore(big_delta_a, big_delta_params)
    # to be written: id 23 was found with the stricter parameters but not with larger ones.
