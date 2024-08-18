import pandas as pd
from extractor.gnps import GnpsAnnotations
from extractor.gnps import GnpsParametersFile
from extractor.gnps import GnpsInchiScore
from extractor import support
from pathlib import PurePath
import json
import xml.etree.ElementTree as ET
import requests
from github import Github
from github import Auth
import urllib
import base64
from github import ContentFile

def test_read_df():
    df = GnpsAnnotations.from_file("tests/resources/26a5cbca3e844cc0b126f992c69df832.json").df()
    assert df.shape == (71, 48)
    
def test_read_df_manual():
    json_file = "tests/resources/26a5cbca3e844cc0b126f992c69df832.json"
    with open(json_file) as json_data:
        js = json.load(json_data)
        assert len(js) == 1
        (k, v), = js.items()
        df = pd.DataFrame(v)
    assert df.shape == (71, 48)

def test_parse_parameters():
    xml_file = "tests/resources/26a5cbca3e844cc0b126f992c69df832.xml"
    ps = GnpsParametersFile(xml_file).params()
    assert len(ps) == 37

def test_parse_parameters_manual():
    xml_file = "tests/resources/26a5cbca3e844cc0b126f992c69df832.xml"
    with open(xml_file) as f:
        xml = f.read()
        tree = ET.parse(xml_file)
        root = tree.getroot()
        p0 = root[0]
        assert p0.tag == "parameter"
        assert p0.attrib["name"] == "ANALOG_SEARCH"

def test_gnps_non_monotonous():
    small_delta_a = GnpsAnnotations.from_file("tests/resources/1c3cbef069874546aca7cca96cb01a05.json")
    small_delta_params = GnpsParametersFile("tests/resources/1c3cbef069874546aca7cca96cb01a05.xml")
    big_delta_a = GnpsAnnotations.from_file("tests/resources/4dc38f3e2d0a4536b2ddb91d6526fca7.json")
    big_delta_params = GnpsParametersFile("tests/resources/4dc38f3e2d0a4536b2ddb91d6526fca7.xml")
    small_delta = GnpsInchiScore(small_delta_a, small_delta_params)
    big_delta = GnpsInchiScore(big_delta_a, big_delta_params)
    # to be written: id 23 was found with the stricter parameters but not with larger ones.

def test_direct_dl_fails():
    credsd = support.creds(PurePath("."))
    sha = "7ff8db2500f278770845be3fa70283242d1054ac"
    url = f"https://github.com/MejriY/MS2DECIDE_article_data/raw/{sha}/Generated/Manufactured%20case/output%20gnps%20with%20new%20file/Manufactured%20annotation.tsv"
    answer = requests.get(url, auth=credsd)
    assert answer.status_code == 404, answer.status_code

def test_list_root():
    auth = Auth.Token(support.creds(PurePath("."))[1])
    with Github(auth=auth) as g:
        repo = g.get_repo("MejriY/MS2DECIDE_article_data")
        sha = "7ff8db2500f278770845be3fa70283242d1054ac"
        server_path = "/"
        contents = repo.get_contents(urllib.parse.quote(server_path), ref=sha)
        assert len(contents) == 6, contents

def test_list_dir():
    auth = Auth.Token(support.creds(PurePath("."))[1])
    with Github(auth=auth) as g:
        repo = g.get_repo("MejriY/MS2DECIDE_article_data")
        sha = "7ff8db2500f278770845be3fa70283242d1054ac"
        server_path = "/Generated"
        contents = repo.get_contents(urllib.parse.quote(server_path), ref=sha)
        assert len(contents) == 1, contents

def test_list_subdir():
    auth = Auth.Token(support.creds(PurePath("."))[1])
    with Github(auth=auth) as g:
        repo = g.get_repo("MejriY/MS2DECIDE_article_data")
        sha = "7ff8db2500f278770845be3fa70283242d1054ac"
        server_path = "/Generated/Manufactured case/output gnps with new file"
        contents = repo.get_contents(server_path, ref=sha)
        assert len(contents) == 2, contents

def test_dl():
    auth = Auth.Token(support.creds(PurePath("."))[1])
    with Github(auth=auth) as g:
        repo = g.get_repo("MejriY/MS2DECIDE_article_data")
        sha = "7ff8db2500f278770845be3fa70283242d1054ac"
        server_path = "/Generated/Manufactured case/output gnps with new file/Manufactured annotation.tsv"
        contents : ContentFile = repo.get_contents(server_path, ref=sha)
        assert contents.name == "Manufactured annotation.tsv", contents
        content_encoded = contents.content
        content_decoded = contents.decoded_content.decode("utf-8")
        assert content_decoded.startswith("ID\t"), content_decoded