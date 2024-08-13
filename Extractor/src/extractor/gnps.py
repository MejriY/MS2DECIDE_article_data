import pandas as pd
from zipfile import ZipFile 
import zipfile
from typing import BinaryIO
import os
import requests

class GnpsFile:
    def __init__(self, compressed: str | os.PathLike | BinaryIO):
        self.compressed = compressed

    def df(self):
        with ZipFile(self.compressed) as myzip:
            p = zipfile.Path(myzip, at="DB_result/")
            l = set(p.iterdir())
            assert len(l) == 1
            (e, ) = l
            with e.open() as myfile:
                return pd.read_csv(myfile, sep="\t")
    
class GnpsFetcher:
    def fetch(task_id: str):
        # url = f"https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={task_id}&view=view_all_annotations_DB"
        url = f"https://gnps.ucsd.edu/ProteoSAFe/result_json.jsp?task={task_id}&view=view_all_annotations_DB"
        with requests.get(url) as r:
            r.raise_for_status()
            return r.content

    def fetch_and_save(task_id: str, path: str | os.PathLike):
        with open(path, "wb") as f:
            f.write(GnpsFetcher.fetch(task_id))
