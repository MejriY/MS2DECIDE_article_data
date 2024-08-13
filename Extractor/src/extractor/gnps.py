import pandas as pd
import os
import requests
import json

class GnpsFile:
    def __init__(self, json_file: str | os.PathLike):
        self.json_file = json_file

    def df(self):
        with open(self.json_file) as json_data:
            js = json.load(json_data)
            assert len(js) == 1
            (k, v), = js.items()
            df = pd.DataFrame(v)
            return df
    
class GnpsFetcher:
    def fetch(task_id: str):
        url = f"https://gnps.ucsd.edu/ProteoSAFe/result_json.jsp?task={task_id}&view=view_all_annotations_DB"
        with requests.get(url) as r:
            r.raise_for_status()
            return r.content

    def fetch_and_save(task_id: str, path: str | os.PathLike):
        with open(path, "wb") as f:
            f.write(GnpsFetcher.fetch(task_id))
