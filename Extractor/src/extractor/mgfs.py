from pathlib import Path
import matchms
import pandas as pd

class MgfFiles:
    def __init__(self, dir: Path):
        self.dir = dir
        self.files = set(dir.glob("*.mgf"))
        self._sp_dict = {f.stem: self._spectrum(f) for f in self.files}

    def _spectrum(self, file):
        ss = list(matchms.importing.load_from_mgf(str(file)))
        assert len(ss) == 1
        return ss[0]
    
    def by_name(self, name):
        return self._sp_dict.get(name)
    
    @property
    def prefix_file_names(self):
        return self._sp_dict.keys()
    
    @property
    def d(self):
        return self._sp_dict
    
    @property
    def precursors(self):
        return pd.Series({n: s.get("precursor_mz") for n, s in self.d.items()})
    
    @property
    def retentions_sec(self):
        return pd.Series({n: s.get("retention_time") for n, s in self.d.items()})