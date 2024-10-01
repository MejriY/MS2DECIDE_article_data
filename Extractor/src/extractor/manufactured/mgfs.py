from pathlib import Path
import matchms
import pandas as pd
from functools import cache

class MgfFiles:
    def __init__(self, dir: Path):
        self.dir = dir
        self.files = set(dir.glob("*.mgf"))
        self._sp_dicts_dict = {f.stem: self._spectra_dict(f) for f in self.files}
        names = self._sp_dicts_dict.keys()
        self._by_name = name_id_df(names)

    def _spectra_dict(self, file):
        ss = list(matchms.importing.load_from_mgf(str(file)))
        assert len(ss) == 2
        return {MgfFiles._mslevel(s): s for s in ss}
    
    def _value(self, name, metadata_key):
        spectras = self.spectra_from_name(name)
        values = set([s.get(metadata_key) for s in spectras])
        assert len(values) == 1, values
        (value, ) = values
        return value

    def _mslevel(spectrum):
        str_mslevel = spectrum.get("ms_level")
        assert str_mslevel is not None, spectrum.metadata
        return int(str_mslevel)
    
    def spectra_dict_from_name(self, name):
        return self._sp_dicts_dict.get(name)
    
    def spectra_from_name(self, name):
        return self._sp_dicts_dict.get(name).values()
    
    @property
    def names(self):
        return self._sp_dicts_dict.keys()
    
    @property
    def dicts_dict(self):
        return self._sp_dicts_dict
    
    def precursors_series(self):
        return pd.Series({self._by_name.at[n, "Id"]: self._value(n, "precursor_mz") for n in self.names})
    
    def retentions_seconds_series(self):
        return pd.Series({self._by_name.at[n, "Id"]: self._value(n, "retention_time") for n in self.names})
    
    def all_ids(self):
        return self._by_name["Id"].to_list()
    
    def all_spectra(self, levels = [1, 2]):
        return self.spectra(self.all_ids(), levels)
    
    def spectra(self, ids, levels = [1, 2]):
        all_spectra = list()
        by_id = self._by_name.reset_index().set_index("Id")
        for id in ids:
            name = by_id.loc[id, "Chemical name"]
            spectrum_dict = self.spectra_dict_from_name(name)
            for level in levels:
                spectrum = spectrum_dict.get(level)
                assert spectrum is not None, name
                spectrum.set("scans", id)
                spectrum.set("feature_id", id)
                all_spectra.append(spectrum)
        return all_spectra
    
    def export_all_level2(self, path, export_style = "matchms"):
        self.export_level2(self.all_ids(), path, export_style)
    
    def export_level2(self, ids, path, export_style = "matchms"):
        path.unlink(missing_ok=True)
        matchms.exporting.save_spectra(self.spectra(ids, levels = [2]), str(path), export_style=export_style, append=False)

    def export_all_sirius(self, path):
        MgfFiles.export_sirius(self.all_spectra(), path)
    
    def export_sirius(spectra, path):
        # MZmine asks to use MASCOT generic format (https://mzmine.github.io/mzmine_documentation/module_docs/io/data-export.html) and https://www.matrixscience.com/help/data_file_help.html says: CHARGE=2+ (matchms) and PEPMASS (gnps) and RTINSECONDS (neither)…
        path.unlink(missing_ok=True)
        matchms.exporting.save_as_mgf(spectra, str(path))
        patched = path.read_text().replace("PRECURSOR_MZ=", "PEPMASS=").replace("RETENTION_TIME=", "RTINSECONDS=").replace("MS_LEVEL=", "MSLEVEL=").replace("NUM_PEAKS=", "Num peaks=")
        path.write_text(patched)

def name_id_df(names):
    assert len(names) == len(set(names))
    sorted_names = sorted(names, key=lambda n: (n.startswith("Unreported "), n))
    return pd.DataFrame({"Chemical name": sorted_names, "Id": range(1, len(sorted_names) + 1)}).set_index("Chemical name")
