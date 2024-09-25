from pathlib import Path
import matchms
import pandas as pd
from functools import cache

class MgfFiles:
    def __init__(self, dir: Path, name_id_df: pd.DataFrame | pd.Series):
        self.dir = dir
        self.files = set(dir.glob("*.mgf"))
        self._sp_dict = {f.stem: self._spectrum(f) for f in self.files}
        name_id_df_2cols = name_id_df.reset_index()
        assert set(name_id_df_2cols.columns) == {"Chemical name", "Id"}
        self._by_name = name_id_df_2cols.set_index("Chemical name")
        names = name_id_df_2cols["Chemical name"].to_list()
        assert self.names == set(names), set(names) - set(self.names)

    def _spectrum(self, file):
        ss = list(matchms.importing.load_from_mgf(str(file)))
        assert len(ss) == 1
        return ss[0]
    
    def from_name(self, name):
        return self._sp_dict.get(name)
    
    @property
    def names(self):
        return self._sp_dict.keys()
    
    @property
    def d(self):
        return self._sp_dict
    
    def precursors_series(self):
        return pd.Series({self._by_name.at[n, "Id"]: s.get("precursor_mz") for n, s in self.d.items()})
    
    def retentions_seconds_series(self):
        return pd.Series({self._by_name.at[n, "Id"]: s.get("retention_time") for n, s in self.d.items()})
    
    @cache
    def all_spectra(self):
        all_spectra = list()
        by_id = self._by_name.reset_index().set_index("Id")
        for id in by_id.index:
            name = by_id.loc[id, "Chemical name"]
            spectrum = self.from_name(name)
            spectrum.set("scans", id)
            # No apparent effect when exported; seems that we need to build the spectrum using Spectrum(mz=sp.peaks.mz,intensities=sp.peaks.intensities,metadata=m).
            # spectrum.set("MSLEVEL", 2)
            all_spectra.append(spectrum)
        return all_spectra
    
    @cache
    def all_spectra_cut(self):
        all_spectra = list()
        by_id = self._by_name.reset_index().set_index("Id")
        for id in by_id.index:
            name = by_id.loc[id, "Chemical name"]
            spectrum = self.from_name(name)
            spectrum.set("scans", id)
            mzs = spectrum.mz
            kept = mzs <= (spectrum.metadata_dict()["precursor_mz"] + 4)
            spectrum_copy = matchms.Spectrum(spectrum.mz[kept], spectrum.intensities[kept], spectrum.metadata, metadata_harmonization=False)
            assert spectrum_copy.mz.size <= spectrum.mz.size
            all_spectra.append(spectrum_copy)
        return all_spectra
    
    def export_all_spectra(self, path, export_style = "matchms"):
        # Delete file first as matchms appends to it.
        path.unlink(missing_ok=True)
        matchms.exporting.save_as_mgf(self.all_spectra(), str(path), export_style=export_style)
    
    def export_all_spectra_cut(self, path, export_style = "matchms"):
        path.unlink(missing_ok=True)
        matchms.exporting.save_as_mgf(self.all_spectra_cut(), str(path), export_style=export_style)