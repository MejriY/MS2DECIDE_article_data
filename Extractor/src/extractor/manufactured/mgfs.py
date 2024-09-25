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
    def all_spectra_sirius(self):
        all_spectra = list()
        by_id = self._by_name.reset_index().set_index("Id")
        for id in by_id.index:
            name = by_id.loc[id, "Chemical name"]
            spectrum = self.from_name(name)
            spectrum.set("scans", id)
            # No apparent effect when exported; seems that we need to build the spectrum using Spectrum(mz=sp.peaks.mz,intensities=sp.peaks.intensities,metadata=m).
            matchms.Spectrum()
            spectrum.set("MSLEVEL", 2)
            print(spectrum.metadata)
            metadata_copy = matchms.Metadata(spectrum.metadata_dict(), False)
            metadata_expanded = metadata_copy.set("RTINSECONDS2", 3)
            spectrum_copy = matchms.Spectrum(spectrum.mz, spectrum.intensities, metadata_expanded.data, False)
            print(metadata_expanded.data)
            print(spectrum_copy.metadata)
            all_spectra.append(spectrum)
        return all_spectra
    
    def export_all_spectra(self, path, export_style = "matchms"):
        # Delete file first as matchms appends to it.
        path.unlink(missing_ok=True)
        matchms.exporting.save_as_mgf(self.all_spectra(), str(path), export_style=export_style)

    def export_all_spectra_sirius(self, path):
        # Delete file first as matchms appends to it.
        path.unlink(missing_ok=True)
        matchms.exporting.save_spectra(self.all_spectra_sirius(), str(path), export_style="matchms", append=False)