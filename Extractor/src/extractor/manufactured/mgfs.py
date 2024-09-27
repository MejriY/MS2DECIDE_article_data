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
    
    def _level2(spectrum):
        l2 = matchms.Spectrum(spectrum.mz, spectrum.intensities, spectrum.metadata, metadata_harmonization=False)
        l2.set("ms_level", 2)
        # l2.set("COLLISION_ENERGY", 0)
        mzs = l2.mz
        mz_parent = l2.metadata_dict()["precursor_mz"]
        kept = mzs <= (mz_parent + 4)
        l2.set("num_peaks", sum(kept))
        kept_intensities = l2.intensities[kept]
        max_intensity = max(kept_intensities)
        normalized_intensities = l2.intensities[kept] / max_intensity * 100
        return matchms.Spectrum(l2.mz[kept], normalized_intensities, l2.metadata, metadata_harmonization=False)
    
    def _level1(spectrum):
        l1 = matchms.Spectrum(spectrum.mz, spectrum.intensities, spectrum.metadata, metadata_harmonization=False)
        mzs = l1.mz
        mz_parent = l1.metadata_dict()["precursor_mz"]
        kept = mzs >= mz_parent - 0.01
        l1.set("ms_level", 1)
        l1.set("num_peaks", sum(kept))
        kept_intensities = l1.intensities[kept]
        max_intensity = max(kept_intensities, default=0)
        normalized_intensities = l1.intensities[kept] / max_intensity * 100
        return matchms.Spectrum(l1.mz[kept], normalized_intensities, l1.metadata, metadata_harmonization=False)
    
    @cache
    def all_spectra_cut(self):
        all_spectra = list()
        by_id = self._by_name.reset_index().set_index("Id")
        for id in by_id.index:
            name = by_id.loc[id, "Chemical name"]
            spectrum = self.from_name(name)
            spectrum.set("scans", id)
            spectrum.set("feature_id", id)
            spectrum_copy = MgfFiles._level2(spectrum)
            assert spectrum_copy.mz.size <= spectrum.mz.size
            all_spectra.append(MgfFiles._level1(spectrum_copy))
            all_spectra.append(spectrum_copy)
        return all_spectra
    
    def export_all_spectra(self, path, export_style = "matchms"):
        # Delete file first as matchms appends to it.
        path.unlink(missing_ok=True)
        matchms.exporting.save_as_mgf(self.all_spectra(), str(path), export_style=export_style)
    
    def export_each_spectra_cut(self, dir):
        dir.mkdir(parents=True, exist_ok=True)
        by_id = self._by_name.reset_index().set_index("Id")
        for id in by_id.index:
            name = by_id.loc[id, "Chemical name"]
            spectrum = self.from_name(name)
            spectrum.set("scans", id)
            spectrum.set("feature_id", id)
            spectrum_l2 = MgfFiles._level2(spectrum)
            assert spectrum_l2.mz.size <= spectrum.mz.size
            path = dir / (name + ".mgf")
            path.unlink(missing_ok=True)
            matchms.exporting.save_as_mgf([MgfFiles._level1(spectrum_l2), spectrum_l2], str(path))
    
    def export_all_spectra_cut(self, path):
        # MZmine asks to use MASCOT generic format (https://mzmine.github.io/mzmine_documentation/module_docs/io/data-export.html) and https://www.matrixscience.com/help/data_file_help.html says: CHARGE=2+ (matchms) and PEPMASS (gnps) and RTINSECONDS (neither)…
        spectra = self.all_spectra_cut()
        MgfFiles.export_sirius(spectra, path)

    def export_sirius(spectra, path):
        path.unlink(missing_ok=True)
        matchms.exporting.save_as_mgf(spectra, str(path))
        patched = path.read_text().replace("PRECURSOR_MZ=", "PEPMASS=").replace("RETENTION_TIME=", "RTINSECONDS=").replace("MS_LEVEL=", "MSLEVEL=").replace("NUM_PEAKS=", "Num peaks=")
        path.write_text(patched)
