import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

class Compounds:
    def __init__(self, compounds):
        self.df = compounds
    
    @classmethod
    def from_tsv(cls, compounds_file):
        return cls(pd.read_csv(compounds_file, sep="\t", dtype={"Id": int}).set_index("Id"))
    
    def add_precursors(self, precursors):
        self.df["Precursor m/z"] = precursors
    
    def add_retention_times(self, retention_seconds):
        self.df["Retention time (sec)"] = retention_seconds
    
    def add_relative_molecular_weights(self):
        self.df["Relative molecular weight"] = self.df["InChI"].apply(lambda i: Descriptors.MolWt(Chem.inchi.MolFromInchi(i)) if i is not None else None)

    def add_diffs(self):
        self.df["Precursor m/z − relative molecular weight"] = (
        self.df["Precursor m/z"] - self.df["Relative molecular weight"]
    )
        
    def quantification_table_minutes(self):
        qt = pd.DataFrame(index=self._compounds.index).rename_axis("row ID")
        qt["row m/z"] = self._compounds["Precursor m/z"]
        qt["row retention time"] = self._compounds["Retention time (sec)"] / 60.0
        qt["1.mzXML Peak area"] = 0
        return qt
    