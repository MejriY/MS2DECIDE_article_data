import os
import pandas as pd
import json
import matchms
from pathlib import Path
import py4cytoscape as p4c
from zipfile import ZipFile
import io
import fs
from IPython.display import SVG, display
import networkx as nx
from functools import cached_property
from functools import cache

# p4c.set_summary_logger(False)
# p4c._SUMMARY_LOG_LEVEL = "NOTSET"
# p4c.cytoscape_ping()
# zip_file = Path("Data") / "Pleiocarpa case" / "8 - Cytoscape" / "network_k.cys"
# s = p4c.open_session(str(zip_file))
# dir = "../../../Local/"
# allsvg = dir + "All.pdf"
# p4c.export_image(filename=allsvg, type="PDF", overwrite_file=True)
# # display(SVG(filename=allsvg))
# nets = p4c.get_network_list()
# assert len(nets) == 1
# net = nets[0]
# p4c.get_edge_count(net)
# es = p4c.get_all_edges(net)
# netx = p4c.create_networkx_from_network()
# netx
# import networkx as nx

# netx = nx.Graph(netx)
# cons = sorted(nx.connected_components(netx), key=len, reverse=True)
# [len(c) for c in cons]
# con = cons[0]
# len(netx.nodes)
# sorted(netx.nodes)

# e = es[0]
# type(e)
# cols = p4c.get_table_columns()
# cols["name"].sort_values()
# p4c.get_node_count(net)
# ns = p4c.get_all_nodes(net)
# p4c.clear_selection()
# n = ns[0]
# p4c.node_name_to_node_suid(n)


# p4c.select_nodes(48129)
# p4c.get_selected_node_count()

# p4c.fit_content(selected_only=True)
# p4c.export_image(filename=dir + "tmp.pdf", type="PDF", overwrite_file=True)

class Network:
    def __init__(self, cys_file, compounds: pd.DataFrame):
        self.session = p4c.open_session(str(cys_file))
        self.compounds = compounds
    
    @cache
    def netx(self):
        return nx.Graph(p4c.create_networkx_from_network())
    
    @cache
    def cytoscape_df(self):
        cdf = p4c.get_table_columns().set_index("SUID")
        cdf = cdf[cdf["ranking by k"].notna()].astype(dtype={"ranking by k": int})
        assert (cdf["parent mass"] == cdf["precursor mass"]).all()
        assert (cdf["shared name"] == cdf["name"]).all()
        assert (cdf["G1"] == 0).all()
        cdf = cdf.drop(columns=["parent mass", "shared name", "selected", "G1"])
        assert cdf["name"].duplicated().sum() == 0
        assert cdf.index.duplicated().sum() == 0
        return cdf
    
    @cache
    def _connected_component_names(self):
        netx = self.netx()
        cons = sorted(nx.connected_components(netx), key=len, reverse=True)
        return cons
    
    def _export_suids(self, suids, output_file):
        given_extension = output_file.split(".")[-1]
        p4c.select_nodes(suids, preserve_current_selection=False)
        p4c.fit_content(selected_only=True)
        p4c.select_nodes(None, preserve_current_selection=False)
        p4c.export_image(filename=output_file, type=given_extension.upper(), overwrite_file=True)

    def export(self):
        cons = self._connected_component_names()
        for i, con in enumerate(cons[:5]):
            suids = p4c.node_name_to_node_suid(list(con))
            self._export_suids(suids, f"Com {i}.svg")
        self._export_suids(self.cytoscape_df().index.to_list(), "All.svg")
        self._export_suids(self.cytoscape_df().index.to_list(), "All.pdf")