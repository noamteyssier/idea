import logging
from typing import Optional
import numpy as np
import pandas as pd
import networkx as nx
from pyvis.network import Network
from ._utils import _array_to_hex


class IDEA:
    def __init__(
        self,
        degs: pd.DataFrame,
        go: pd.DataFrame,
        deg_size_name: str = "padj",
        deg_gene_name: str = "gene",
        go_size_name: str = "adj_pvalue",
        go_term_name: str = "term_name",
        go_overlap_name: str = "overlapping_genes",
        neg_log_xform_degs: bool = True,
        neg_log_xform_go: bool = True,
        set_deg_mass: bool = True,
        set_go_mass: bool = True,
        edge_color: str = "grey",
        gene_palette: str = "Reds",
        term_palette: str = "Blues",
        gene_color: Optional[str] = None,
        term_color: Optional[str] = None,
    ):
        """
        Initialize the IDEA class.

        Parameters
        ----------
        degs : pd.DataFrame
            A dataframe of differentially expressed genes.
        go : pd.DataFrame
            A dataframe of gene ontology terms.
        deg_size_name : str, optional
            The name of the column in `degs` that should be used for
            the size of the differentially expressed genes. By default,
            this is `"padj"`.
        deg_gene_name : str, optional
            The name of the column in `degs` that should be used for
            the gene names of the differentially expressed genes. By
            default, this is `"gene"`.
        go_size_name : str, optional
            The name of the column in `go` that should be used for
            the size of the gene ontology terms. By default, this is
            `"adj_pvalue"`.
        go_term_name : str, optional
            The name of the column in `go` that should be used for
            the gene ontology terms. By default, this is `"term_name"`.
        go_overlap_name : str, optional
            The name of the column in `go` that should be used for
            the overlapping genes of the gene ontology terms. By
            default, this is `"overlapping_genes"`.
        neg_log_xform_degs : bool, optional
            Whether or not to invert the attributes of the differentially
            expressed genes. By default, this is `True`.
        neg_log_xform_go : bool, optional
            Whether or not to invert the attributes of the gene ontology
            terms. By default, this is `True`.
        set_deg_mass : bool, optional
            Whether or not to set the mass of the differentially expressed
            genes with the deg size. By default, this is `True`.
        set_go_mass : bool, optional
            Whether or not to set the mass of the gene ontology terms
            with the go size. By default, this is `True`.
        edge_color : str, optional
            The color of the edges in the bipartite graph. By default,
            this is `"black"`.
        term_palette : str, optional
            The color palette to use for the gene ontology terms. By
            default, this is `"Blues"`. Palettes can be found at
            matplotlib.org/stable/tutorials/colors/colormaps.html.
            If `None`, the color will be set to the `edge_color`.
        gene_palette : str, optional
            The color palette to use for the differentially expressed
            genes. By default, this is `"Reds"`. Palettes can be found
            at matplotlib.org/stable/tutorials/colors/colormaps.html.
            If `None`, the color will be set to the `edge_color`.
        gene_color : str, optional
            The color of the differentially expressed genes. By default,
            this is `None`. If `None`, the color will be set to the
            `gene_palette`. Otherwise, the color will be set to the
            `gene_color`.
        term_color : str, optional
            The color of the gene ontology terms. By default, this is
            `None`. If `None`, the color will be set to the `term_palette`.
            Otherwise, the color will be set to the `term_color`.
        """
        self._degs = degs
        self._go = go
        self._deg_size_name = deg_size_name
        self._deg_gene_name = deg_gene_name
        self._go_size_name = go_size_name
        self._go_term_name = go_term_name
        self._go_overlap_name = go_overlap_name
        self._neg_log_xform_degs = neg_log_xform_degs
        self._neg_log_xform_go = neg_log_xform_go
        self._set_deg_mass = set_deg_mass
        self._set_go_mass = set_go_mass
        self._edge_color = edge_color
        self._gene_palette = gene_palette
        self._term_palette = term_palette
        self._gene_color = gene_color
        self._term_color = term_color

        self._validate_degs()
        self._validate_go()
        self._build_bipartite_graph()

    def _validate_degs(self):
        for col in [self._deg_size_name, self._deg_gene_name]:
            if col not in self._degs.columns:
                raise ValueError(f"Column '{col}' not found in DEGs.")
        logging.info(f"Found {self._degs.shape[0]} differentially expressed genes.")

    def _validate_go(self):
        for col in [self._go_size_name, self._go_term_name, self._go_overlap_name]:
            if col not in self._go.columns:
                raise ValueError(f"Column '{col}' not found in GO terms.")
        logging.info(f"Found {self._go.shape[0]} gene ontology terms.")

    def _build_deg_attributes(self):
        genes = self._degs[self._deg_gene_name].values
        attributes = self._degs[self._deg_size_name].values
        if self._neg_log_xform_degs:
            attributes = -np.log10(attributes)
        self._gene_attributes = dict(zip(genes, attributes))

    def _build_go_attributes(self):
        terms = self._go[self._go_term_name].values
        overlaps = self._go[self._go_overlap_name].values
        attributes = self._go[self._go_size_name].values

        if self._neg_log_xform_go:
            attributes = -np.log10(attributes)

        self._term_attributes = {
            term: {"attribute": attribute, "overlap": overlap}
            for term, attribute, overlap in zip(terms, attributes, overlaps)
        }

    def _insert_go_term(self, term: str):
        """
        Inserts a GO term into the graph.
        """
        self.graph.add_node(
            term,
            cond="term",
            shape="square",
            attr=self._term_attributes[term]["attribute"],
            size=self._term_attributes[term]["attribute"],
            mass=self._term_attributes[term]["attribute"] if self._set_go_mass else 1.0,
            color=self._term_color
            if self._term_color is not None
            else self._edge_color,
            color_border="black",
        )

    def _insert_deg_gene(self, gene: str):
        """
        Inserts a DEG gene into the graph.
        """
        self.graph.add_node(
            gene,
            cond="gene",
            shape="circle",
            attr=self._gene_attributes[gene],
            size=self._gene_attributes[gene],
            mass=self._gene_attributes[gene] if self._set_deg_mass else 1.0,
        )

    def _insert_edge(self, term: str, gene: str):
        """
        Inserts an edge between a GO term and a DEG gene into the graph.
        """
        self.graph.add_edge(term, gene, color=self._edge_color)

    def _set_color(self, cond: str):
        distribution = []

        for node in self.graph.nodes(data=True):
            if node[1]["cond"] == cond:
                distribution.append(node[1]["attr"])
        distribution = np.array(distribution)

        palette = self._term_palette if cond == "term" else self._gene_palette
        hex_colors = _array_to_hex(distribution, palette)

        idx = 0
        for node in self.graph.nodes(data=True):
            if node[1]["cond"] == cond:
                node[1]["color"] = hex_colors[idx]
                idx += 1

    def _build_bipartite_graph(self):
        self._build_deg_attributes()
        self._build_go_attributes()

        self.graph = nx.Graph()
        for term in self._term_attributes:
            self._insert_go_term(term)
            for gene in self._term_attributes[term]["overlap"]:
                self._insert_deg_gene(gene)
                self._insert_edge(term, gene)

        if self._term_palette is not None and self._term_color is None:
            self._set_color("term")
        if self._gene_palette is not None and self._gene_color is None:
            self._set_color("gene")

        self._build_statistics()

    def _build_statistics(self):
        self._num_terms = len(
            [n for n in self.graph.nodes(data=True) if n[1]["cond"] == "term"]
        )
        self._num_genes = len(
            [n for n in self.graph.nodes(data=True) if n[1]["cond"] == "gene"]
        )
        self._num_edges = len(self.graph.edges)

        logging.info(
            f"Built bipartite graph with {self._num_terms} terms, "
            f"{self._num_genes} genes, and {self._num_edges} edges."
        )

    def visualize(
        self,
        filepath: str = "network.html",
        height: str = "1000px",
        width: str = "100%",
        physics_toggle: bool = True,
        physics: bool = True,
        notebook: bool = False,
        **kwargs,
    ):
        """
        Creates the pyvis visualization. For more information on the
        additional keyword arguments, see `pyvis.network.Network`.

        Parameters
        ----------
        height : str, optional
            The height of the visualization. By default, this is `"1000px"`.
        width : str, optional
            The width of the visualization. By default, this is `"100%"`.
        kwargs
            Additional keyword arguments to pass to `pyvis.network.Network`.
        """
        net = Network(height=height, width=width, notebook=notebook, **kwargs)
        net.from_nx(self.graph)
        net.toggle_physics(physics_toggle)
        if physics_toggle:
            net.show_buttons(filter_=["physics"])
        net.write_html(filepath)
        logging.info(f"Visualization saved to {filepath}.")
