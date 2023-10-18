import logging
from typing import Optional, Union
import numpy as np
import pandas as pd
import networkx as nx
from ._utils import _array_to_hex, _is_dark
from ._constants import DIVERGING_PALETTES
from pvsvg import Network


class IDEA:
    def __init__(
        self,
        degs: pd.DataFrame,
        go: pd.DataFrame,
        deg_size_name: str = "padj",
        deg_color_name: str = "padj",
        deg_gene_name: str = "gene",
        go_size_name: str = "adj_pvalue",
        go_term_name: str = "term_name",
        go_overlap_name: str = "overlapping_genes",
        neg_log_xform_degs_size: bool = True,
        neg_log_xform_degs_color: bool = True,
        absolute_degs_color: bool = True,
        neg_log_xform_go: bool = True,
        set_deg_mass: bool = True,
        set_go_mass: bool = True,
        edge_color: str = "grey",
        edge_width: float = 1.0,
        gene_palette: str = "Reds",
        term_palette: str = "Blues",
        gene_color: Optional[str] = None,
        term_color: Optional[str] = None,
        center: Optional[float] = None,
        fontsize: int = 14,
        fontface: str = "arial",
        fontcolor: str = "auto",
        deg_node_scalar: float = 1.0,
        go_node_scalar: float = 1.0,
        force_options: bool = False,
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
        deg_color_name : str, optional
            The name of the column in `degs` that should be used for
            the color of the differentially expressed genes. By default,
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
        neg_log_xform_degs_size : bool, optional
            Whether or not to invert the attributes of the differentially
            expressed genes. By default, this is `True`.
        neg_log_xform_degs_color : bool, optional
            Whether or not to invert the attributes of the differentially
            expressed genes. By default, this is `True`.
        absolute_degs_color: bool, optional
            Whether or not to take the absolute value of the color attribute
            for the differentially expressed genes. By default, this is `True`.
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
        center : float, optional
            The center of the DEG color scale. Set to zero for diverging
            color scales. By default, this is `None`.
        fontsize : int, optional
            The fontsize of the labels. By default, this is `14`.
        edge_width : float, optional
            The width of the edges. By default, this is `1.0`.
        fontface : str, optional
            The fontface of the labels. By default, this is `"arial"`.
        fontcolor : str, optional
            The fontcolor of the labels. By default, this is `"black"`.
        deg_node_scalar : float, optional
            The scalar to multiply the size of the DEG nodes by. By default,
            this is `1.0`.
        go_node_scalar : float, optional
            The scalar to multiply the size of the GO nodes by. By default,
            this is `1.0`.
        force_options : bool, optional
            Whether or not to force the options provided even if they don't match
            expectations. By default, this is `False`.
        """
        self._degs = degs
        self._go = go
        self._deg_size_name = deg_size_name
        self._deg_color_name = deg_color_name
        self._deg_gene_name = deg_gene_name
        self._go_size_name = go_size_name
        self._go_term_name = go_term_name
        self._go_overlap_name = go_overlap_name
        self._neg_log_xform_degs_size = neg_log_xform_degs_size
        self._neg_log_xform_degs_color = neg_log_xform_degs_color
        self._absolute_degs_color = absolute_degs_color
        self._neg_log_xform_go = neg_log_xform_go
        self._set_deg_mass = set_deg_mass
        self._set_go_mass = set_go_mass
        self._edge_color = edge_color
        self._edge_width = edge_width
        self._gene_palette = gene_palette
        self._term_palette = term_palette
        self._gene_color = gene_color
        self._term_color = term_color
        self._center = center
        self._fontsize = fontsize
        self._fontface = fontface
        self._fontcolor = fontcolor
        self._deg_node_scalar = deg_node_scalar
        self._go_node_scalar = go_node_scalar
        self._force_options = force_options
        if self._fontcolor == "auto":
            self._font_shorthand = f"{self._fontsize}px {self._fontface} black"
        else:
            self._font_shorthand = (
                f"{self._fontsize}px {self._fontface} {self._fontcolor}"
            )

        self._validate_options()
        self._validate_degs()
        self._validate_go()
        self._warn_options()
        self._build_bipartite_graph()

    def _validate_degs(self):
        for col in [self._deg_size_name, self._deg_color_name, self._deg_gene_name]:
            if col not in self._degs.columns:
                raise ValueError(f"Column '{col}' not found in DEGs.")
        logging.info(f"Found {self._degs.shape[0]} differentially expressed genes.")

    def _validate_go(self):
        for col in [self._go_size_name, self._go_term_name, self._go_overlap_name]:
            if col not in self._go.columns:
                raise ValueError(f"Column '{col}' not found in GO terms.")
        logging.info(f"Found {self._go.shape[0]} gene ontology terms.")

    def _validate_options(self):
        if self._deg_node_scalar <= 0:
            raise ValueError("deg_node_scalar must be greater than zero.")
        if self._go_node_scalar <= 0:
            raise ValueError("go_node_scalar must be greater than zero.")

    def _warn_options(self):
        if (
            not np.all(self._degs[self._deg_size_name] > 0)
            and self._neg_log_xform_degs_size
        ):
            logging.warning("Not all DEG sizes are positive. ")
            if not self._force_options:
                logging.warning("Setting neg_log_xform_degs_size to False.")
                self._neg_log_xform_degs_size = False
            else:
                logging.warning("Continuing with neg_log_xform_degs_size set to True.")

        if (
            not np.all(self._degs[self._deg_color_name] > 0)
            and self._neg_log_xform_degs_color
        ):
            logging.warning("Not all DEG colors are positive. ")
            if not self._force_options:
                logging.warning("Setting neg_log_xform_degs_color to False.")
                self._neg_log_xform_degs_color = False
            else:
                logging.warning("Continuing with neg_log_xform_degs_color set to True.")

        if (
            not np.all(self._degs[self._deg_color_name] > 0)
            and self._absolute_degs_color
        ):
            logging.warning("Not all DEG colors are positive")
            if not self._force_options:
                logging.warning("Setting absolute_degs_color to False.")
                self._absolute_degs_color = False

                if self._gene_palette not in DIVERGING_PALETTES:
                    logging.warning("gene_palette is not a diverging palette.")
                    logging.warning(
                        "Setting palette to 'seismic'. for diverging color scale."
                    )
                    self._gene_palette = "seismic"
            else:
                logging.warning("Continuing with absolute_degs_color set to True.")

        if not np.all(self._go[self._go_size_name] > 0) and self._neg_log_xform_go:
            logging.warning("Not all GO sizes are positive. ")
            if not self._force_options:
                logging.warning("Setting neg_log_xform_go to False.")
                self._neg_log_xform_go = False
            else:
                logging.warning("Continuing with neg_log_xform_go set to True.")

        if self._gene_color is not None and self._gene_palette is not None:
            logging.warning("Both gene_color and gene_palette are set. ")
            if not self._force_options:
                logging.warning("Setting gene_color to None.")
                self._gene_color = None
            else:
                logging.warning("Continuing with gene_color set to not None.")

    def _build_deg_attributes(self):
        genes = self._degs[self._deg_gene_name].values.astype(str)
        sizes = self._degs[self._deg_size_name].values.astype(float)
        colors = self._degs[self._deg_color_name].values.astype(float)
        if self._neg_log_xform_degs_size:
            sizes = -np.log10(sizes)
        if self._absolute_degs_color:
            colors = np.abs(colors)
        if self._neg_log_xform_degs_color:
            colors = -np.log10(colors)
        self._gene_attributes = {
            gene: {"attr_size": size, "attr_color": color}
            for gene, size, color in zip(genes, sizes, colors)
        }

    def _build_go_attributes(self):
        terms = self._go[self._go_term_name].values.astype(str)
        overlaps = self._go[self._go_overlap_name].values
        attributes = self._go[self._go_size_name].values.astype(float)

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
            title=term,
            label=term,
            attr=self._term_attributes[term]["attribute"],
            size=self._term_attributes[term]["attribute"] * self._go_node_scalar,
            mass=self._term_attributes[term]["attribute"] if self._set_go_mass else 1.0,
            color=self._term_color
            if self._term_color is not None
            else self._edge_color,
            color_border="black",
            font=self._font_shorthand,
        )

    def _insert_deg_gene(self, gene: str):
        """
        Inserts a DEG gene into the graph.
        """
        self.graph.add_node(
            gene,
            cond="gene",
            shape="ellipse",
            title=gene,
            label=gene,
            color_border="black",
            attr=self._gene_attributes[gene]["attr_color"],
            size=self._gene_attributes[gene]["attr_size"] * self._deg_node_scalar,
            mass=self._gene_attributes[gene]["attr_size"]
            if self._set_deg_mass
            else 1.0,
            font={
                "face": self._fontface,
                "color": self._fontcolor if self._fontcolor != "auto" else "black",
            },
        )

    def _insert_edge(self, term: str, gene: str):
        """
        Inserts an edge between a GO term and a DEG gene into the graph.
        """
        self.graph.add_edge(
            term,
            gene,
            color=self._edge_color,
            weight=self._edge_width,
        )

    def _valid_gene(self, gene: str) -> bool:
        """
        Validates that a gene is in the DEGs.
        """
        return gene in self._gene_attributes

    def _set_color(self, cond: str):
        distribution = []

        for node in self.graph.nodes(data=True):
            if node[1]["cond"] == cond:
                distribution.append(node[1]["attr"])
        distribution = np.array(distribution)

        palette = self._term_palette if cond == "term" else self._gene_palette

        if cond == "term":
            hex_colors = _array_to_hex(distribution, palette)
        else:
            hex_colors = _array_to_hex(distribution, palette, self._center)

        idx = 0
        for node in self.graph.nodes(data=True):
            if node[1]["cond"] == cond:
                node_color = hex_colors[idx]
                node[1]["color"] = node_color
                if cond == "gene":
                    if self._fontcolor == "auto":
                        node[1]["font"]["color"] = (
                            "white" if _is_dark(node_color) else "black"
                        )
                idx += 1

    def _build_bipartite_graph(self):
        self._build_deg_attributes()
        self._build_go_attributes()

        self.graph = nx.Graph()
        for term in self._term_attributes:
            self._insert_go_term(term)
            for gene in self._term_attributes[term]["overlap"]:
                if not self._valid_gene(gene):
                    continue
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
        height: Union[int, str] = "800px",
        width: Union[int, str] = "100%",
        **kwargs,
    ):
        """
        Creates the vis.js visualization. For more information on the
        additional keyword arguments, see `https://github.com/noamteyssier/pvsvg`.

        Parameters
        ----------
        height : str, optional
            The height of the visualization. By default, this is `"1000px"`.
        width : str, optional
            The width of the visualization. By default, this is `"100%"`.
        kwargs
            Additional keyword arguments to pass to `pvsvg.Network`.
        """
        net = Network(
            graph=self.graph,
            height=height,
            width=width,
        )
        net.draw(filepath)
        logging.info(f"Visualization saved to {filepath}.")
