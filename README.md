# IDEA

## Integrated Differential Expression and Annotation

<p align="center">
    <img src="https://github.com/noamteyssier/idea/blob/main/assets/logo.png?raw=true" width="50%" />
</p>

This is a python module to perform GO analysis using [Enrichr](https://maayanlab.cloud/Enrichr/)
and visualize the [bipartite graph](https://en.wikipedia.org/wiki/Bipartite_graph)
of terms and genes as an interactive force-directed graph.

This uses [`pyvis`](https://pyvis.readthedocs.io/en/latest/tutorial.html) as the
force-directed graph backend and [`ggetrs`](https://noamteyssier.github.io/ggetrs)
to perform the gene set enrichment using the [`Enrichr`](https://maayanlab.cloud/Enrichr/) API.

## Installation and Usage

Checkout my [documentation](https://idea-bio.readthedocs.io/en/latest/index.html)
for information on installation, usage, and recipes for analysis.
