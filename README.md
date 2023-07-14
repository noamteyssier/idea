# IDEA

![logo](https://github.com/noamteyssier/idea/blob/main/assets/logo.png?raw=true)

**I**ntegrated **D**ifferential **E**xpression and **A**nnotation

## Introduction

This is a python module to perform GO analysis using [Enrichr](https://maayanlab.cloud/Enrichr/)
and visualize the [bipartite graph](https://en.wikipedia.org/wiki/Bipartite_graph)
of terms and genes as an interactive force-directed graph.

This uses [`pyvis`](https://pyvis.readthedocs.io/en/latest/tutorial.html) as the
force-directed graph backend and [`ggetrs`](https://noamteyssier.github.io/ggetrs)
to perform the gene set enrichment using the [`Enrichr`](https://maayanlab.cloud/Enrichr/) API.
