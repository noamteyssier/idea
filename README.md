# IDEA

![logo](https://github.com/noamteyssier/idea/blob/main/assets/logo.png?raw=true)

**I**ntegrated **D**ifferential **E**xpression and **A**nnotation

## Introduction

This is a python module to perform GO analysis using [Enrichr](https://maayanlab.cloud/Enrichr/)
and visualize the [bipartite graph](https://en.wikipedia.org/wiki/Bipartite_graph)
of terms and genes as an interactive force-directed graph.

This uses [`pyvis`](https://pyvis.readthedocs.io/en/latest/tutorial.html) as the
force-directed graph backend and [`ggetrs`](https://noamteyssier.github.io/ggetrs)
to perform the gene set enrichment using `Enrichr`'s API.

## Installation

You can install this like other python packages using `pip`:

``` bash
pip install idea-bio
```

## Usage

The basic workflow for this tool is made up of 3 steps:

1. Performing the gene set enrichment analysis
2. Constructing the network
3. Visualizing the network

``` python

import idea
import pandas as pd

#################
# Preprocessing #
#################

# Load in our example dataframe
url = "https://github.com/noamteyssier/idea/raw/main/example_data/AP2S1.tab.gz"
deg_frame = pd.read_csv(url, sep="\t")

# Filter to significant enrichments
sig_degs = deg_frame[
    (deg_frame.log2FoldChange > 0) &
    (deg_frame.padj < 0.05)
]

# Select the gene names
geneset = sig_degs.gene.values

########################################
# Perform Gene Set Enrichment Analysis #
########################################

gsea = idea.run_go(
    geneset,
    threshold=0.05,
    library="BP",
)

###############################
# Build and Visualize Network #
###############################

# Build Network
id = idea.IDEA(
    sig_degs,
    gsea.head(30),
)

# Write HTML of network to `network.html`
id.visualize("network.html")
```

### Visualization

We can then visualize and interact with our network by opening
the created `*.html` in our favorite browser.

Here is a static image of an example network:

![network.png](https://github.com/noamteyssier/idea/blob/main/assets/example_network.png?raw=true)
