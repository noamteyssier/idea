# Usage

You can import the module using its shorthand name: `idea`

```
import idea
```

## Gene Set Enrichment

### Introduction to Gene Set Enrichment Analysis

[Gene Set Enrichment Analysis (GSEA)](https://en.wikipedia.org/wiki/Gene_set_enrichment_analysis)
is a a method to identify classes of genes or proteins that are **overrepresented** in a large set
of genes and proteins and may have an association with specific phenotypes.

The method generally involves statistically testing whether a set of genes/proteins is overrepresented
in a collection of terms, each of which has an associated set of genes/proteins, and returns the set
of terms with which the original gene set is significantly overrepresented.

A common tool people use for GSEA is the [`Enrichr`](https://maayanlab.cloud/Enrichr) cloud,
which is an online tool which lets you drop in a gene set and perform the GSEA quickly and easily.

### Performing GSEA

This module allows you to perform a GSEA directly from your python session.
The GSEA is performed using [`Enrichr`](https://maayanlab.cloud/Enrichr/)
whose API is accessed using [`ggetrs`](https://noamteyssier.github.io/ggetrs).

The function for running the GSEA is `idea.run_gsea`:

#### GSEA Usage

```python
gsea = idea.run_gsea(
    geneset,
    threshold=0.05,
    library="BP"
)
```

In this example:

- `geneset` is assumed to be a list, or list-like, of gene symbols.
- `threshold` is used to filter the resulting GSEA dataframe on the adjusted p-value
- `library` is used to select the [`Enrichr library`](https://maayanlab.cloud/Enrichr/#libraries)
  to run the GSEA with.
    - this option has some shorthands built in to avoid typing out long names for common libraries.

The result of this function is a pandas `DataFrame` which can be manipulated for whatever downstream
purposes you wish.

## Network Visualization

The core of this library is the generation of the [`Bipartite Graph`](https://en.wikipedia.org/wiki/Bipartite_graph)
of:

1. gene set **terms** 
2. differentially expressed **genes**.

> **Note:**
>
> All edges in this graph will **only** be between the terms and genes and not within each group.

Each term in the GSEA has a set of overlapping genes from the geneset, and the graph is created
by adding an undirected edge between all terms and their overlapping genes.

The resulting networks are informative to seeing how your gene set is distributed across the
enriched terms, and helps pick out the genes and gene clusters that are shared members of terms.

### Building the Network

The `idea.IDEA` class performs the network construction and also handles attribute selection
for visualization purposes.

We can build our network by instantiating the class:

```python
id = idea.IDEA(degs, gsea)
```

That's it! 

### Visualizing the Network

We can then create an interactive HTML visualization (powered by [pyvis](https://pyvis.readthedocs.io/en/latest/#)) 
of our network using the `visualize` method of the `IDEA` class:

```python
id.visualize("network.html")
```

This will create a file `network.html` in the current working directory of the python
session. 
You can open this file up in your favorite browser to interact with the network that is
created.

By default the nodes are sized and colored by their significance and shaped by their
term/gene class.
See further instructions on configuration or check out the `IDEA` API for its configuration
options.

```{raw} html
:file: ../assets/example_network.html
```

