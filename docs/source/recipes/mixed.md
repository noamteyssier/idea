## Mixed DEGs

Sometimes more data is more better - let's run an analysis where we do a GSEA using **both**
enriched and depleted DEGs.

We'll need to configure the visualization to color based on the `log2FoldChange` so we can
differentiate the two groups more easily in the network.
We'll also change the term palette to gradiate through shades of green to avoid color collisions
between the two classes.

We'll also turn off the default DEG color transformations (which are expecting p-values
by default), and set our color-scale center to 0.

```python
import idea
import pandas as pd

# URL to our example DEG dataframe
url = "https://github.com/noamteyssier/idea/raw/main/example_data/AP2S1.tab.gz"

# Load in our DEG analysis
frame = pd.read_csv(url, sep="\t")

# Filter for significance and depletions
sig_degs = frame[ (frame.padj < 0.05) ]

# Select the gene names
geneset = sig_degs.gene.values

# Run the GSEA
gsea = idea.run_go(
    geneset,
    threshold=0.05,
    library="BP",
)

# Build and Visualize Network
id = idea.IDEA(
    sig_degs,
    gsea.head(15), # only showing top 15 terms to show minimal example
    gene_palette="bwr", # gradiate through blue-white-red
    term_palette="Greens", # gradiate through shades of green
    deg_color_name="log2FoldChange", # select the column name for FC in the DEGs
    neg_log_xform_degs_color=False, # disable negative log transformation of deg color attribute
    absolute_degs_color=False, # disable absolute value transformation of deg color attribute
    center=0, # set the center of our diverging color scale to 0
)
id.visualize("mixed_network.html")
```

```{raw} html
:file: ../../../assets/mixed_network.html
```

