# Depleted DEGs

The converse of the previous recipe is the genes that are **depleted** in the test
case compared to the control case.

We can perform the exact same analysis, but let's change the color palette so that
it's less misleading (red means up, blue means down?)

Let's run through an analysis and visualization only considering the enrichments.

```python
import idea
import pandas as pd

# URL to our example DEG dataframe
url = "https://github.com/noamteyssier/idea/raw/main/example_data/AP2S1.tab.gz"

# Load in our DEG analysi#s
frame = pd.read_csv(url, sep="\t")

# Filter for significance and depletions
sig_degs = frame[
    (frame.padj < 0.05) &
    (frame.log2FoldChange < 0)
]

# Select the gene names
geneset = sig_degs.gene.values

# Run the GSEA
gsea = idea.run_gsea(
    geneset,
    threshold=0.05,
    library="BP",
)

# Build and Visualize Network
id = idea.IDEA(
    sig_degs,
    gsea.head(10), # only showing top 10 terms for minimal example
    gene_palette="Blues", # Updating the gene palette to gradiate through shades of blue
    term_palette="Reds", # Updating the term palette to gradiate through shades of red
)
id.visualize("depleted_network.html")
```

```{raw} html
:file: ../../assets/depleted_network.html
```
