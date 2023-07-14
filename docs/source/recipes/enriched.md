# Enriched DEGs

What most people think of when you throw the term differentially expressed genes (DEGs)
is the genes that are enriched.
Specifically these are the genes that are overrepresented in the test case compared to
the controls.

This module was built with these as the default, and the color palettes reflect that.

Let's run through an analysis and visualization only considering the enrichments.

```python
import idea
import pandas as pd

# URL to our example DEG dataframe
url = "https://github.com/noamteyssier/idea/raw/main/example_data/AP2S1.tab.gz"

# Load in our DEG analysis
frame = pd.read_csv(url, sep="\t")

# Filter for significance and enrichments
sig_degs = frame[
    (frame.padj < 0.05) &
    (frame.log2FoldChange > 0)
]

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
    gsea.head(10), # only showing top 10 terms for minimal example
)
id.visualize("enriched_network.html")
```

```{raw} html
:file: ../../../assets/enriched_network.html
```
