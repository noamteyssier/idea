from typing import Optional, Union
import numpy as np
import pandas as pd
import ggetrs as gg


def run_gsea(
    genes: Union[list, np.ndarray],
    library: str = "BP",
    threshold: Optional[float] = 0.05,
    background: Optional[Union[list, np.ndarray]] = None,
    use_pvalue: bool = False,
) -> pd.DataFrame:
    """
    Run gene ontology analysis.

    Parameters
    ----------
    genes : Union[list, np.ndarray]
        A list of genes to run gene ontology analysis on.
    library : str, optional
        The gene ontology library to use. One of `Enrichr.LIBRARIES`.
        Shortcuts for:
            - `BP`: `GO_Biological_Process_2023`
            - `CC`: `GO_Cellular_Component_2023`
            - `MF`: `GO_Molecular_Function_2023`
            - `msigdb`: `MSigDB_Hallmark_2020`
            - `kegg`: `KEGG_2021_Human`
    threshold : Optional[float], optional
        The threshold to use for filtering the gene ontology terms
        by p-value. If `None`, no filtering is done.
    background : Optional[Union[list, np.ndarray]], optional
        The background genes to use for the gene ontology analysis.
    use_pvalue : bool, optional
        Whether to use the p-value or the adjusted p-value for filtering
        the gene ontology terms. (Default: using the adjusted p-value)

    Returns
    -------
    pd.DataFrame
    """
    if library == "BP":
        library = "GO_Biological_Process_2023"
    elif library == "CC":
        library = "GO_Cellular_Component_2023"
    elif library == "MF":
        library = "GO_Molecular_Function_2023"
    elif library == "msigdb":
        library = "MSigDB_Hallmark_2020"
    elif library == "kegg":
        library = "KEGG_2021_Human"

    if background is None:
        response = gg.enrichr(
            library_name=library,
            gene_list=[g for g in genes],
        )
    else:
        response = gg.enrichr_background(
            library_name=library,
            gene_list=[g for g in genes],
            background_list=[b for b in background],
        )
    frame = pd.DataFrame(response[library])

    if threshold is not None:
        if use_pvalue:
            frame = frame[frame["pvalue"] < threshold]
        else:
            frame = frame[frame["adj_pvalue"] < threshold]

    return frame
