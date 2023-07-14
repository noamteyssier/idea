from typing import Optional
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import CenteredNorm, to_hex


def _array_to_hex(
    array: np.ndarray,
    palette: str = "viridis",
    center: Optional[float] = None,
) -> np.ndarray:
    """
    Convert an array of values to a hex color palette.
    """
    cmap = plt.get_cmap(palette)
    if center is None:
        norm = plt.Normalize(array.min(), array.max())
    else:
        norm = CenteredNorm(vcenter=center)
    rgba = cmap(norm(array))
    hex_colors = np.array([to_hex(c) for c in rgba])
    return hex_colors
