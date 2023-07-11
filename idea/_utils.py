import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex


def _array_to_hex(
    array: np.ndarray,
    palette: str = "viridis",
) -> np.ndarray:
    """
    Convert an array of values to a hex color palette.
    """
    cmap = plt.get_cmap(palette)
    norm = plt.Normalize(array.min(), array.max())
    rgba = cmap(norm(array))
    hex_colors = np.array([to_hex(c) for c in rgba])
    return hex_colors
