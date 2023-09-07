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


# Function to calculate perceived brightness from a hex color
def _calculate_brightness(hex_color):
    r, g, b = tuple(int(hex_color[i : i + 2], 16) for i in (1, 3, 5))
    brightness = (r * 299 + g * 587 + b * 114) / 1000
    return brightness


def _is_dark(hex_color):
    return _calculate_brightness(hex_color) < 128
