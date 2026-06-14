"""
Visualization classes.
"""
import warnings
import logging
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from typing import Any
from collections.abc import Sequence

from wmbe.symbols import *

warnings.filterwarnings("ignore")

__all__ = [
    "VizBase",
]

class VizBase:
    """
    Provide a visualization class.
    """
    # Args:
    #     dpi:
    #         set resolution for rasterized images
    #     font_size:
    #         set mpl default font size
    #     font_family:
    #         set mpl default font family

    # Attributes:
    #     dpi (int):
    #         rasterization resolution
    #     fdict  (dict):
    #         dictionary to which each figure is appended as it is generated
    dpi: int
    fdict: dict[Any, Any]

    def __init__(
            self, 
            dpi: int=150, 
            font_size: int=11, 
            font_family: str="Arial",
        ) -> None:
        self.dpi = dpi
        self.fdict = {}
        try:
            mpl.rc("font", size=font_size, family=font_family)
        except:
            mpl.rc("font", size=font_size, family="")

        self.markers = {
            'o': 'circle', 
            'v': 'triangle_down', 
            '^': 'triangle_up', 
            '<': 'triangle_left', 
            '>': 'triangle_right', 
            #  '1': 'tri_down', '2': 'tri_up', 
            #  '3': 'tri_left', '4': 'tri_right',
            #  's': 'square', 
            '8': 'octagon', 
            'p': 'pentagon', 
            # '*': 'star', 
            'h': 'hexagon1',
            'H': 'hexagon2', 
            #  '+': 'plus', 'x': 'x', 
            'D': 'diamond', 
            'd': 'thin_diamond', 
            '|': 'vline', 
            '_': 'hline', 
            'P': 'plus_filled', 
            'X': 'x_filled', 
            #  0: 'tickleft', 1: 'tickright', 
            #  2: 'tickup', 3: 'tickdown', 
            #  4: 'caretleft', 5: 'caretright', 6: 'caretup',
            #  7: 'caretdown', 8: 'caretleftbase', 
            #  9: 'caretrightbase', 10: 'caretupbase', 11: 'caretdownbase'
        }

    def create_figure(
        self,
        name: str,
        size: tuple[float, float] | None = None,
        dpi: int | None = None,
    ) -> Figure:
        """
        Initialize a Pyplot figure.

        Set its size and DPI. Append it to the figures dictionary.
        """
        # Args:
        #     fig_name:
        #         name of figure; used as key in figures dictionary
        #     fig_size:
        #         optional width and height of figure in inches
        #     dpi:
        #         rasterization resolution

        # Returns:
        #     figure:
        #         reference to MatPlotLib/Pyplot figure

        fig_size_: tuple[float, float] = (
            (6, 4,) if size is None else size
        )
        dpi_: float = self.dpi if dpi is None else dpi
        logging.info(
            "Viz:\n   "
            + f"Creating plot: {name} size={fig_size_} @ {dpi_} dpi"
        )
        fig = plt.figure()
        self.fdict.update({name: fig})
        if fig_size_ is not None:
            fig.set_size_inches(*fig_size_)
        fig.set_dpi(dpi_)
        return fig
