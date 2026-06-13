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

    """
    Model application visualization class.
    """
    def channel_generic(
            self, 
            name: str,
            title: str|None=None,
            zys: Sequence|None=None, 
            do_equal_aspect=False,
            text_labels: Sequence|None=None,
            fig_size: tuple[float,float]=(6, 4,),
        ) -> None:            
        """
        Plot numerical solutions applied to channel cross-section model (vertical profiles).
        """            
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure
        #     zy_list (list): 
        #         set of numerical solutions to plot
        #     text_label (list): text annotation as list of form (x-y coordinate, string, 
        #                     font size)
        #     do_equal_aspect (bool): 
        #         flag whether to use force equal sizing of x and y axis scales
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        zy=zys[0]
        plt.plot(
            zy[2], zy[0], label=zy[4], color="k",
        )
        plt.ylabel(zy[1])
        plt.xlabel(zy[3])
        axes = plt.gca()
        if do_equal_aspect:
            axes.set_aspect("equal")
        else:
            pass
        if text_labels is not None:
            for text_label in text_labels:
                plt.text(
                    *text_label[0], 
                    text_label[1], 
                    color=text_label[3], 
                    size=text_label[2],
                    verticalalignment="center", 
                    horizontalalignment="center",
                    transform=axes.transAxes, 
                    rotation=text_label[4],
                )
        plt.grid("on",ls=":")
        if len(zys)>=2:
            zy=zys[1]
            plt.plot(0,0, label=zy[4], color="forestgreen",)
        plt.legend()
        
        if len(zys)>=2:
            zy=zys[1]
            alt_axes = axes.twiny()
            alt_axes.plot(zy[2], zy[0], label=zy[4], color="forestgreen",)
            alt_axes.set_xlabel(zy[3], color="forestgreen",)
        plt.grid(ls=":")

    def channel_refweatheringrate_referosionrate_W(
            self, 
            name: str,
            title: str|None=None,
            model: Any|None=None, 
            text_label: Sequence|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:   
        """
        Plot numerical solutions applied to channel cross-section model 
        (vertical profiles).
        """            
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        #     cw (:class:`~.solve1p1d.ChannelWall`): instance of :mod:`~.solve1p1d` model 
        #                         class that simulates channel cross-sectional geometry
        #     text_label (list): text annotation as list of form (x-y coordinate, string, 
        #                     font size)

        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        plt.plot(
            model.w0_array/model.physical_parameters[k], 
            model.z_array, 
            label=r"$w_0/k$",
        )
        plt.plot(
            model.v0_array, 
            model.z_array, 
            label=r"$u_0$",
        )
        x_limits = plt.xlim()
        y_limits = plt.ylim()
        plt.plot(
            (-1, -1,), 
            label=r"${\mathcal{W}}=w_0/{u_0} k$", 
            color="forestgreen",
        )
        plt.plot(
            model.vs_array,
            model.z_array, 
            label=r"$v_s$", 
            color="k", 
            lw=2,
        )
        plt.xlim(x_limits)
        plt.ylim(y_limits)
        plt.xlabel(r"Speeds $w_0(z)/k$, ${u_0}(z)$, $v_s(z)$")
        plt.ylabel(r"Height above bed  $z$")
        plt.legend(loc="upper center")
        plt.grid(ls=":")

        axes = plt.gca()
        alt_axes = axes.twiny()
        alt_axes.plot(
            model.W_array,  
            model.z_array, 
            label=r"${\mathcal{W}}$", 
            color="forestgreen",
        )
        alt_axes.set_xlabel(
            r"Weathering number  ${\mathcal{W}}(z)$", 
            color="forestgreen",
        )
        x_limits = axes.get_xlim()
        axes.set_xlim(x_limits[0], x_limits[1]*1.05)

        if text_label is not None:
            plt.text(
                *text_label[0], 
                text_label[1], 
                color="k", 
                size=text_label[2],
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )