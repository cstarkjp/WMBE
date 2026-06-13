"""
Visualization classes.
"""
import warnings
import matplotlib.pyplot as plt

from typing import Any
from collections.abc import Sequence

from wmbe.symbols import *
from wmbe.viz import BaseViz

warnings.filterwarnings("ignore")

__all__ = [
    "AppViz",
]

class AppViz(BaseViz):
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
            
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure
            zy_list (list): 
                set of numerical solutions to plot
            text_label (list): text annotation as list of form (x-y coordinate, string, 
                            font size)
            do_equal_aspect (bool): 
                flag whether to use force equal sizing of x and y axis scales
        """
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
        Plot numerical solutions applied to channel cross-section model (vertical profiles).
            
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
            cw (:class:`~.solve1p1d.ChannelWall`): instance of :mod:`~.solve1p1d` model 
                                class that simulates channel cross-sectional geometry
            text_label (list): text annotation as list of form (x-y coordinate, string, 
                            font size)
        """
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        plt.plot(
            model.w0_array/model.physical_parameters[k], 
            model.z_array, 
            label="$w_0/k$",
        )
        plt.plot(
            model.v0_array, 
            model.z_array, 
            label="${u_0}$",
        )
        x_limits = plt.xlim()
        y_limits = plt.ylim()
        plt.plot(
            (-1, -1,), 
            label="${\mathcal{W}}=w_0/{u_0} k$", 
            color="forestgreen",
        )
        plt.plot(
            model.vs_array,
            model.z_array, 
            label="$v_s$", 
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
            label="${\mathcal{W}}$", 
            color="forestgreen",
        )
        alt_axes.set_xlabel(
            "Weathering number  ${\mathcal{W}}(z)$", color="forestgreen",
        )
        x_limits = axes.get_xlim()
        axes.set_xlim(x_limits[0],x_limits[1]*1.05)

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