"""
Data visualization.
"""
import warnings
import logging
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import matplotlib.ticker as ticker

from typing import Any, Callable
from collections.abc import Sequence
from numpy.typing import NDArray

from wmbe.theory import Theory
from wmbe.solve1d import eta_chi_tau, NumericalModel
from wmbe.data import linear_model, ExperimentalData
from wmbe.symbols import *

warnings.filterwarnings("ignore")

__all__ = ["Viz"]

class Viz:
    """
    Provide a visualization class.

    Args:
        dpi:
            set resolution for rasterized images
        font_size:
            set mpl default font size
        font_family:
            set mpl default font family

    Attributes:
        dpi (int):
            rasterization resolution
        fdict  (dict):
            dictionary to which each figure is appended as it is generated
    """

    dpi: int
    fdict: dict[Any, Any]

    def __init__(
            self, dpi: int=150, font_size: int=11, font_family: str="Arial",
        ) -> None:
        self.dpi = dpi
        self.fdict = {}
        try:
            mpl.rc("font", size=font_size, family=font_family)
        except:
            mpl.rc("font", size=font_size, family="")

        self.markers = {
            'o': 'circle', 'v': 'triangle_down', '^': 'triangle_up', '<': 'triangle_left', 
            '>': 'triangle_right', 
            #  '1': 'tri_down', '2': 'tri_up', '3': 'tri_left', '4': 'tri_right',
            #  's': 'square', 
            '8': 'octagon', 'p': 'pentagon', 
            # '*': 'star', 
            'h': 'hexagon1',
            'H': 'hexagon2', 
            #  '+': 'plus', 'x': 'x', 
            'D': 'diamond', 
            'd': 'thin_diamond', '|': 'vline', '_': 'hline', 
            'P': 'plus_filled', 'X': 'x_filled', 
            #  0: 'tickleft', 1: 'tickright', 
            #  2: 'tickup', 3: 'tickdown', 
            #  4: 'caretleft', 5: 'caretright', 6: 'caretup',
            #    7: 'caretdown', 8: 'caretleftbase', 9: 'caretrightbase', 10: 'caretupbase', 11: 'caretdownbase'
        }

    def create_figure(
        self,
        fig_name: str,
        fig_size: tuple[float, float] | None = None,
        dpi: int | None = None,
    ) -> Figure:
        """
        Initialize a Pyplot figure.

        Set its size and DPI. Append it to the figures dictionary.


        Args:
            fig_name:
                name of figure; used as key in figures dictionary
            fig_size:
                optional width and height of figure in inches
            dpi:
                rasterization resolution

        Returns:
            figure:
                reference to MatPlotLib/Pyplot figure
        """
        fig_size_: tuple[float, float] = (
            (6, 4) if fig_size is None else fig_size
        )
        dpi_: float = self.dpi if dpi is None else dpi
        logging.info(
            "Viz:\n   "
            + f"Creating plot: {fig_name} size={fig_size_} @ {dpi_} dpi"
        )
        fig = plt.figure()
        self.fdict.update({fig_name: fig})
        if fig_size_ is not None:
            fig.set_size_inches(*fig_size_)
        fig.set_dpi(dpi_)
        return fig


    def inoue_w_wetdryN(
            self, 
            name: str,
            title: str|None=None,
            expt_data: ExperimentalData|None=None, 
            text_label: str|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        r"""
        Plot Inoue et al data on weakness :math:`w` versus 
        proxy time :math:`N` (number of wet/dry cycles).
        
        Generate graph of rock weakness :math:`w` versus proxy time inferred from 
        `Inoue et al (2017)`_
        data on tensile strength :math:`\sigma_T` 
        (normalized by a reference tensile strength :math:`\sigma_\mathrm{ref}`)
        after :math:`N` wetting and drying cycles.
            
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
            ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
                                        containing data sets as :mod:`pandas` dataframes
            text_label (list): text annotation as list of form (x-y coordinate, string, 
                            font size)
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is not None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        df = expt_data.ddict.get("inoue")
        # sigmaT  = df.sigmaT
        wetdryN = df.wetdryN
        erodibility_sigma2   = df.w_sigma2
        # erodibility_sigma1p5 = df.w_sigma1p5
        erodibility_sigma2_fit = expt_data.fdict["inoue"][0]
            
        plt.plot(wetdryN,linear_model(wetdryN,*erodibility_sigma2_fit),
                color="k",
                label=r"$w \sim 1/\sigma_T^2$")
        plt.errorbar(np.unique(wetdryN),expt_data.w_s2_means, label="",
                    xerr=None,
                    yerr=expt_data.w_s2_stds*1,
                    ecolor="k", mec="w", 
                    color="w", fillstyle="full", 
                    alpha=1,
                    fmt="o", 
                    markersize=0, markeredgewidth=2,
                    elinewidth=1.5,capthick=3,capsize=7)
        plt.plot(
            np.unique(wetdryN),
            expt_data.w_s2_means, 
            label="mean data",
            color="lightgray", 
            fillstyle="full", 
            ls="", 
            marker="o",
            markeredgecolor="k",
            ms=13,
        ) 
        plt.plot(
            wetdryN,
            erodibility_sigma2, 
            label="raw data",  
            color="orange", alpha=0.7,
            ls="", 
            marker="s",
            markeredgecolor="k",
            ms=5,
        )

        plt.legend(loc="upper left")
    #     plt.ylim(0,)
        plt.xlabel(r"Proxy time (no. wet/dry cycles)  $N$  [-]")
        plt.ylabel(r"Weakness  $w=(\sigma_T\left[N\right]/\sigma_{\mathrm{ref}})^{-2}$  [-]")
        # plt.ylabel(r"Weakness  [-]")

        if text_label is not None:
            axes = plt.gca()
            plt.text(*text_label[0], text_label[1], 
                    color="k", size=text_label[2],
                    verticalalignment="center", horizontalalignment="center",
                    transform=axes.transAxes)
        plt.grid(ls=":")
            
    def li_w_wetdryN(
            self, 
            name: str,
            title: str|None=None,
            expt_data: ExperimentalData|None=None, 
            text_label: str|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        r"""
        Plot Li et al  data on weakness :math:`w` versus proxy time :math:`N` 
        (number of wet/dry cycles).
        
        Generate graph of rock weakness :math:`w` versus proxy time :math:`N` 
        inferred from  `Li et al (2016)`_
        data on compressive strength :math:`\sigma_C` 
        (normalized by a reference compressive strength :math:`\sigma_\mathrm{ref}`)
        after :math:`N` wetting and drying cycles
        at a range of confining pressures :math:`P`.
        
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
            ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
                                        containing data sets as :mod:`pandas` dataframes
            text_label (list): text annotation as list of form (x-y coordinate, string, 
                            font size)
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is not None:
            plt.title(title, fontdict={"fontsize": 11.5})

        df = expt_data.ddict["li"]
        color_list = ("darkblue","darkmagenta","firebrick","red","orange","green")
        marker_list = ("o","s","v","^","8","+")
        n_cols = len(color_list)
        for idx,P__ in enumerate(np.unique(df.P)):
            selection_name = "{0}_{1}_{2}".format("li","P",P__)
            # sigma  = df.sigmaC[df.P==P__]/100
            wetdryN = df.wetdryN[df.P==P__]
            erodibility_sigma   = df.w_sigma2[df.P==P__]
            erodibility_sigma_fit = expt_data.fdict[selection_name][0]
            
            plt.plot(wetdryN,linear_model(wetdryN,*erodibility_sigma_fit),
                    color=color_list[idx%n_cols],label="")
            plt.errorbar(wetdryN,erodibility_sigma, 
                        label=r"$P = ${}$\,$MPa".format(P__),
                        xerr=None,
                        yerr=None,
                        ecolor="k", mec="k", 
                        color=color_list[idx%n_cols], fillstyle="full", 
                        alpha=0.7,
                        fmt=marker_list[idx%n_cols], 
                        markersize=7, markeredgewidth=0.5,
                        elinewidth=1.5,capthick=3,capsize=7)    
        
        plt.legend(loc="upper left")
        plt.ylim(0,)
        plt.xlabel("Proxy time (no. wet/dry cycles)  $N$  [-]")
        plt.ylabel(r"Weakness  $w=(\sigma_C\left[N\right]/\sigma_\mathrm{ref})^{-2}$  [-]")
        # plt.ylabel(r"Weakness  $w=(\sigma_C\left[N\right])$  [-]")
        # plt.ylabel(r"Weakness   [-]")

        if text_label is not None:
            axes = plt.gca()
            plt.text(*text_label[0], text_label[1], 
                    color="k", size=text_label[2],
                    verticalalignment="center", horizontalalignment="center",
                    transform=axes.transAxes)
        plt.grid(ls=":")
            
    def li_w_P(
            self, 
            name: str,
            title: str|None=None,
            expt_data: ExperimentalData|None=None, 
            text_label: str|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        r"""
        Plot Li et al  data on weakness :math:`w` versus proxy depth.
        
        Generate graph of rock weakness :math:`w` versus proxy depth 
        inferred from  `Li et al (2016)`_ data on compressive strength :math:`\sigma_C` 
        at a range of confining pressures :math:`P` after :math:`N` 
        wetting and drying cycles.
        
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
            ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
                                        containing data sets as :mod:`pandas` dataframes
            text_label (list): text annotation as list of form (x-y coordinate, string, 
                            font size)
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is not None:
            plt.title(title, fontdict={"fontsize": 11.5})

        color_list = ("darkblue","darkmagenta","firebrick","red","orange","green")
        marker_list = ("o","s","v","^","P","p")
        n_cols = len(color_list)
        
        df = expt_data.ddict["li"]
        # fit = expt_data.fdict["li"]
        sampled_fit = expt_data.sdict["li"]
        # w_ref_vec = np.flipud(sampled_fit[2].T)[:,0]-1
        
        for idx,wetdryN in enumerate(np.flip(np.unique(df.wetdryN))):
            wdN__ = df.wetdryN[df.wetdryN==wetdryN]
            P__ = df.P[df.wetdryN==wetdryN]
            w__ = df.w_sigma2[df.wetdryN==wetdryN]
            P_fit = sampled_fit[1]
            w_fit = np.flipud(sampled_fit[2].T)[idx]
            plt.plot(P_fit, w_fit, color=color_list[idx%n_cols])
            plt.errorbar(P__, w__,
                        label="N = {}".format(wetdryN),
                        xerr=None,
                        yerr=None,
                        ecolor="k", mec="k", 
                        color=color_list[idx%n_cols], fillstyle="full", 
                        alpha=0.7,
                        fmt=marker_list[idx%n_cols], 
                        markersize=7, markeredgewidth=0.5,
                        elinewidth=1.5,capthick=3,capsize=7)   
            
        plt.legend(loc="upper right")
        plt.ylim(0.,)
        plt.autoscale(enable=True, tight=True, axis="x")
        x_limits = plt.xlim()
        plt.plot(x_limits,(1,1), color="gray", ls=":")
        plt.xlabel("Proxy depth (confining pressure $P$)  [MPa]")
        plt.ylabel(r"Weakness  $w=(\sigma_C[N]/\sigma_\mathrm{ref})^{-2}$  [-]")

        if text_label is not None:
            axes = plt.gca()
            plt.text(*text_label[0], text_label[1], 
                    color="k", size=text_label[2],
                    verticalalignment="center", horizontalalignment="center",
                    transform=axes.transAxes)
        plt.grid(ls=":")
                
    def li_w_normed_P(
            self, 
            name: str,
            title: str|None=None,
            expt_data: ExperimentalData|None=None, 
            text_label: str|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        r"""
        Plot Li et al data on weakness :math:`w` (normalized using 2D model) 
        versus proxy depth.
        
        Generate graph of rock weakness versus proxy depth
        inferred from  `Li et al (2016)`_
        data on compressive strength :math:`\sigma_C` 
        at a range of confining pressures :math:`P` after :math:`N` wetting and drying cycles.
        
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
            ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
                                        containing data sets as :mod:`pandas` dataframes
            text_label (list): text annotation as list of form (x-y coordinate, string, 
                            font size)
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is not None:
            plt.title(title, fontdict={"fontsize": 11.5})

        color_list = ("darkblue","darkmagenta","firebrick","red","orange","green")
        marker_list = ("o","s","v","^","P","p")
        n_cols = len(color_list)
        
        df = expt_data.ddict["li"]
        # fit = expt_data.fdict["li"]
        sampled_fit = expt_data.sdict["li"]
        w_ref_vec = np.flipud(sampled_fit[2].T)[:,0]-1
        
        P_fit = sampled_fit[1]
        w_fit = (np.flipud(sampled_fit[2].T)[-1]-1)/w_ref_vec[-1]+1

        plt.errorbar(
            np.unique(df.P),
            expt_data.w_s2normed_means, 
            label="mean data",
            xerr=None,
            yerr=expt_data.w_s2normed_stds*2,
            ecolor="k", mec="k", 
            color="lightgray", 
            fillstyle="full", 
            alpha=0.7,
            fmt="o", 
            markersize=13, 
            markeredgewidth=1.5,
            elinewidth=1.5,
            capthick=3,
            capsize=7,
        )
        
        for idx,wetdryN in enumerate(np.flip(np.unique(df.wetdryN))):
            # wdN__ = df.wetdryN[df.wetdryN==wetdryN]
            P__ = df.P[df.wetdryN==wetdryN]
            w_normed = df.w_s2normed[df.wetdryN==wetdryN]
            plt.errorbar(
                P__, 
                w_normed,
                label="N = {}".format(wetdryN),
                xerr=None,
                yerr=None,
                ecolor="k", 
                mec="k", 
                color=color_list[idx%n_cols], 
                fillstyle="full", 
                alpha=0.7,
                fmt=marker_list[idx%n_cols], 
                markersize=7, 
                markeredgewidth=0.5,
                elinewidth=1.5,
                capthick=3,
                capsize=7,
            )   
                    
        plt.plot(P_fit, w_fit, color="darkblue", label=r"$w \sim \exp(-k\chi)$")

        plt.legend(loc="upper right")
        plt.ylim(0.8,)
        plt.autoscale(enable=True, tight=True, axis="x")
        x_limits = plt.xlim()
        plt.plot(x_limits,(1,1), color="gray", ls=":")
        plt.xlabel("Proxy depth (confining pressure $P$)  [MPa]")
        plt.ylabel("Normalized weakness  $w(\\tau,\\chi)/w(0,\\tau)$  [-]")
        
        if text_label is not None:
            axes = plt.gca()
            plt.text(*text_label[0], text_label[1], 
                    color="k", size=text_label[2],
                    verticalalignment="center", horizontalalignment="center",
                    transform=axes.transAxes)
        plt.grid(ls=":")
            
    def li_w_wetdryN_P(
            self, 
            name: str,
            title: str|None=None,
            expt_data: ExperimentalData|None=None, 
            model_surface: str|None=None,
            text_label: str|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        r"""
        Plot Li et al data in 3D on weakness versus proxy depth and confining pressure.
        
        Generate 3D view of 2D surface model & data
        of rock weakness versus proxy depth inferred from 
        `Li et al (2016)`_ data on compressive strength :math:`\sigma_C` 
        at a range of confining pressures :math:`P` after :math:`N` wetting and drying cycles.
        
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
            ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
                                        containing data sets as :mod:`pandas` dataframes
            model_surface (str): 
                    key to dict (stored in :attr:`ed`) reference in  to 
                    2D  function regressed in 
                    :meth:`~.data.ExptData.fit_weathering_model` that models
                    joint dependence of weakness on proxy depth and proxy time
            text_label (list): text annotation as list of form (x-y coordinate, string, 
                            font size)
        """
        fig = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is not None:
            plt.title(title, fontdict={"fontsize": 11.5})

        axes = fig.add_subplot(1, 1, 1, projection = "3d",)
        
        df      = expt_data.ddict["li"]
        # TBD: need an exception here if model_surface not specified
        X,Y,Z   = expt_data.fdict[model_surface]

        axes.scatter(df.wetdryN, df.P, df.w_sigma2, color="k", s=40,)
        
        colors_list = ("y", "b")
        colors_array = np.empty(X.shape, dtype=str)
        rf = 3
        for y in range(len(X)):
            for x in range(len(Y)):
                colors_array[x, y] = colors_list[(x//rf+y//rf) % 2]
        axes.plot_surface(X,Y,Z,
                        facecolors=colors_array,  alpha=0.35,
                        linewidth=0, antialiased=True)
        plt.xlim(0,)
        plt.ylim(0,)
        axes.xaxis.pane.set_edgecolor("k")
        axes.yaxis.pane.set_edgecolor("k")
        axes.zaxis.pane.set_edgecolor("k")
        # axes.view_init(15, 50)
        axes.view_init(25, 72)
        axes.set_xlabel("Proxy time  $N$  [-]")
        axes.set_ylabel("Proxy depth  $P$ [MPa]")
        axes.set_zlabel("Weakness  $w=$  [-]")
        plt.grid(ls=":")



    def frontspeed_evolution(
            self, 
            name: str,
            title: str|None=None,
            nm: NumericalModel|None=None, 
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        """
        Plot time-evolution of dimensionless erosion rate :math:`\\nu(\\tau)`.
        
        Generate graph of numerical solutions :math:`j`  of dimensionless rock-surface erosion
        rate :math:`\\nu_i^j` over dimensionless time :math:`\\tau^j`. 
        These solutions are provided in the
        class instance :class:`~.solve1d.ErosionWeathering`.
        Legend-label by weathering number for this instance.
            
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
            ew (:class:`~.solve1d.ErosionWeathering`): 
                    instance of 1d weathering-mediated erosion model :mod:`~.solve1d` class
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        plt.plot(
            nm.tau_array, 
            nm.nu_array,  
            color="k", 
            lw=1, 
            label=r"${\mathcal{W}}=$"+f"{nm.W}",
        )
        plt.legend(loc="center right")
        plt.xlabel(r"Time  $\tau$  [-]")
        plt.ylabel(r"Front speed  $\partial\varphi/\partial\tau$  [-]")
        plt.grid(ls=":")

    def weakness_evolution(
            self, 
            name: str,
            title: str|None=None,
            nm: NumericalModel|None=None, 
            tc: float=40, 
            nd: int=2,
            text_label: str|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        r"""
        Plot 1d evolution of weathering profile :math:`{\omega}(\chi,\\tau)` 
        undergoing erosion.
        
        Generate graph of selected solutions :math:`j` 
        of 1d weathering-mediated erosion model
        weakness :math:`{\omega}_i^j = {\omega}(\chi_i,\\tau^j)`.
        This series of time slices 
        shows the propagation and development of a steady-state form
        of the weathering depth-profile
        as the rock surface is erodexpt_data. Quantities are all dimensionless.
        
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
            ew (:class:`~.solve1d.ErosionWeathering`): 
                    instance of 1d weathering-mediated erosion model :mod:`~.solve1d` class
            tc (int): 
                    :math:`\\tau^j`  slicing "rate"
            nd (int): 
                    number of decimal places in :math:`\\tau`  legend label
            text_label (list): text annotation as list of form (x-y coordinate, string, 
                            font size)
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        chi__ = nm.chi_array
        tau__ = nm.tau_array
        eta__ = nm.eta_array
        j__   = nm.j

        tau_slices1 = np.linspace(0,(nm.tau_n_steps-1)//tc,
                                num=5,endpoint=True,dtype=np.int64)
        tau_slices2 = np.linspace((nm.tau_n_steps-1)//tc*2 ,j__+1,
                                num=5,endpoint=True,dtype=np.int64)
        tau_slices = np.concatenate((tau_slices1,tau_slices2))
        cmap = plt.cm.brg.reversed()
        label=r"$\tau$={0:3."+str(nd)+"f}"
        for idx,tau_slice in enumerate(tau_slices):
            eta_slice = eta__[tau_slice]
            chi_front = chi__[eta_slice==0]
            chi_front = 0 if chi_front.shape[0]==0 else chi_front[-1]
            plt.plot(chi__[chi__>=chi_front],eta_slice[chi__>=chi_front],
                    color=cmap(idx/tau_slices.size),
                    label=label.format(tau__[tau_slice]))
            
        axes = plt.gca()
        plt.xlim(chi__[0],chi__[-1])
        plt.xlim(-(chi__[-1]-chi__[0])/30,chi__[-1]*1.08)
        x_limits = plt.xlim()
        y_limits = plt.ylim()
        plt.ylim(y_limits[0]/3,y_limits[1])
        bbox_props = dict(boxstyle="rarrow,pad=0.3", lw=1.5, 
                        fc="white", ec="k")
        t = axes.text((x_limits[1]-x_limits[0])*0.45, (y_limits[1]-y_limits[0])*0.5, 
                    "erosion", ha="right", va="center", rotation=0, color="k",
                    size=12, bbox=bbox_props)
        
        if text_label is not None:
            plt.text(*text_label[0], text_label[1], 
                    color="k", size=text_label[2],
                    verticalalignment="center", horizontalalignment="center",
                    transform=axes.transAxes)

        plt.legend(loc="upper right", fontsize=10,)
        plt.xlabel(r"Distance  $\chi$  [-]")
        plt.ylabel(r"Weakness  ${\omega}(\chi,\tau)$  [-]")
        plt.grid(ls=":")

    def stability_check(
            self, 
            name: str,
            title: str|None=None,
            tau: NDArray|None=None, 
            nu: NDArray|None=None, 
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:  
        """
        Visualize stability of numerical solution.
        
        Check numerical stability of time-stepping by plotting a close-up of
        the erosion front speed over time.
        
        Args:
            tau (numpy.ndarray): time slices :math:`\\tau^j` of numerical solution 
            nu (numpy.ndarray): speed of erosion front :math:`\\nu^j` 
                                at each time slice :math:`\\tau^j` 
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})

        plt.plot(tau,nu,"o-")
        plt.ylim(1.204,1.21);
        plt.xlim(4,4.03);

    def weakness_steadystate(
            self, 
            name: str,
            title: str|None=None,
            nm: NumericalModel|None=None, 
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        r"""
        Plot steady-state solution of weakness ${\omega}_s$.
        
        Graph the numerical solution 
        of the 1d weathering-mediated erosion model
        for weakness ${\omega}_s(\chi_i | {\mathcal{W}})$ 
        as a function of depth from the rock surface $\chi_i$
        for a given value of the weathering number 4{\mathcal{W}}$.
            
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>`
                                            figure 
            ew (:class:`~.solve1d.ErosionWeathering`): 
                    instance of 1d weathering-mediated erosion model :mod:`~.solve1d` class
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        j__ = (nm.j*3)//4
        i_offset = nm.n_chi_domain//7
        phi__  = nm.phi_array[j__]
        phi0__ = int(phi__/nm.Delta_chi)-i_offset
        chi_s__= nm.chi_array[phi0__:]-nm.chi_array[phi0__+i_offset]
        tau__  = nm.tau_array[-1]

        eta_s_numerical  = nm.eta_array[j__, phi0__:]
        eta_s_analytical = eta_chi_tau(chi_s__, tau__, nm.W)
        
        chi_front = chi_s__[eta_s_numerical==0][-1]
        plt.plot(
            chi_s__[chi_s__>=chi_front],
            eta_s_numerical[chi_s__>=chi_front],
            color="k", 
            lw=1, 
            label="numerical",
        )
        plt.plot(
            chi_s__[chi_s__>=chi_front],
            eta_s_analytical[chi_s__>=chi_front], 
            color="r", 
            lw=2, 
            label="analytical", 
            ls=(0, (4, 5)),
        )
        plt.xlim(nm.tau_array[0],nm.tau_array[-1])
        
        axes = plt.gca()
        plt.xlim((chi_s__[0],chi_s__[-1]))
        y_limits = plt.ylim()
        plt.ylim(y_limits[0]/3,y_limits[1])
        bbox_props = dict(
            boxstyle="rarrow,pad=0.3", 
            lw=1.5, 
            fc="white", 
            ec="DarkGreen",
        )
        _ = axes.text(
            -0.5, 
            (y_limits[1]-y_limits[0])/2, 
            "front motion", 
            ha="right", 
            va="center", 
            rotation=0, 
            color="DarkGreen",
            size=10, 
            bbox=bbox_props,
        )
        
        plt.legend()
        plt.xlabel(r"Distance relative to front  $\chi_s=\chi-\varphi_s$  [-]")
        plt.ylabel(r"Weakness  ${\omega}_s(\chi_s)$  [-]")
        plt.grid(ls=":")

    def weakness_steadystate_setW(
            self, 
            name: str,
            title: str|None=None,
            nms: Sequence|None=None, 
            chi_max: float=8,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:            
        r"""
        Plot a set of steady-state solutions of weakness ${\omega}_s$.
        
        Graph a set of numerical solution 
        of the 1d weathering-mediated erosion model
        for weakness ${\omega}_s(\chi_i | {\mathcal{W}})$
        as a function of depth from the rock surface $\chi_i$
        for a set of weathering numbers ${\mathcal{W}}$.
        
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>`
                                            figure 
            ew_list (list): 
                    list of instances of 1d weathering-mediated erosion model 
                    :mod:`~.solve1d` class :class:`~.solve1d.ErosionWeathering`
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})

        cmap = plt.cm.brg
        
        chi_min = 0
        for (idx, (nm, label,),) in enumerate(nms):
            j__ = (nm.j*3)//4
            i_offset = nm.n_chi_domain//10
            phi__  = nm.phi_array[j__]
            phi0__ = int(phi__/nm.Delta_chi)-i_offset
            chi_s__= nm.chi_array[phi0__:]-nm.chi_array[phi0__+i_offset]
            tau__  = nm.tau_array[-1]
            chi_min = min(chi_min,np.min(chi_s__))
            color=cmap(idx/len(nms))
            eta_s_numerical  = nm.eta_array[j__,phi0__:]
            chi_front = chi_s__[eta_s_numerical==0][-1]
            plt.plot(chi_s__[chi_s__>=chi_front],eta_s_numerical[chi_s__>=chi_front],
                    lw=1, color=color, label=label)
            
        axes = plt.gca()
        x_limits = plt.xlim()
        y_limits = plt.ylim()
        plt.ylim(y_limits[0]/3,y_limits[1])
        bbox_props = dict(
            boxstyle="rarrow,pad=0.3", 
            lw=1.5, 
            fc="white", 
            ec="DarkGreen",
        )
        _ = axes.text(
            -0.5, 
            (y_limits[1]-y_limits[0])/2, 
            "front motion", 
            ha="right", 
            va="center", 
            rotation=0, 
            color="DarkGreen",
            size=12, 
            bbox=bbox_props,
        )
        plt.text(
            chi_min/3,0.5, 
            "air", 
            color="k", 
            alpha=0.7, 
            size=14, 
            verticalalignment="center", 
            horizontalalignment="right",
        )
        plt.text(
            -chi_min/3,0.5, 
            "rock", 
            color="k", 
            alpha=0.7, 
            size=14,
            verticalalignment="center", 
            horizontalalignment="left",
        )
            
        plt.xlim((chi_min,chi_max))
        plt.legend(fontsize=11,)
        plt.xlabel(r"Distance relative to front  $\chi_s=\chi-\varphi_s$  [-]")
        plt.ylabel(r"Weakness  ${\omega}_s(\chi_s)$  [-]")
        plt.grid(ls=":")

    def weakness_steadystate_W(
            self, 
            name: str,
            title: str|None=None,
            eqns: Theory|None=None, 
            do_loglog: bool=True, 
            nms: Sequence|None=None, 
            text_label: str|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        r"""
        Plot the 1d model steady-state erosion rate $\nu_s$ 
        versus weathering number ${\mathcal{W}}$.
        
        Graph the functional dependence of dimensionless steady-state 
        erosion rate $\nu_s$ as a function of versus weathering number ${\mathcal{W}}$
        for the 1d weathering-mediated erosion model.
        The analytical solution is plotted as a black curve; numerical solutions are
        plotted as black circles; asymptotic behavior for low and high ${\mathcal{W}}$
        are shown as dashed lines. Explanatory annotations are includexpt_data.
        
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
            em (:class:`~.theory.WeatheringMediatedErosion`): 
                    instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
            do_loglog (bool): 
                flag whether to use log scales on both axes
            nus_solns_list (list): 
                set of numerical solutions of $\nu_s$
            text_label (list): text annotation as list of form (x-y coordinate, string, 
                            font size)
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        n_W_pts = 200
        W_array   = np.exp(np.linspace(np.log(0.01),np.log(50),n_W_pts))
        nus_array = np.array([eqns.nus_eqn_W.rhs.subs({W:W__}) for W__ in W_array])
        y_limits = (nus_array[0]*0.95,nus_array[-1])
    
        plt.plot(W_array, nus_array, color="k", lw=1.5, label="analytical")
        plt.autoscale(enable=True, tight=True)
        axes = plt.gca()
        if do_loglog:
            # Updated 2040/5/29 - turned off clipping since arg has been deprecated?
            axes.set_xscale("log",) # nonposx="clip")
            axes.set_yscale("log",) # nonposy="clip")
            axes.yaxis.set_major_formatter(ticker.FormatStrFormatter("%0.0f"))
            axes.yaxis.set_minor_formatter(ticker.FormatStrFormatter("%0.0f"))
            plt.plot(
                (0.25,0.25), (y_limits[0],y_limits[1]/3), color="gray", ls=":",
            )
            plt.plot(
                (2.63,2.63), (y_limits[0],y_limits[1]), color="gray", ls=":",
            )
    #         plt.ylim(*y_limits)
            plt.text(
                0.2, 0.23, 
                r"low ${\mathcal{W}}$", 
                color="brown",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            plt.text(
                0.2, 0.13, 
                r"$\nu_\mathsf{{s}} \approx 1 + {\mathcal{W}}$", 
                color="brown",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            # plt.text(0.2, 0.15, r"$v_s \approx {u_0} + \dfrac{w_0}{k}$", color="brown",
            #          verticalalignment="center", horizontalalignment="center",
            #          transform=axes.transAxes)
            plt.text(
                0.52, 0.43, 
                r"transitional ${\mathcal{W}}$", 
                color="gray",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            plt.text(
                0.78, 0.9, 
                r"high ${\mathcal{W}}$", 
                color="blue",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            plt.text(
                0.78, 0.8, 
                r"$\nu_\mathsf{{s}} \approx \dfrac{1}{2}+\sqrt{{\mathcal{W}}}$", 
                color="blue",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            # plt.text(0.81,0.92, r"$v_s \approx \dfrac{{u_0}}{2}+\sqrt{\dfrac{{u_0} w_0}{k}}$", 
            #          color="blue",
            #          verticalalignment="center", horizontalalignment="center",
            #          transform=axes.transAxes)
            
        if nms is not None:
            for idx,nus_soln in enumerate(nms):
                plt.plot(
                    nus_soln.W,
                    nus_soln.nu_s,
                    "o",
                    c="k",
                    label=("numerical" if idx==0 else None),
                )
        
        plt.plot(
            W_array[W_array<0.7],
            1+W_array[W_array<0.7], 
            label=r"low ${\mathcal{W}}$ approx",  
            ls="--",
            c="brown",
        )
        plt.plot(
            W_array[W_array>1],
            0.5+np.sqrt(W_array[W_array>1]), 
            label=r"high ${\mathcal{W}}$ approx", 
            ls="--",
            c="blue",
        )
        if text_label is not None:
            plt.text(
                0.82,0.25, 
                text_label, 
                color="k", 
                size=14,
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
        
        plt.legend(loc="upper left", fontsize=10,)
        plt.xlabel(r"Weathering number  ${\mathcal{W}}$  [-]")
        plt.ylabel(r"Erosion rate  ${\omega}_\mathsf{{s}}$  [-]")
        # plt.grid(ls=":")

    def weakness_steadystate_W_transition(
            self, 
            name: str,
            title: str|None=None,
            eqns: Theory|None=None, 
            text_label: str|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        """
        Plot the steady-state erosion rate :math:`\\nu_s` relative to its asymptotic behavior.
        
        Visualize the transitional behavior of the steady-state erosion rate :math:`\\nu_s`
        at intermediate weathering numbers :math:`W` by plotting how asymptotes at
        low and high :math:`W` deviate from the full model.
        
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
            em (:class:`~.theory.WeatheringMediatedErosion`): 
                    instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
            text_label (list): text annotation as list of form (x-y coordinate, string, 
                            font size)
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        n_W_pts = 200
        W_array = np.exp(np.linspace(np.log(0.015),np.log(50),n_W_pts))
        nus_array = np.array([eqns.nus_eqn_W.rhs.subs({W:W__}) for W__ in W_array])
        plt.plot(W_array, 1+0*W_array, color="k", lw=1.5)
        axes = plt.gca()
        plt.plot(W_array[W_array<0.48], 
            nus_array[W_array<0.48]/(1+(W_array[W_array<0.48])),
                color="brown", ls="-", lw=1.5, label="low ${\mathcal{W}}$ approx")
        plt.plot(W_array[W_array>1], 
                nus_array[W_array>1]/(0.5+np.sqrt(W_array[W_array>1])),
                color="b", ls="-", lw=1.5, label="high ${\mathcal{W}}$ approx")
        plt.text(0.18,0.18, r"low ${\mathcal{W}}$", color="brown",
                verticalalignment="center", horizontalalignment="center",
                transform=axes.transAxes)
        plt.text(0.18,0.26, r"$\nu_\mathsf{{s}} \approx 1 + {\mathcal{W}}$", color="brown",
                verticalalignment="center", horizontalalignment="center",
                transform=axes.transAxes)
        plt.text(0.5,0.53, r"transitional ${\mathcal{W}}$", color="gray",
                verticalalignment="center", horizontalalignment="center",
                transform=axes.transAxes)
        plt.text(0.82,0.85, r"high ${\mathcal{W}}$", color="blue",
                verticalalignment="center", horizontalalignment="center",
                transform=axes.transAxes)
        plt.text(0.82,0.75, r"$\nu_\mathsf{{s}} \approx \dfrac{1}{2}+\sqrt{{\mathcal{W}}}$", 
                color="blue",
                verticalalignment="center", horizontalalignment="center",
                transform=axes.transAxes)
        if text_label is not None:
            plt.text(0.82,0.25, text_label, 
                    color="k", size=14,
                    verticalalignment="center", horizontalalignment="center",
                    transform=axes.transAxes)
        plt.autoscale(enable=True, tight=True)
        axes.set_xscale("log",) #nonposx="clip")
        y_limits = (0.9,1.1)
        plt.ylim(*y_limits)
        axes.set_yticks(np.linspace(*y_limits,5))
        axes.yaxis.set_major_formatter(ticker.FormatStrFormatter("%0.2f"))
        plt.plot((0.25,0.25),y_limits, color="gray", ls=":")
        plt.plot((2.63,2.63),y_limits, color="gray", ls=":")
        plt.legend(loc="upper left")
        plt.xlabel("Weathering number  ${\mathcal{W}}$  [-]")
        plt.ylabel(
            r"Approx erosion rate deviation  "
            + r"$\nu_\mathsf{{s}}^\mathrm{apx}/\nu_\mathsf{{s}}$  [-]"
        )
        # plt.grid(ls=":")



    def referosionrate_vs_refweatheringrate(
            self, 
            name: str,
            title: str|None=None,
            eqns: Theory|None=None, 
            k__: float=1, 
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        """
        Plot baseline erosion rate versus baseline weathering rate.
            
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
            em (:class:`~.theory.WeatheringMediatedErosion`): 
                    instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
            k__ (float): 
                weathering-weakening rate profile reciprocal e-folding depth
            text_label (list): text annotation as list of form (x-y coordinate, string, 
                            font size)
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        w0_array = 10**np.linspace(-3,+1,100)
        for vs__ in [0.1,1,3]:
            v0_array = np.array(
                [sy.N(eqns.v0_eqn_vs_w0.rhs.subs({v_s:vs__, w_0:w0__,k:k__})) 
                                    for w0__ in w0_array])
            plt.plot(v0_array,w0_array, label="$v_s=${}".format(vs__))
        plt.ylabel(r"$w_0$")
        plt.xlabel(r"${u_0}$")
        plt.legend(loc="lower right")
        plt.grid(ls=":")

    def referosionrate_vs_refweatheringrate_steadystate(
            self, 
            name: str,
            title: str|None=None,
            eqns: Theory|None=None, 
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        """
        Plot baseline erosion rate versus surface weakness at steady-state.
            
        Args:
            fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
                reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
            em (:class:`~.theory.WeatheringMediatedErosion`): 
                    instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
            text_label (list): text annotation as list of form (x-y coordinate, string, 
                            font size)
        """
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        etas0_array = 10**np.linspace(-2,+1,100)
        for vs__ in [0.5,1,1.5]:
            v0_array = np.array(
                [sy.N(eqns.v0_eqn_etas0_vs.rhs.subs({v_s:vs__,eta_s0:etas0__})) 
                                    for etas0__ in etas0_array])
            plt.plot(etas0_array,v0_array, label="$v_s=${}".format(vs__))
        plt.xlabel(r"Surface weakness (degree of weathering)  ${\omega}_{s0}$")
        plt.ylabel(r"Baseline (potential) erosion rate  ${u_0}$")
        plt.gca().invert_yaxis()
        plt.xlim(0,3)
        plt.ylim(5,0)
        plt.legend(loc="lower right")
        plt.grid(ls=":")



    def channel_generic(
            self, 
            name: str,
            title: str|None=None,
            zys: Sequence|None=None, 
            do_equal_aspect=False,
            text_labels: Sequence|None=None,
            fig_size: tuple[float,float]=(6,4,),
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
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        zy=zys[0]
        plt.plot(zy[2], zy[0], label=zy[4], color="k")
        plt.ylabel(zy[1])
        plt.xlabel(zy[3])
        axes = plt.gca()
        if do_equal_aspect:
            axes.set_aspect("equal")
        else:
            pass
        if text_labels is not None:
            for text_label in text_labels:
                plt.text(*text_label[0], text_label[1], 
                        color=text_label[3], size=text_label[2],
                        verticalalignment="center", horizontalalignment="center",
                        transform=axes.transAxes, rotation=text_label[4])
        plt.grid("on",ls=":")
        if len(zys)>=2:
            zy=zys[1]
            plt.plot(0,0, label=zy[4], color="forestgreen")
        plt.legend()
        
        if len(zys)>=2:
            zy=zys[1]
            alt_axes = axes.twiny()
            alt_axes.plot(zy[2], zy[0], label=zy[4], color="forestgreen")
            alt_axes.set_xlabel(zy[3], color="forestgreen")
        plt.grid(ls=":")

    def channel_refweatheringrate_referosionrate_W(
            self, 
            name: str,
            title: str|None=None,
            model: Any|None=None, 
            text_label: str|None=None,
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
        _ = self.create_figure(fig_name=name, fig_size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        

        plt.plot(model.w0_array/model.pdict[k], model.z_array, label="$w_0/k$")
        plt.plot(model.v0_array, model.z_array, label="${u_0}$")
        x_limits = plt.xlim()
        y_limits = plt.ylim()
        plt.plot(
            [-1,-1], label="${\mathcal{W}}=w_0/{u_0} k$", color="forestgreen",
        )
        plt.plot(
            model.vs_array, model.z_array, label="$v_s$", color="k", lw=2,
        )
        plt.xlim(x_limits)
        plt.ylim(y_limits)
        plt.xlabel(r"Speeds $w_0(z)/k$, ${u_0}(z)$, $v_s(z)$")
        plt.ylabel(r"Height above bed  $z$")
        plt.legend(loc="upper center")
        plt.grid("on",ls=":")

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
            plt.text(*text_label[0], text_label[1], 
                    color="k", size=text_label[2],
                    verticalalignment="center", horizontalalignment="center",
                    transform=axes.transAxes)

