"""
Visualization classes.
"""
import warnings
import numpy as np
import matplotlib.pyplot as plt

from typing import Any
from collections.abc import Sequence
from pandas import DataFrame

from wmbe.model import WeatheringMediatedWeakness, linear_model
from wmbe.viz_base import VizBase

warnings.filterwarnings("ignore")

__all__ = [
    "DataViz",
]

class VizData(VizBase):
    """
    Data visualization class.
    """
    def inoue_w_wetdryN(
            self, 
            name: str,
            title: str|None=None,
            data: DataFrame|None=None, 
            model: WeatheringMediatedWeakness|None=None,
            text_label: Sequence|None=None,
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
        """
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        #     ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
        #                                 containing data sets as :mod:`pandas` dataframes
        #     text_label (list): text annotation as list of form (x-y coordinate, string, 
        #                     font size)

        _ = self.create_figure(name=name, size=fig_size,)
        if title is not None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        # sigmaT  = df.sigmaT
        wetdryN = data.wetdryN
        erodibility_sigma2   = data.w_sigma2
        # erodibility_sigma1p5 = df.w_sigma1p5
        erodibility_sigma2_fit = model.fits["default"][0]
            
        plt.plot(
            wetdryN,
            linear_model(wetdryN, *erodibility_sigma2_fit,),
            color="k",
            label=r"$w \sim 1/\sigma_T^2$",
        )
        plt.errorbar(
            np.unique(wetdryN),
            model.w_s2_means, 
            label="",
            xerr=None,
            yerr=model.w_s2_stds*1,
            ecolor="k", 
            mec="w", 
            color="w", 
            fillstyle="full", 
            alpha=1,
            fmt="o", 
            markersize=0, 
            markeredgewidth=2,
            elinewidth=1.5,
            capthick=3,
            capsize=7,
        )
        plt.plot(
            np.unique(wetdryN),
            model.w_s2_means, 
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
            color="orange", 
            alpha=0.7,
            ls="", 
            marker="s",
            markeredgecolor="k",
            ms=5,
        )

        plt.legend(loc="upper left")
    #     plt.ylim(0,)
        plt.xlabel(
            r"Proxy time (no. wet/dry cycles)  $N$  [-]"
        )
        plt.ylabel(
            r"Weakness  $w=(\sigma_T\left[N\right]/\sigma_{\mathrm{ref}})^{-2}$  [-]"
        )
        # plt.ylabel(r"Weakness  [-]")

        if text_label is not None:
            axes = plt.gca()
            plt.text(
                *text_label[0], 
                text_label[1], 
                color="k", 
                size=text_label[2],
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
        plt.grid(ls=":")
            
    def li_w_wetdryN(
            self, 
            name: str,
            title: str|None=None,
            data: DataFrame|None=None, 
            model: WeatheringMediatedWeakness|None=None,
            text_label: Sequence|None=None,
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
        """
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        #     ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
        #                                 containing data sets as :mod:`pandas` dataframes
        #     text_label (list): text annotation as list of form (x-y coordinate, string, 
        #                     font size)

        _ = self.create_figure(name=name, size=fig_size,)
        if title is not None:
            plt.title(title, fontdict={"fontsize": 11.5})

        color_list = (
            "darkblue",
            "darkmagenta",
            "firebrick",
            "red",
            "orange",
            "green",
        )
        marker_list = (
            "o", "s", "v", "^", "8", "+",
        )
        n_cols = len(color_list)
        for (idx, P_,) in enumerate(np.unique(data.P)):
            selection_name = f"{'P'}_{P_}"
            # sigma  = df.sigmaC[df.P==P_]/100
            wetdryN = data.wetdryN[data.P==P_]
            erodibility_sigma   = data.w_sigma2[data.P==P_]
            erodibility_sigma_fit = model.fits[selection_name][0]
            
            plt.plot(
                wetdryN,
                linear_model(wetdryN, *erodibility_sigma_fit,),
                color=color_list[idx%n_cols],
                label="",
            )
            plt.errorbar(
                wetdryN,
                erodibility_sigma, 
                label=r"$P = ${}$\,$MPa".format(P_),
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
        
        plt.legend(loc="upper left")
        plt.ylim(0,)
        plt.xlabel(
            "Proxy time (no. wet/dry cycles)  $N$  [-]"
        )
        plt.ylabel(
            r"Weakness  $w=(\sigma_C\left[N\right]/\sigma_\mathrm{ref})^{-2}$  [-]"
        )

        if text_label is not None:
            axes = plt.gca()
            plt.text(
                *text_label[0], text_label[1], 
                color="k", 
                size=text_label[2],
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
        plt.grid(ls=":")
            
    def li_w_P(
            self, 
            name: str,
            title: str|None=None,
            data: DataFrame|None=None, 
            model: WeatheringMediatedWeakness|None=None,
            text_label: Sequence|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        r"""
        Plot Li et al  data on weakness :math:`w` versus proxy depth.
        
        Generate graph of rock weakness :math:`w` versus proxy depth 
        inferred from  `Li et al (2016)`_ data on compressive strength :math:`\sigma_C` 
        at a range of confining pressures :math:`P` after :math:`N` 
        wetting and drying cycles.
        """
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        #     ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
        #                                 containing data sets as :mod:`pandas` dataframes
        #     text_label (list): text annotation as list of form (x-y coordinate, string, 
        #                     font size)
        _ = self.create_figure(name=name, size=fig_size,)
        if title is not None:
            plt.title(title, fontdict={"fontsize": 11.5})

        color_list = (
            "darkblue",
            "darkmagenta",
            "firebrick",
            "red",
            "orange",
            "green",
        )
        marker_list = (
            "o", "s", "v", "^", "P", "p",
        )
        n_cols = len(color_list)
        
        sampled_fit = model.fits2d["default"]
        
        for idx,wetdryN in enumerate(np.flip(np.unique(data.wetdryN))):
            # wdN_ = df.wetdryN[df.wetdryN==wetdryN]
            P_ = data.P[data.wetdryN==wetdryN]
            w_ = data.w_sigma2[data.wetdryN==wetdryN]
            P_fit = sampled_fit[1]
            w_fit = np.flipud(sampled_fit[2].T)[idx]
            plt.plot(
                P_fit, w_fit, color=color_list[idx%n_cols],
            )
            plt.errorbar(
                P_, 
                w_,
                label=f"N = {wetdryN}",
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
            
        plt.legend(loc="upper right")
        plt.ylim(0.,)
        plt.autoscale(enable=True, tight=True, axis="x",)
        x_limits = plt.xlim()
        plt.plot(
            x_limits, (1,1), color="gray", ls=":",
        )
        plt.xlabel("Proxy depth (confining pressure $P$)  [MPa]")
        plt.ylabel(r"Weakness  $w=(\sigma_C[N]/\sigma_\mathrm{ref})^{-2}$  [-]")

        if text_label is not None:
            axes = plt.gca()
            plt.text(
                *text_label[0], 
                text_label[1], 
                color="k", 
                size=text_label[2],
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
        plt.grid(ls=":")
                
    def li_w_normed_P(
            self, 
            name: str,
            title: str|None=None,
            data: DataFrame|None=None, 
            model: WeatheringMediatedWeakness|None=None,
            text_label: Sequence|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        r"""
        Plot Li et al data on weakness :math:`w` (normalized using 2D model) 
        versus proxy depth.
        
        Generate graph of rock weakness versus proxy depth
        inferred from  `Li et al (2016)`_
        data on compressive strength :math:`\sigma_C` 
        at a range of confining pressures :math:`P` after :math:`N` wetting and drying cycles.
        """        
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        #     ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
        #                                 containing data sets as :mod:`pandas` dataframes
        #     text_label (list): text annotation as list of form (x-y coordinate, string, 
        #                     font size)

        _ = self.create_figure(name=name, size=fig_size,)
        if title is not None:
            plt.title(title, fontdict={"fontsize": 11.5})

        color_list = (
            "darkblue",
            "darkmagenta",
            "firebrick",
            "red",
            "orange",
            "green",
        )
        marker_list = (
            "o", "s", "v", "^", "P", "p",
        )
        n_cols = len(color_list)
        
        sampled_fit = model.fits2d["default"]
        w_ref_vec = np.flipud(sampled_fit[2].T)[:,0]-1
        
        P_fit = sampled_fit[1]
        w_fit = (np.flipud(sampled_fit[2].T)[-1]-1)/w_ref_vec[-1]+1

        plt.errorbar(
            np.unique(data.P),
            model.w_s2normed_means, 
            label="mean data",
            xerr=None,
            yerr=model.w_s2normed_stds*2,
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
        
        for idx,wetdryN in enumerate(np.flip(np.unique(data.wetdryN))):
            # wdN_ = df.wetdryN[df.wetdryN==wetdryN]
            P_ = data.P[data.wetdryN==wetdryN]
            w_normed = data.w_s2normed[data.wetdryN==wetdryN]
            plt.errorbar(
                P_, 
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
                    
        plt.plot(
            P_fit,
            w_fit, 
            color="darkblue", 
            label=r"$w \sim \exp(-k\chi)$",
        )

        plt.legend(loc="upper right")
        plt.ylim(0.8,)
        plt.autoscale(enable=True, tight=True, axis="x")
        x_limits = plt.xlim()
        plt.plot(
            x_limits,(1,1), color="gray", ls=":",
        )
        plt.xlabel("Proxy depth (confining pressure $P$)  [MPa]")
        plt.ylabel("Normalized weakness  $w(\\tau,\\chi)/w(0,\\tau)$  [-]")
        
        if text_label is not None:
            axes = plt.gca()
            plt.text(
                *text_label[0], 
                text_label[1], 
                color="k", 
                size=text_label[2],
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
        plt.grid(ls=":")
            
    def li_w_wetdryN_P(
            self, 
            name: str,
            title: str|None=None,
            data: DataFrame|None=None, 
            model: WeatheringMediatedWeakness|None=None,
            surface: str|None=None,
            # text_label: Sequence|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        r"""
        Plot Li et al data in 3D on weakness versus proxy depth and confining pressure.
        
        Generate 3D view of 2D surface model & data
        of rock weakness versus proxy depth inferred from 
        `Li et al (2016)`_ data on compressive strength :math:`\sigma_C` 
        at a range of confining pressures :math:`P` after :math:`N` wetting and drying cycles.
        """        
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        #     ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
        #                                 containing data sets as :mod:`pandas` dataframes
        #     model_surface (str): 
        #             key to dict (stored in :attr:`ed`) reference in  to 
        #             2D  function regressed in 
        #             :meth:`~.data.ExptData.fit_weathering_model` that models
        #             joint dependence of weakness on proxy depth and proxy time
        #     text_label (list): text annotation as list of form (x-y coordinate, string, 
        #                     font size)
        fig = self.create_figure(name=name, size=fig_size,)
        if title is not None:
            plt.title(title, fontdict={"fontsize": 11.5})

        axes = fig.add_subplot(1, 1, 1, projection = "3d",)
        
        # TBD: need an exception here if model_surface not specified
        (X, Y, Z,) = model.fits[surface]

        axes.scatter(
            data.wetdryN, data.P, data.w_sigma2, color="k", s=40,
        )
        
        colors_list = ("y", "b",)
        colors_array = np.empty(X.shape, dtype=str,)
        rf = 3
        for y in range(len(X)):
            for x in range(len(Y)):
                colors_array[x, y] = colors_list[(x//rf+y//rf) % 2]
        axes.plot_surface(
            X,
            Y,
            Z,
            facecolors=colors_array,  
            alpha=0.35,
            linewidth=0, 
            antialiased=True,
        )
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