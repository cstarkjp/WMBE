"""
Visualization classes.
"""
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from typing import Any
from collections.abc import Sequence
from numpy.typing import NDArray

from wmbe.numerical.solve1d import weakness_chi_tau, NumericalModel
from wmbe.theory.symbols import *
from wmbe.theory.equations import Equations
from wmbe.viz.base import VizBase

warnings.filterwarnings("ignore")

__all__ = [
    "VizSim",
]

class VizSimulations(VizBase):
    """
    Numerical simulation visualization class.
    """
    def front_speed_evolution(
            self, 
            name: str,
            title: str|None=None,
            numerical_model: NumericalModel|None=None, 
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        """
        Plot time-evolution of dimensionless erosion rate $\\nu(\\tau)$.
        
        Generate graph of numerical solutions $j$ of dimensionless rock-surface 
        erosion rate $\\nu_i^j$ over dimensionless time $\\tau^j$. 
        These solutions are provided in the class instance 
        `ErosionWeathering`.
        Legend-label by weathering number for this instance.

        Args:
            name: key for figure dictionary
            title: optional title for figure
            nm: instance of 1d weathering-mediated erosion NumericalModel class
            fig_size: tuple (x,y) in inches (!)
        """
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        plt.plot(
            numerical_model.τ_array, 
            numerical_model.ν_array,  
            color="k", 
            lw=1, 
            label=r"${\mathcal{W}}=$"+f"{numerical_model.W}",
        )
        plt.legend(loc="center right")
        plt.xlabel(r"Time  $\tau$  [-]")
        plt.ylabel(r"Front speed  $\partial\phi/\partial\tau$  [-]")
        plt.grid(ls=":")

    def weakness_evolution(
            self, 
            name: str,
            title: str|None=None,
            numerical_model: NumericalModel|None=None, 
            tc: float=40, 
            nd: int=2,
            do_moving_frame: bool=False,
            do_zoom: bool=False,
            text_label: Sequence|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        """
        Plot 1d evolution of weathering profile ${\\omega}(\\chi,\\tau)$
        undergoing erosion.
        
        Generate graph of selected solutions $j$ of 1d weathering-mediated 
        erosion model weakness ${\\omega}_i^j = {\\omega}(\\chi_i,\\tau^j)$. 
        This series of time slices shows the propagation and development of 
        a steady-state form of the weathering depth-profile as the rock surface
        is eroded. Quantities are all dimensionless.

        Args:
            name: key for figure dictionary
            title: optional title for figure
            nm: instance of 1d weathering-mediated erosion NumericalModel class
            tc: $\\tau^j$  slicing "rate"
            nd: number of decimal places in $\\tau$ legend label
            do_moving_frame: shift curves into frame of moving front
            do_zoom: reduce y axis range
            text_label: 
                text annotation as tuple of form 
                (x-y coordinate, string, font size)
            fig_size: tuple (x,y) in inches (!)
        """
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        χ: NDArray = numerical_model.χ_array
        τ: NDArray = numerical_model.τ_array
        ω: NDArray = numerical_model.ω_array
        j: int = numerical_model.j

        τ_slices1: NDArray = np.linspace(
            0, 
            (numerical_model.τ_n_steps-1)//tc,
            num=5,
            endpoint=True,
            dtype=np.int64,
        )
        τ_slices2: NDArray = np.linspace(
            (numerical_model.τ_n_steps-1)//tc*2,
            j+1,
            num=5,
            endpoint=True,
            dtype=np.int64,
        )
        τ_slices: NDArray = np.concatenate((τ_slices1, τ_slices2,))
        cmap: Any = plt.cm.brg.reversed()
        label: str = r"$\tau$={0:3."+str(nd)+"f}"
        for (idx, τ_slice,) in enumerate(τ_slices):
            ω_slice: NDArray = ω[τ_slice]
            i_front: int = (ω_slice>0).argmax()
            i_front = i_front if i_front==0 else i_front-1
            χ_front: NDArray = (
                χ - χ[i_front] if do_moving_frame
                else χ
            )
            plt.plot(
                χ_front[i_front:],
                ω_slice[i_front:],
                # χ[χ>=χ_front],
                # ω_slice[χ>=χ_front],
                color=cmap(idx/τ_slices.size),
                label=label.format(τ[τ_slice]),
            )
            
        axes = plt.gca()
        if do_moving_frame:
            plt.xlim(
                -(χ[-1]-χ[0])/45 if do_zoom else -(χ[-1]-χ[0])/30, 
                χ[-1]*1.08/6 if do_zoom else χ[-1]*1.08/4,
            )
        else:
            plt.xlim(-(χ[-1]-χ[0])/30, χ[-1]*1.08,)
        x_limits = plt.xlim()
        y_limits = plt.ylim()
        plt.ylim(
            0.99 if do_moving_frame and do_zoom else y_limits[0]/3, 
            1.027 if do_moving_frame and do_zoom else y_limits[1],
        )

        if not do_moving_frame:
            bbox_props = dict(
                boxstyle="rarrow,pad=0.3", 
                lw=1.5, 
                fc="white", 
                ec="k",
            )
            t = axes.text(
                (x_limits[1]-x_limits[0])*0.45, 
                (y_limits[1]-y_limits[0])*0.5, 
                "erosion", 
                ha="right", 
                va="center", 
                rotation=0, 
                color="k",
                size=12, 
                bbox=bbox_props,
            )
        
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

        plt.legend(loc="upper right", fontsize=10, framealpha=1,)
        if do_moving_frame:
            x_label = r"Distance from moving surface  $\chi -{\phi}$  [-]"
            y_label = r"Weakness  ${\omega}(\chi, \tau)$  [-]"
        else:
            x_label = r"Distance  $\chi$  [-]"
            y_label = r"Weakness  ${\omega}(\chi,\tau)$  [-]"
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.grid(ls=":")

    def stability_check(
            self, 
            name: str,
            title: str|None=None,
            τ: NDArray|None=None, 
            ν: NDArray|None=None, 
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:  
        """
        Visualize stability of numerical solution.
        
        Check numerical stability of time-stepping by plotting a close-up of
        the erosion front speed over time.

        Args:
            name: key for figure dictionary
            title: optional title for figure
            tau: time slices $\\tau^j$ of numerical solution
            nu: speed of erosion front $\\nu^j$ at each 
                        time slice $\\tau^j$
            fig_size: tuple (x,y) in inches (!)
        """
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})

        plt.plot(τ, ν, "o-",)
        plt.ylim(1.204, 1.21,);
        plt.xlim(4, 4.03,);

    def weakness_steadystate(
            self, 
            name: str,
            title: str|None=None,
            numerical_model: NumericalModel|None=None, 
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        """
        Plot steady-state solution of weakness ${\\omega}_s$.
        
        Graph the numerical solution of the 1d weathering-mediated erosion model
        for weakness ${\\omega}_s(\\chi_i | {\\mathcal{W}})$ as a function of 
        depth from the rock surface $\\chi_i$ for a given value of the 
        weathering number ${\\mathcal{W}}$.

        Args:
            name: key for figure dictionary
            title: optional title for figure
            nm: instance of 1d weathering-mediated erosion NumericalModel class
            fig_size: tuple (x,y) in inches (!)
        """
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        j_ = (numerical_model.j*3)//4
        i_offset = numerical_model.n_χ_domain//7
        phi_  = numerical_model.φ_array[j_]
        phi0_ = int(phi_/numerical_model.Δχ)-i_offset
        chi_s_= (
            numerical_model.χ_array[phi0_:]
            - numerical_model.χ_array[phi0_+i_offset]
        )
        tau_  = numerical_model.τ_array[-1]

        eta_s_numerical  = numerical_model.ω_array[j_, phi0_:]
        eta_s_analytical = weakness_chi_tau(chi_s_, tau_, numerical_model.W)
        
        chi_front = chi_s_[eta_s_numerical==0][-1]
        plt.plot(
            chi_s_[chi_s_>=chi_front],
            eta_s_numerical[chi_s_>=chi_front],
            color="k", 
            lw=1, 
            label="numerical",
        )
        plt.plot(
            chi_s_[chi_s_>=chi_front],
            eta_s_analytical[chi_s_>=chi_front], 
            color="r", 
            lw=2, 
            label="analytical", 
            ls=(0, (4, 5),),
        )
        plt.xlim(numerical_model.τ_array[0], numerical_model.τ_array[-1],)
        
        axes = plt.gca()
        plt.xlim((chi_s_[0],chi_s_[-1]))
        y_limits = plt.ylim()
        plt.ylim(y_limits[0]/3,y_limits[1])
        bbox_props = dict(
            boxstyle="rarrow,pad=0.3", 
            lw=1.5, 
            fc="white", 
            ec="DarkGreen",
        )
        axes.text(
            -0.5, 
            (y_limits[1]-y_limits[0])/2, 
            "surface erosion", 
            ha="right", 
            va="center", 
            rotation=0, 
            color="DarkGreen",
            size=10, 
            bbox=bbox_props,
        )
        plt.legend()
        plt.xlabel(
            r"Distance relative to erosion surface  $\chi_s=\chi-\phi_s$  [-]"
        )
        plt.ylabel(r"Weakness  ${\omega}_s(\chi_s)$  [-]")
        plt.grid(ls=":")

    def weakness_steadystate_setW(
            self, 
            name: str,
            title: str|None=None,
            numerical_models: Sequence|None=None, 
            chi_max: float=8,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:            
        """
        Plot a set of steady-state solutions of weakness ${\\omega}_s$.
        
        Graph a set of numerical solution of the 1d weathering-mediated erosion
        model for weakness ${\\omega}_s(\\chi_i | {\\mathcal{W}})$ as a function
        of depth from the rock surface $\\chi_i$ for a set of weathering numbers
        ${\\mathcal{W}}$.

        Args:
            name: key for figure dictionary
            title: optional title for figure
            numerical_models: 
                sequence of instances of 1d weathering-mediated erosion 
                NumericalModel class
            chi_max: maximum $\\chi$ on x axis
            fig_size: tuple (x,y) in inches (!)
        """
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})

        cmap = plt.cm.brg
        
        chi_min = 0
        for (idx, (nm, label,),) in enumerate(numerical_models):
            j_ = (nm.j*3)//4
            i_offset = nm.n_χ_domain//15
            phi_  = nm.φ_array[j_]
            phi0_ = int(phi_/nm.Δχ)-i_offset
            chi_s_= nm.χ_array[phi0_:]-nm.χ_array[phi0_+i_offset]
            # tau_  = nm.τ_array[-1]
            chi_min = min(chi_min, np.min(chi_s_))
            color=cmap(idx/len(numerical_models))
            eta_s_numerical  = nm.ω_array[j_, phi0_:]
            chi_front = chi_s_[eta_s_numerical==0][-1]
            plt.plot(
                chi_s_[chi_s_>=chi_front],
                eta_s_numerical[chi_s_>=chi_front],
                lw=1, 
                color=color, 
                label=label,
            )
            
        axes = plt.gca()
        # x_limits = plt.xlim()
        y_limits = plt.ylim()
        plt.ylim(y_limits[0]/3, y_limits[1],)
        bbox_props = dict(
            boxstyle="rarrow,pad=0.3", 
            lw=1.5, 
            fc="white", 
            ec="DarkGreen",
        )
        axes.text(
            -0.5, 
            (y_limits[1]-y_limits[0])/2, 
            "erosion", 
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
            color="dimgrey", 
            alpha=0.7, 
            size=14, 
            verticalalignment="center", 
            horizontalalignment="right",
        )
        plt.text(
            -chi_min/3,0.5, 
            "rock", 
            color="dimgrey", 
            alpha=0.7, 
            size=14,
            verticalalignment="center", 
            horizontalalignment="left",
        )
            
        plt.xlim(chi_min, chi_max,)
        plt.legend(fontsize=11,)
        plt.xlabel(
            r"Distance relative to erosion surface  $\chi_s=\chi-\phi_s$  [-]"
        )
        plt.ylabel(r"Weakness  ${\omega}_s(\chi_s)$  [-]")
        plt.grid(ls=":")

    def erosion_rate_steadystate_W(
            self, 
            name: str,
            title: str|None=None,
            eqns: Equations|None=None, 
            do_loglog: bool=True, 
            numerical_models: Sequence|None=None, 
            text_label: Sequence|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        """
        Plot the 1d model steady-state erosion rate $\\nu_s$  versus weathering
        number ${\\mathcal{W}}$.
        
        Graph the functional dependence of dimensionless steady-state erosion 
        rate $\\nu_s$ as a function of versus weathering number ${\\mathcal{W}}$
        for the 1d weathering-mediated erosion model. The analytical solution 
        is plotted as a black curve; numerical solutions are plotted as black
        circles; asymptotic behavior for low and high ${\\mathcal{W}}$
        are shown as dashed lines. 

        Args:
            name: key for figure dictionary
            title: optional title for figure
            eqns: 
                instance of 1d weathering-mediated erosion theory class
            do_loglog: 
                flag whether to use log scales on both axes
            numerical_models: 
                sequence of numerical solutions of $\\nu_s$
            text_label: 
                text annotation as tuple of form 
                (x-y coordinate, string, font size)
            fig_size: tuple (x,y) in inches (!)
        """
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        n_W_pts: int = 200
        W_array: NDArray = np.exp(
            np.linspace(np.log(0.01), np.log(50), n_W_pts,)
        )
        ν_s_array: NDArray = np.array([
            eqns.nus_eqn_W.rhs.subs({W: W_}) 
            for W_ in W_array
        ])
        y_limits: tuple[float,float] = (ν_s_array[0]*0.95, ν_s_array[-1],)
    
        plt.plot(W_array, ν_s_array, color="k", lw=1.5, label="analytical",)
        axes = plt.gca()
        axes.autoscale(enable=True, tight=True)
        if do_loglog:
            axes.set_xscale("log",)
            axes.set_yscale("log",)
            axes.yaxis.set_major_formatter(ticker.FormatStrFormatter("%0.0f"))
            axes.yaxis.set_minor_formatter(ticker.FormatStrFormatter("%0.0f"))
            # transition line lower W
            plt.plot(
                (0.25, 0.25,), 
                (y_limits[0], y_limits[1]/5,), 
                lw=2,
                color="dimgrey", 
                ls=":",
                # alpha=0.7,
            )
            # plt.plot(
            #     # (1e-2, 2e-1,), 
            #     (1.5e-1, 5e1,), 
            #     (eqns.nus_eqn_W.rhs.subs({W: 0.25}), 
            #      eqns.nus_eqn_W.rhs.subs({W: 0.25}),),
            #     lw=2,
            #     color="dimgrey", 
            #     ls=":",
            #     alpha=0.7,
            # )
            # transition line higher W
            plt.plot(
                (2.5, 2.5,), 
                (y_limits[0], y_limits[1]/2.9,), 
                lw=2,
                color="dimgrey", 
                ls=":",
                # alpha=0.7,
            )
            # plt.plot(
            #     (1.5e0, 5e1,), 
            #     (eqns.nus_eqn_W.rhs.subs({W: 2.5}), 
            #      eqns.nus_eqn_W.rhs.subs({W: 2.5}),),
            #     lw=2,
            #     color="dimgrey", 
            #     ls=":",
            #     alpha=0.7,
            # )
            plt.text(
                0.18, 
                0.23, 
                r"low ${\mathcal{W}}$", 
                color="brown",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            plt.text(
                0.18, 
                0.13, 
                r"$\nu_\mathsf{{s}} \approx 1 + {\mathcal{W}}$", 
                color="brown",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            plt.text(
                0.51, 
                0.39, 
                r"transitional", 
                color="dimgrey",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            plt.text(
                0.51, 
                0.3, 
                r"${\mathcal{W}}$", 
                color="dimgrey",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            plt.text(
                0.78, 
                0.9, 
                r"high ${\mathcal{W}}$", 
                color="blue",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            plt.text(
                0.78, 
                0.8, 
                r"$\nu_\mathsf{{s}} \approx \dfrac{1}{2}+\sqrt{{\mathcal{W}}}$", 
                color="blue",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            
        if numerical_models is not None:
            for (idx, nus_soln,) in enumerate(numerical_models):
                plt.plot(
                    nus_soln.W,
                    nus_soln.ν_s,
                    "o",
                    c="k",
                    label=("numerical" if idx==0 else None),
                )
        
        plt.plot(
            W_array[W_array>1.5],
            0.5+np.sqrt(W_array[W_array>1.5]), 
            label=r"high ${\mathcal{W}}$ approx", 
            ls="--",
            c="blue",
        )
        plt.plot(
            W_array[W_array<0.45],
            1+W_array[W_array<0.45], 
            label=r"low ${\mathcal{W}}$ approx",  
            ls="--",
            c="brown",
        )
        if text_label is not None:
            plt.text(
                0.82,
                0.25, 
                text_label, 
                color="k", 
                size=14,
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
        axes.set_ylim(0.98, 7.5,)
        plt.legend(loc="upper left", fontsize=10,)
        plt.xlabel(r"Weathering number  ${\mathcal{W}}$  [-]")
        axes.set_ylabel(
            r"Dimensionless erosion rate  ${\nu}_\mathsf{{s}}$  [-]"
        )
        plt.grid(which="major", lw=1, ls="-", alpha=0.3,)
        plt.grid(which="minor", lw=1, ls="-", alpha=0.3,)

        alt_axes = axes.twinx()  
        alt_axes.autoscale(enable=True, tight=True)
        if do_loglog:
            alt_axes.set_yscale("log",)
            alt_axes.yaxis.set_major_formatter(ticker.FormatStrFormatter(""))
            alt_axes.yaxis.set_minor_formatter(ticker.FormatStrFormatter(""))
            alt_axes.set_ylim(0.98, 7.5,)

    def erosion_rate_steadystate_W_transition(
            self, 
            name: str,
            title: str|None=None,
            eqns: Equations|None=None, 
            text_label: Sequence|None=None,
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        """
        Plot the steady-state erosion rate $\\nu_s$ relative to its asymptotic 
        behavior.
        
        Visualize the transitional behavior of the steady-state erosion rate 
        $\\nu_s$ at intermediate weathering numbers $W$ by plotting how 
        asymptotes at low and high $W$ deviate from the full model.

        Args:
            name: key for figure dictionary
            title: optional title for figure
            eqns: 
                instance of 1d weathering-mediated erosion theory class
            text_label: 
                text annotation as tuple of form 
                (x-y coordinate, string, font size)
            fig_size: tuple (x,y) in inches (!)
        """
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        n_W_pts = 200
        W_array = np.exp(np.linspace(np.log(0.015), np.log(50), n_W_pts,))
        nus_array = np.array([
            eqns.nus_eqn_W.rhs.subs({W: W_}) for W_ in W_array
        ])
        plt.plot(
            W_array, 1+0*W_array, color="k", lw=1.5,
        )
        axes = plt.gca()
        plt.plot(
            W_array[W_array>1], 
            nus_array[W_array>1]/(0.5+np.sqrt(W_array[W_array>1])),
            color="b", 
            ls="-", 
            lw=1.5, 
            label=r"high ${\mathcal{W}}$ approx",
        )
        plt.plot(
            W_array[W_array<0.48], 
            nus_array[W_array<0.48]/(1+(W_array[W_array<0.48])),
            color="brown", 
            ls="-", 
            lw=1.5, 
            label=r"low ${\mathcal{W}}$ approx",
        )
        plt.text(
            0.18,
            0.18, 
            r"low ${\mathcal{W}}$", 
            color="brown",
            verticalalignment="center", 
            horizontalalignment="center",
            transform=axes.transAxes,
        )
        plt.text(
            0.18,
            0.26, 
            r"$\nu_\mathsf{{s}} \approx 1 + {\mathcal{W}}$", 
            color="brown",
            verticalalignment="center", 
            horizontalalignment="center",
            transform=axes.transAxes,
        )
        plt.text(
            0.5,
            0.53, 
            r"transitional ${\mathcal{W}}$", 
            color="gray",
            verticalalignment="center", 
            horizontalalignment="center",
            transform=axes.transAxes,
        )
        plt.text(
            0.82,
            0.85, 
            r"high ${\mathcal{W}}$", 
            color="blue",
            verticalalignment="center", 
            horizontalalignment="center",
            transform=axes.transAxes,
        )
        plt.text(
            0.82,
            0.75, 
            r"$\nu_\mathsf{{s}} \approx \dfrac{1}{2}+\sqrt{{\mathcal{W}}}$", 
            color="blue",
            verticalalignment="center", 
            horizontalalignment="center",
            transform=axes.transAxes,
        )
        if text_label is not None:
            plt.text(
                0.82,
                0.25, 
                text_label, 
                color="k", 
                size=14,
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
        plt.autoscale(enable=True, tight=True,)
        axes.set_xscale("log",) #nonposx="clip")
        y_limits = (0.9, 1.1,)
        plt.ylim(*y_limits)
        axes.set_yticks(np.linspace(*y_limits, 5,))
        axes.yaxis.set_major_formatter(ticker.FormatStrFormatter("%0.2f"))
        plt.plot((0.25, 0.25,), y_limits, color="gray", ls=":",)
        plt.plot((2.63, 2.63,), y_limits, color="gray", ls=":",)
        plt.legend(loc="upper left")
        plt.xlabel(r"Weathering number  ${\mathcal{W}}$  [-]")
        plt.ylabel(
            r"Approx erosion rate deviation  "
            + r"$\nu_\mathsf{{s}}^\mathrm{apx}/\nu_\mathsf{{s}}$  [-]"
        )
        # plt.grid(ls=":")

    def ref_erosion_rate_vs_ref_weathering_rate(
            self, 
            name: str,
            title: str|None=None,
            eqns: Equations|None=None, 
            k: float=1, 
            fig_size: tuple[float,float]=(6, 4,),
        ) -> None:
        """
        Plot baseline erosion rate versus baseline weathering rate.

        Args:
            name: key for figure dictionary
            title: optional title for figure
            eqns: 
                instance of 1d weathering-mediated erosion theory class
            k: weathering-depth model $k$ value
            fig_size: tuple (x,y) in inches (!)
        """
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        w0_array = 10**np.linspace(-3,+1,100)
        for vs in (0.1, 1, 3):
            v0_array = np.array([
                sy.N(eqns.v0_eqn_vs_w0.rhs.subs({v_s:vs, w_0:w0, k:k})) 
                for w0 in w0_array
            ])
            plt.plot(
                v0_array,
                w0_array, 
                label="$v_s=${}".format(vs),
            )
        plt.ylabel(r"$w_0$")
        plt.xlabel(r"${u_0}$")
        plt.legend(loc="lower right")
        plt.grid(ls=":")

    def ref_erosion_rate_vs_ref_weathering_rate_steadystate(
            self, 
            name: str,
            title: str|None=None,
            eqns: Equations|None=None, 
            fig_size: tuple[float,float]=(6, 4,),
        ) -> None:
        """
        Plot baseline erosion rate versus surface weakness at steady-state.

        Args:
            name: key for figure dictionary
            title: optional title for figure
            eqns: 
                instance of 1d weathering-mediated erosion theory class
            fig_size: tuple (x,y) in inches (!)
        """
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        etas0_array = 10**np.linspace(-2, +1, 100,)
        for vs in (0.5, 1, 1.5,):
            v0_array = np.array([
                sy.N(eqns.v0_eqn_etas0_vs.rhs.subs({v_s:vs,eta_s0:etas0_})) 
                for etas0_ in etas0_array
            ])
            plt.plot(
                etas0_array,
                v0_array, 
                label="$v_s=${}".format(vs),
            )
        plt.xlabel(r"Surface weakness (degree of weathering)  ${\omega}_{s0}$")
        plt.ylabel(r"Baseline (potential) erosion rate  ${u_0}$")
        plt.gca().invert_yaxis()
        plt.xlim(0,3)
        plt.ylim(5,0)
        plt.legend(loc="lower right")
        plt.grid(ls=":")
