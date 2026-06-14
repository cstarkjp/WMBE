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

from wmbe.solve1d import eta_chi_tau, NumericalModel
from wmbe.symbols import *
from wmbe.theory import Equations
from wmbe.viz_base import VizBase

warnings.filterwarnings("ignore")

__all__ = [
    "VizSim",
]

class VizSimulations(VizBase):
    """
    Numerical simulation visualization class.
    """
    def frontspeed_evolution(
            self, 
            name: str,
            title: str|None=None,
            nm: NumericalModel|None=None, 
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        """
        Plot time-evolution of dimensionless erosion rate $\\nu(\\tau)$.
        
        Generate graph of numerical solutions $j$ of dimensionless rock-surface 
        erosion rate $\\nu_i^j$ over dimensionless time $\\tau^j$. 
        These solutions are provided in the class instance 
        :class:`~.solve1d.ErosionWeathering`.
        Legend-label by weathering number for this instance.
        """            
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        #     ew (:class:`~.solve1d.ErosionWeathering`): 
        #             instance of 1d weathering-mediated erosion model :mod:`~.solve1d` class
        _ = self.create_figure(name=name, size=fig_size,)
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
        """
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        #     ew (:class:`~.solve1d.ErosionWeathering`): 
        #             instance of 1d weathering-mediated erosion model :mod:`~.solve1d` class
        #     tc (int): 
        #             $\\tau^j`  slicing "rate"
        #     nd (int): 
        #             number of decimal places in $\\tau`  legend label
        #     text_label (list): text annotation as list of form (x-y coordinate, string, 
        #                     font size)
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        chi_ = nm.chi_array
        tau_ = nm.tau_array
        eta_ = nm.eta_array
        j_ = nm.j

        tau_slices1 = np.linspace(
            0, 
            (nm.tau_n_steps-1)//tc,
            num=5,
            endpoint=True,
            dtype=np.int64,
        )
        tau_slices2 = np.linspace(
            (nm.tau_n_steps-1)//tc*2,
            j_+1,
            num=5,
            endpoint=True,
            dtype=np.int64,
        )
        tau_slices = np.concatenate((tau_slices1, tau_slices2,))
        cmap = plt.cm.brg.reversed()
        label=r"$\tau$={0:3."+str(nd)+"f}"
        for idx,tau_slice in enumerate(tau_slices):
            eta_slice = eta_[tau_slice]
            chi_front = chi_[eta_slice==0]
            chi_front = 0 if chi_front.shape[0]==0 else chi_front[-1]
            plt.plot(
                chi_[chi_>=chi_front],
                eta_slice[chi_>=chi_front],
                color=cmap(idx/tau_slices.size),
                label=label.format(tau_[tau_slice]),
            )
            
        axes = plt.gca()
        plt.xlim(chi_[0],chi_[-1])
        plt.xlim(-(chi_[-1]-chi_[0])/30,chi_[-1]*1.08)
        x_limits = plt.xlim()
        y_limits = plt.ylim()
        plt.ylim(y_limits[0]/3,y_limits[1])
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
        """
        # Args:
        #     tau (numpy.ndarray): time slices $\\tau^j` of numerical solution 
        #     nu (numpy.ndarray): speed of erosion front $\\nu^j` 
        #                         at each time slice $\\tau^j` 
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})

        plt.plot(tau, nu, "o-",)
        plt.ylim(1.204, 1.21,);
        plt.xlim(4, 4.03,);

    def weakness_steadystate(
            self, 
            name: str,
            title: str|None=None,
            nm: NumericalModel|None=None, 
            fig_size: tuple[float,float]=(6,4,),
        ) -> None:
        """
        Plot steady-state solution of weakness ${\\omega}_s$.
        
        Graph the numerical solution of the 1d weathering-mediated erosion model
        for weakness ${\\omega}_s(\\chi_i | {\\mathcal{W}})$ as a function of 
        depth from the rock surface $\\chi_i$ for a given value of the 
        weathering number ${\\mathcal{W}}$.
        """
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>`
        #                                     figure 
        #     ew (:class:`~.solve1d.ErosionWeathering`): 
        #             instance of 1d weathering-mediated erosion model :mod:`~.solve1d` class
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        j_ = (nm.j*3)//4
        i_offset = nm.n_chi_domain//7
        phi_  = nm.phi_array[j_]
        phi0_ = int(phi_/nm.Delta_chi)-i_offset
        chi_s_= nm.chi_array[phi0_:]-nm.chi_array[phi0_+i_offset]
        tau_  = nm.tau_array[-1]

        eta_s_numerical  = nm.eta_array[j_, phi0_:]
        eta_s_analytical = eta_chi_tau(chi_s_, tau_, nm.W)
        
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
        plt.xlim(nm.tau_array[0], nm.tau_array[-1],)
        
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
        """
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>`
        #                                     figure 
        #     ew_list (list): 
        #             list of instances of 1d weathering-mediated erosion model 
        #             :mod:`~.solve1d` class :class:`~.solve1d.ErosionWeathering`
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})

        cmap = plt.cm.brg
        
        chi_min = 0
        for (idx, (nm, label,),) in enumerate(numerical_models):
            j_ = (nm.j*3)//4
            i_offset = nm.n_chi_domain//10
            phi_  = nm.phi_array[j_]
            phi0_ = int(phi_/nm.Delta_chi)-i_offset
            chi_s_= nm.chi_array[phi0_:]-nm.chi_array[phi0_+i_offset]
            # tau_  = nm.tau_array[-1]
            chi_min = min(chi_min,np.min(chi_s_))
            color=cmap(idx/len(numerical_models))
            eta_s_numerical  = nm.eta_array[j_,phi0_:]
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
            
        plt.xlim(chi_min, chi_max,)
        plt.legend(fontsize=11,)
        plt.xlabel(r"Distance relative to front  $\chi_s=\chi-\varphi_s$  [-]")
        plt.ylabel(r"Weakness  ${\omega}_s(\chi_s)$  [-]")
        plt.grid(ls=":")

    def erosionrate_steadystate_W(
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
        """
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        #     em (:class:`~.theory.WeatheringMediatedErosion`): 
        #             instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
        #     do_loglog (bool): 
        #         flag whether to use log scales on both axes
        #     nus_solns_list (list): 
        #         set of numerical solutions of $\nu_s$
        #     text_label (list): text annotation as list of form (x-y coordinate, string, 
        #                     font size)
        _ = self.create_figure(name=name, size=fig_size,)
        if title is None:
            plt.title(title, fontdict={"fontsize": 11.5})
        
        n_W_pts = 200
        W_array   = np.exp(np.linspace(np.log(0.01), np.log(50), n_W_pts,))
        nus_array = np.array([
            eqns.nus_eqn_W.rhs.subs({W: W_}) 
            for W_ in W_array
        ])
        y_limits = (nus_array[0]*0.95, nus_array[-1],)
    
        plt.plot(W_array, nus_array, color="k", lw=1.5, label="analytical",)
        plt.autoscale(enable=True, tight=True)
        axes = plt.gca()
        if do_loglog:
            axes.set_xscale("log",)
            axes.set_yscale("log",)
            axes.yaxis.set_major_formatter(ticker.FormatStrFormatter("%0.0f"))
            axes.yaxis.set_minor_formatter(ticker.FormatStrFormatter("%0.0f"))
            plt.plot(
                (0.25, 0.25,), 
                (y_limits[0], y_limits[1]/3,), 
                color="gray", 
                ls=":",
            )
            plt.plot(
                (2.63, 2.63,), 
                (y_limits[0], y_limits[1],), 
                color="gray", 
                ls=":",
            )
            plt.text(
                0.2, 
                0.23, 
                r"low ${\mathcal{W}}$", 
                color="brown",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            plt.text(
                0.2, 
                0.13, 
                r"$\nu_\mathsf{{s}} \approx 1 + {\mathcal{W}}$", 
                color="brown",
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
            plt.text(
                0.52, 
                0.43, 
                r"transitional ${\mathcal{W}}$", 
                color="gray",
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
                0.82,
                0.25, 
                text_label, 
                color="k", 
                size=14,
                verticalalignment="center", 
                horizontalalignment="center",
                transform=axes.transAxes,
            )
        
        plt.legend(loc="upper left", fontsize=10,)
        plt.xlabel(r"Weathering number  ${\mathcal{W}}$  [-]")
        plt.ylabel(r"Erosion rate  ${\nu}_\mathsf{{s}}$  [-]")
        # plt.grid(ls=":")

    def erosionrate_steadystate_W_transition(
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
        """
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        #     em (:class:`~.theory.WeatheringMediatedErosion`): 
        #             instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
        #     text_label (list): text annotation as list of form (x-y coordinate, string, 
        #                     font size)
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
            W_array[W_array<0.48], 
            nus_array[W_array<0.48]/(1+(W_array[W_array<0.48])),
            color="brown", 
            ls="-", 
            lw=1.5, 
            label=r"low ${\mathcal{W}}$ approx",
        )
        plt.plot(
            W_array[W_array>1], 
            nus_array[W_array>1]/(0.5+np.sqrt(W_array[W_array>1])),
            color="b", 
            ls="-", 
            lw=1.5, 
            label=r"high ${\mathcal{W}}$ approx",
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

    def referosionrate_vs_refweatheringrate(
            self, 
            name: str,
            title: str|None=None,
            eqns: Equations|None=None, 
            k: float=1, 
            fig_size: tuple[float,float]=(6, 4,),
        ) -> None:
        """
        Plot baseline erosion rate versus baseline weathering rate.
        """
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        #     em (:class:`~.theory.WeatheringMediatedErosion`): 
        #             instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
        #     k_ (float): 
        #         weathering-weakening rate profile reciprocal e-folding depth
        #     text_label (list): text annotation as list of form (x-y coordinate, string, 
        #                     font size)
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

    def referosionrate_vs_refweatheringrate_steadystate(
            self, 
            name: str,
            title: str|None=None,
            eqns: Equations|None=None, 
            fig_size: tuple[float,float]=(6, 4,),
        ) -> None:
        """
        Plot baseline erosion rate versus surface weakness at steady-state.
        """
        # Args:
        #     fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
        #         reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        #     em (:class:`~.theory.WeatheringMediatedErosion`): 
        #             instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
        #     text_label (list): text annotation as list of form (x-y coordinate, string, 
        #                     font size)
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
        Plot numerical solutions applied to channel cross-section model 
        (vertical profiles).
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
            label=r"${\\mathcal{W}}=w_0/{u_0} k$", 
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
            r"Weathering number  ${\\mathcal{W}}(z)$", 
            color="forestgreen",
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