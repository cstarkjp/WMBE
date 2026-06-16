"""
Module for performing finite-difference solution of weathering-driven
weakening of a 1d bedrock surface half-space and its concomitant erosion.

The moving-boundary problem is solved using a continuous-valued 
distance (from the erosion front) function that tracks motion with sub-grid
resolution. Upwind differencing and explicit Euler methods are employed. 
"""
import warnings
import numpy as np
from numpy.typing import NDArray
from wmbe.theory.symbols import *

warnings.filterwarnings("ignore")

__all__ = [
    "erosionrate_steadystate_W",
    "weakness_chi_tau",
    "NumericalModel",
]

def weakness_chi_tau(chi: NDArray, tau: NDArray, W: float,) -> float|NDArray:
    """
    Weathering-driven weakness function.
    
    Analytic solution for $\\omega(\\chi,\\tau)$ as a function of dimensionless 
    distance (depth into the rock) $\\chi$ and time $\\tau$, and parameterized
    by weathering number $W$, assuming an exponential-decay model for 
    weathering.

    Args:
        chi: dimensionless distance $\\chi$
        tau: dimensionless time $\\tau$
        W:   weathering number $W$

    Returns:
        weakness $\\omega(\\chi,\\tau; W)$
    """    
    return (
        (1+(W/erosionrate_steadystate_W(W))
            *np.exp(-(chi)))
            *np.heaviside(chi,0)
    )

def erosionrate_steadystate_W(W: float) -> float|NDArray:
    """
    Dimensionless steady-state speed of the erosion front $\\nu_s(W)$.

    Assumes: $\\nu_s = \\tfrac{1}{2}\\left(1+\\sqrt{1+4W}\\right)$.
    
    Args:
        W: weathering number $W$
    
    Returns:
        dimensionless erosion rate $\\omega_s$
    """    
    return 0.5*(1+np.sqrt(1+4*W))


class NumericalModel:
    """
    Numerical solution of $\\omega(\\chi,\\tau)$ and $\\phi(\\tau)$ evolution.
    
    Class that provides a finite-difference method for solving the 
    space-time evolution of a weakness profile $\\omega(\\chi,\\tau)$ and 
    its eroding surface position $\\phi(\\tau)$ as 2d array $\\omega_i^j$ 
    and 1d vector $\\phi^j$ respectively, and that provides dictionaries 
    for the model and its numerical solution parameters.

    Args:
        physical_parameters: model parameters dictionary
        model_parameters: numerical method parameters dictionary

    """
    def __init__(
            self, 
            physical_parameters: dict, 
            model_parameters: dict,
        ) -> None:
        """
        Initialize class instance.
        
        Attributes:
            parameters (dict): 
                model parameters dictionary, extended during & after instantiation
            model_parameters (dict): 
                numerical method parameters dictionary
            
            chi_domain_size (float):  
                length of chi solution domain (extracted from model_parameters)
            Delta_chi (float):        
                spacing between discrete chi solution points (extracted from model_parameters)
            n_chi_domain (int):       
                number of solution points in distance chi (extracted from model_parameters)
            tau_domain_size (float):  
                maximum duration of solution (truncated if/when front 
                exits chi domain) (extracted from model_parameters)
            tau_n_steps (int):        
                number of solution points in time $\\tau$
            Delta_tau (float):        
                spacing between discrete tau solution points
    
            chi_array (NDArray): 
                discrete distances $\\chi_i$ 
            tau_array (NDArray): 
                discrete times $\\tau^j$
            eta_array (NDArray): 
                discretized weakness profile $\\omega_i^j$
            phi_array (NDArray): 
                discrete (in time) series of erosion front 
                positions $\\phi^j$ (smoothly resolved as floats)
            nu_array (NDArray):  
                discrete (in time) series of 
                dimensionless erosion rates $\\nu^j$
            
            j (int):  
                final time step index $j$
    
            W (float): 
                weathering number $W$
            nu_s (float): 
                predicted (by analytical solution) dimensionless 
                steady-state erosion rate $\\nu_s$
            v_s (float): 
                predicted (by analytical solution) steady-state erosion rate $v_s$    
        """
        self.physical_parameters = physical_parameters
        self.W = (
            physical_parameters[w_0]
                /
            (physical_parameters[k]*physical_parameters[v_0])
        )
        self.nu_s  = erosionrate_steadystate_W(self.W)
        self.v_s = self.nu_s*physical_parameters[v_0]
        self.physical_parameters.update({W:self.W, nu_s:self.nu_s, v_s:self.v_s})
        
        self.model_parameters = model_parameters
        self.chi_domain_size = model_parameters[chi_domain_size]
        self.tau_domain_size = model_parameters[tau_domain_size]
        self.Delta_chi = model_parameters[Delta_chi]
        self.Delta_tau = model_parameters[Delta_tau]
        self.n_chi_domain = np.int64(self.chi_domain_size/self.Delta_chi)+1
        self.tau_n_steps  = np.int64(self.tau_domain_size/self.Delta_tau)+1

        self.eta_array = np.zeros((self.tau_n_steps,self.n_chi_domain),dtype=np.float64)
        self.eta_array[0] = np.ones(self.n_chi_domain,dtype=np.float64)
        self.phi_array = np.zeros((self.tau_n_steps),dtype=np.float64)
        self.chi_array = np.linspace(0,self.chi_domain_size,self.n_chi_domain,
                                     dtype=np.float64)
        self.tau_array = np.linspace(0,self.tau_domain_size,self.tau_n_steps,
                                     dtype=np.float64)
        self.nu_array  = np.zeros((self.tau_n_steps),dtype=np.float64)
        self.j = 0
        
    def solve(self) -> None:
        """
        Use an explicit finite-difference scheme to solve for evolution of a 
        weakness profile $\\omega_i^j$ and its eroding surface position 
        $\\phi^j$.
        """
        # Attributes:
        #     parameters (dict) : 
        #         model parameters dictionary, extended during & after instantiation
        #     nu_array  (NDArray) :  discrete (in time) series of 
        #                                 dimensionless erosion rates $\\nu^j` 
            
        #     j (int) :  final time step index $j`
        #     nu_s_bar (float) : 
        #         post-hoc estimate (from averaging portion of solutions) of
        #                        dimensionless steady-state erosion rate 
        #                        $\\overline{\\nu}_s`
        W   = self.W
        eta = self.eta_array
        phi = self.phi_array
        chi = self.chi_array
        tau = self.tau_array
        Delta_chi = self.Delta_chi
        Delta_tau = self.Delta_tau
        for j,tau_step in enumerate(tau[:-1]):
            self.j = j
            f = np.int64(phi[j]/Delta_chi)
            f_right = phi[j]/Delta_chi-f
            f_left = 1.0-f_right
            fp1 = f+1
            fp2 = f+2
            if fp2>=self.n_chi_domain:
                break
            self.nu_array[j] = (f_left*eta[j,f]+f_right*eta[j,fp1])
            Delta_phi_j = (self.nu_array[j]*Delta_tau)/(W*2)
            phi[j+1] = phi[j]+Delta_phi_j
            eta[j+1,f:-1] = (
                  eta[j,f:-1] 
                + (Delta_phi_j*(eta[j,fp1:]-eta[j,f:-1]))/(Delta_chi)
                +  Delta_tau*self.neg_exp_Heaviside(chi[f:-1]-phi[j],Delta_chi)
                )
            eta[j+1,-1] = (
                  eta[j,-1] 
                + (Delta_phi_j*(eta[j,-1]-eta[j,-2]))/(Delta_chi)
                +  Delta_tau*self.neg_exp_Heaviside(chi[-1]-phi[j],Delta_chi)
            )
        self.nu_array[j+1] = (f_left*eta[j,f]+f_right*eta[j,fp1])
        di = self.nu_array.shape[0]//10
        self.nu_s_bar = np.mean(self.nu_array[4*di:6*di])
        self.physical_parameters.update({nu_s_bar:self.nu_s_bar})

    @staticmethod
    def neg_exp_Heaviside(
        chi: float|NDArray , dchi: float|NDArray,
    ) -> float|NDArray:
        """
        Negate, exponentiate and Heaviside clip.
        
        Assumes the argument chi= $\\chi$ is the dimensionless distance from the 
        erosion front, and dchi= $\\Delta\\chi$  is the discrete spatial step size.
        Exponentiates $-\\chi H(\\chi+\\Delta\\chi)$ where H=Heaviside function.
        When invoked in computing $d{\\omega}/d{\\chi}$,  adding ${\\Delta}{\\chi}$ means
        a non-clipped value is returned for the sample point just to the left 
        (-ve $\\chi$). Thus the ${\\omega}$ gradient is estimated across the erosion 
        front, where sample points span the moving origin, as well as for all points
        $\\chi\\geq 0$.
        """
        # Args:
        #     chi (float): distance $\chi$ from the erosion front
        #     dchi (float): discrete spacing $\Delta\chi$
        #         between sample points along $\chi$
        
        # Returns:
        #     float:
        #         $\exp(-\chi)$ if $\chi+\Delta\chi \geq 0$, 1 otherwise
        return np.exp(-chi*np.heaviside(chi+dchi,0))