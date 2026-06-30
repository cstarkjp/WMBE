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
            
            χ_domain_size (float):  
                length of chi solution domain (extracted from model_parameters)
            Δχ (float):        
                spacing between discrete chi solution points (extracted from model_parameters)
            n_χ_domain (int):       
                number of solution points in distance chi (extracted from model_parameters)
            τ_domain_size (float):  
                maximum duration of solution (truncated if/when front 
                exits chi domain) (extracted from model_parameters)
            τ_n_steps (int):        
                number of solution points in time $\\tau$
            Δτ (float):        
                spacing between discrete tau solution points
    
            χ_array (NDArray): 
                discrete distances $\\chi_i$ 
            τ_array (NDArray): 
                discrete times $\\tau^j$
            ω_array (NDArray): 
                discretized weakness profile $\\omega_i^j$
            φ_array (NDArray): 
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
            u_s (float): 
                predicted (by analytical solution) steady-state erosion rate $v_s$    
        """
        self.physical_parameters: dict = physical_parameters
        self.W: float = (
            physical_parameters[w_0]
                /
            (physical_parameters[k]*physical_parameters[u_0])
        )
        self.ν_s: float  = erosionrate_steadystate_W(self.W)
        self.u_s: float = self.ν_s*physical_parameters[u_0]
        self.physical_parameters.update({
            W: self.W, 
            ν_s: self.ν_s, 
            u_s: self.u_s,
        })
        
        self.model_parameters: dict = model_parameters
        self.χ_domain_size: float = model_parameters[x_domain_size]
        self.τ_domain_size: float = model_parameters[τ_domain_size]
        self.Δχ: float = model_parameters[Δχ]
        self.Δτ: float = model_parameters[Δτ]
        self.n_χ_domain: int = np.int64(self.χ_domain_size/self.Δχ)+1
        self.τ_n_steps: int  = np.int64(self.τ_domain_size/self.Δτ)+1

        self.ω_array: NDArray = np.zeros(
            (self.τ_n_steps, self.n_χ_domain,),
            dtype=np.float64,
        )
        self.ω_array[0] = np.ones(self.n_χ_domain,dtype=np.float64)
        self.φ_array: NDArray = np.zeros(
            self.τ_n_steps, 
            dtype=np.float64,
        )
        self.χ_array: NDArray = np.linspace(
            0,
            self.χ_domain_size,
            self.n_χ_domain,
            dtype=np.float64,
        )
        self.τ_array: NDArray = np.linspace(
            0,
            self.τ_domain_size,
            self.τ_n_steps,
            dtype=np.float64,
        )
        self.ν_array: NDArray  = np.zeros(
            self.τ_n_steps,
            dtype=np.float64,
        )
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
        ω = self.ω_array
        φ = self.φ_array
        χ = self.χ_array
        τ = self.τ_array
        Δχ = self.Δχ
        Δτ = self.Δτ
        j: int
        τ_step: float
        for (j, τ_step,) in enumerate(τ[:-1]):
            self.j = j
            f = np.int64(φ[j]/Δχ)
            f_right = φ[j]/Δχ-f
            f_left = 1.0-f_right
            fp1 = f+1
            fp2 = f+2
            if fp2>=self.n_χ_domain:
                break
            self.ν_array[j] = (f_left*ω[j,f]+f_right*ω[j,fp1])
            Delta_phi_j = (self.ν_array[j]*Δτ)/(W*2)
            φ[j+1] = φ[j]+Delta_phi_j
            ω[j+1,f:-1] = (
                  ω[j,f:-1] 
                + (Delta_phi_j*(ω[j,fp1:]-ω[j,f:-1]))/(Δχ)
                +  Δτ*self.neg_exp_Heaviside(χ[f:-1]-φ[j], Δχ,)
                )
            ω[j+1,-1] = (
                  ω[j,-1] 
                + (Delta_phi_j*(ω[j,-1]-ω[j,-2]))/(Δχ)
                +  Δτ*self.neg_exp_Heaviside(χ[-1]-φ[j], Δχ,)
            )
        self.ν_array[j+1] = (f_left*ω[j,f]+f_right*ω[j,fp1])
        di = self.ν_array.shape[0]//10
        self.nu_s_bar = np.mean(self.ν_array[4*di:6*di])
        self.physical_parameters.update({nu_s_bar:self.nu_s_bar})

    @staticmethod
    def neg_exp_Heaviside(
            χ: float|NDArray , 
            Δχ: float|NDArray,
        ) -> float|NDArray:
        """
        Negate, exponentiate and Heaviside clip.
        
        Assumes the argument χ = $\\chi$ is the dimensionless distance from the 
        erosion front, and Δχ = $\\Delta\\chi$  is the discrete spatial step size.
        Exponentiates $-\\chi H(\\chi+\\Delta\\chi)$ where H=Heaviside function.
        When invoked in computing $d{\\omega}/d{\\chi}$,  
        adding ${\\Delta}{\\chi}$ means
        a non-clipped value is returned for the sample point just to the left 
        (-ve $\\chi$). Thus the ${\\omega}$ gradient is estimated across the erosion 
        front, where sample points span the moving origin, as well as for all points
        $\\chi\\geq 0$.

        Args:
            χ: distance $\\chi$ from the erosion front
            Δχ: discrete spacing $\\Delta\\chi$
                between sample points along $\\chi$
        
        Returns:
            $\\exp(-\\chi)$ if $\\chi+\\Delta\\chi \\geq 0$, 1 otherwise
        """
        return np.exp(-χ*np.heaviside(χ+Δχ,0))