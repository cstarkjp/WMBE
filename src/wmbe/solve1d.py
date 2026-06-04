"""
---------------------------------------------------------------------

Module for performing finite-difference solution of weathering-driven
weakening of a 1d bedrock surface half-space and its concomitant erosion.

The moving-boundary problem is solved using a continuous-valued 
distance (from the erosion front) function that tracks motion with sub-grid
resolution. Upwind differencing and explicit Euler methods are employed. 

---------------------------------------------------------------------

Requires Python packages/modules:
  -  :mod:`numpy`

Imports symbols from :mod:`.symbols` module

---------------------------------------------------------------------

"""


import numpy as np
from .symbols import *


def negExpH(chi,dchi):
    """
    Negate, exponentiate and Heaviside clip.
    
    Assumes the argument chi= :math:`\\chi`
    is the dimensionless distance from the erosion front, 
    and dchi= :math:`\\Delta\\chi`  is the discrete spatial step size.
    Exponentiates :math:`-\\chi H(\\chi+\\Delta\\chi)` where H=Heaviside function.
    When invoked in computing :math:`d{\\eta}/d{\\chi}`, 
    adding :math:`\\Delta\\chi` means a non-clipped
    value is returned for the sample point just to the left (-ve :math:`\\chi`).
    Thus the :math:`\\eta` gradient is estimated across the erosion front, 
    where sample points
    span the moving origin, as well as for all points :math:`\\chi\\geq 0`.
    
    Args:
        chi (float): distance :math:`\\chi` from the erosion front
        dchi (float): discrete spacing :math:`\\Delta\\chi` 
            between sample points along :math:`\\chi`
    
    Returns:
        float:
            :math:`\\exp(-\\chi)` if :math:`\\chi+\\Delta\\chi \\geq 0`, 1 otherwise
    """
    return np.exp(-chi*np.heaviside(chi+dchi,0))

def nu_s_W(W):
    """
    Dimensionless steady-state speed of the erosion front :math:`\\nu_s(W)`.
    
    Assumes: :math:`\\nu_s = \\tfrac{1}{2}\\left(1+\sqrt{1+4W}\\right)`
    
    Args:
        W (float): weathering number :math:`W`
    
    Returns:
        float: dimensionless erosion rate :math:`\\nu_s`
    """
    return 0.5*(1+np.sqrt(1+4*W))

def eta_chi_tau(chi,tau,W):
    """
    Weathering-driven weakness function.
    
    Analytic solution for :math:`\\eta(\\chi,\\tau)` 
    as a function of dimensionless distance
    (depth into the rock) :math:`\\chi` and time :math:`\\tau`, and parameterized
    by weathering number :math:`W`, assuming an exponential-decay model for weathering.
    
    Args:
        tau (float): dimensionless time :math:`\\tau`
        chi (float): dimensionless distance :math:`\\chi`
        W (float):   weathering number :math:`W`
    
    Returns:
        float: 
            weakness :math:`\\eta(\\chi,\\tau;W)`
    """
    return (1+(W/nu_s_W(W))*np.exp(-(chi)))*np.heaviside(chi,0)


class ErosionWeathering:
    """
    Numerical solution of :math:`\\eta(\\chi,\\tau)` 
    and :math:`\\varphi(\\tau)` evolution 
    
    Class that provides a finite-difference method for solving the 
    :math:`(chi,tau)` evolution
    of a weakness profile :math:`\\eta(\\chi,\\tau)`
    and its eroding surface position :math:`\\varphi(\\tau)`
    as 2d array :math:`\\eta_i^j` and 1d vector :math:`\\varphi^j` respectively.
    and that provides dictionaries for the model and its numerical solution
    parameters.
    
    
    Args:
        pdict (dict): model parameters dictionary
        ndict (dict): numerical method parameters dictionary
        

    """
    def __init__(self, pdict, ndict):
        """
        Initialize class instance.
        
        
        Attributes:
            pdict (:obj:`dict`) : 
                model parameters dictionary, extended during & after instantiation
            ndict (:obj:`dict`) : 
                numerical method parameters dictionary
            
            chi_domain_size (:obj:`float`):  
                length of chi solution domain (extracted from ndict)
            Delta_chi (:obj:`float`):        
                spacing between discrete chi solution points (extracted from ndict)
            n_chi_domain (:obj:`int`):       
                number of solution points in distance chi (extracted from ndict)
            tau_domain_size (:obj:`float`):  
                maximum duration of solution (truncated if/when front 
                exits chi domain) (extracted from ndict)
            tau_n_steps (:obj:`int`):        
                number of solution points in time tau
            Delta_tau (:obj:`float`):        
                spacing between discrete tau solution points
    
            chi_array (:class:`numpy.ndarray`) : 
                discrete distances :math:`\\chi_i` 
            tau_array (:class:`numpy.ndarray`) : 
                discrete times :math:`\\tau^j`
            eta_array (:class:`numpy.ndarray`) : 
                discretized weakness profile :math:`\\eta_i^j`
            phi_array (:class:`numpy.ndarray`) : 
                discrete (in time) series of erosion front 
                positions :math:`\\phi^j` (smoothly resolved as floats)
            nu_array  (:class:`numpy.ndarray`) :  
                discrete (in time) series of 
                dimensionless erosion rates :math:`\\nu^j` 
            
            j (:obj:`int`) :  
                final time step index :math:`j`
    
            W (:obj:`float`)    : 
                weathering number :math:`W`
            nu_s (:obj:`float`) : 
                predicted (by analytical solution) dimensionless 
                steady-state erosion rate :math:`\\nu_s`
            v_s (:obj:`float`)  : 
                predicted (by analytical solution) steady-state erosion rate :math:`v_s`     
        """
        self.pdict = pdict
        self.W     = pdict[w_0]/(pdict[k]*pdict[v_0])
        self.nu_s  = nu_s_W(self.W)
        self.v_s   = self.nu_s*pdict[v_0]
        self.pdict.update({W:self.W, nu_s:self.nu_s, v_s:self.v_s})
        
        self.ndict = ndict
        self.chi_domain_size = ndict[chi_domain_size]
        self.tau_domain_size = ndict[tau_domain_size]
        self.Delta_chi = ndict[Delta_chi]
        self.Delta_tau = ndict[Delta_tau]
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
        self.j         = 0
        

    def solve(self):
        """
        Use an explicit finite-difference scheme to solve for evolution of
        a weakness profile :math:`\\eta_i^j` 
        and its eroding surface position :math:`\\phi^j`.
        
        Attributes:
            pdict (:obj:`dict`) : 
                model parameters dictionary, extended during & after instantiation
            nu_array  (:class:`numpy.ndarray`) :  discrete (in time) series of 
                                        dimensionless erosion rates :math:`\\nu^j` 
            
            j (int) :  final time step index :math:`j`
            nu_s_bar (:obj:`float`) : 
                post-hoc estimate (from averaging portion of solutions) of
                               dimensionless steady-state erosion rate 
                               :math:`\\overline{\\nu}_s`
        """
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
                +  Delta_tau*negExpH(chi[f:-1]-phi[j],Delta_chi)
                )
            eta[j+1,-1] = (
                  eta[j,-1] 
                + (Delta_phi_j*(eta[j,-1]-eta[j,-2]))/(Delta_chi)
                +  Delta_tau*negExpH(chi[-1]-phi[j],Delta_chi)
            )
        self.nu_array[j+1] = (f_left*eta[j,f]+f_right*eta[j,fp1])
        di = self.nu_array.shape[0]//10
        self.nu_s_bar = np.mean(self.nu_array[4*di:6*di])
        self.pdict.update({nu_s_bar:self.nu_s_bar})
