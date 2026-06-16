"""
Module applying the weathering-mediated erosion model to 1d bedrock channel.
"""
import warnings
import numpy as np
from scipy import integrate

from wmbe.theory.symbols import *
from wmbe.theory.equations import Equations

warnings.filterwarnings("ignore")

__all__ = [
    "ChannelWall",
]

class ChannelWall:
    """
    Numerical solution of combined models of weathering-mediated erosion 
    & bedrock channel.
    
    Class that provides numerical solution of a combination model 1d 
    weathering-mediated erosion and 1+1d bedrock channel cross-section 
    (at the channel wall).

    Args:
        eqns: instance of 1d weathering-mediated erosion theory class
        physical_parameters: physical parameters dictionary
    """
    def __init__(self, eqns: Equations, physical_parameters: dict) -> None:
        """
        Initialize class instance.

        Attributes:
            physical_parameters (dict): physical parameters dictionary
            w0_eqn_z_calibrated (sympy.Eq):
                surface weakness calibrated using model parameters: 
                $w_0 =  w_r H_s[z,z_\\mathrm{wc},k_\\mathrm{w}]$
                where 
                $H_s = \\dfrac{1}{2}\\left(1 + \\tanh{[\\kappa_w(z-z_\\mathrm{wc}]}\\right)$
            v0_eqn_z_calibrated (sympy.Eq):
                erosion rate calibrated using model parameters: 
                $v_0 =  v_r \\left\\{ (h-H_s[z,z_\\mathrm{vc},k_\\mathrm{v}])(1-v_b)+v_b \\right\\}$
                where 
                $H_s = \\dfrac{1}{2}\\left(1 + \\tanh{[\\kappa_v(z-z_\\mathrm{vc}]}\\right)$
            vs_eqn (sympy.Eq):
                surface-normal erosion rate (uncalibrated)
                $v_0 =  v_r \\left\\{ (h-H_s[z,z_\\mathrm{vc},k_\\mathrm{v}])(1-v_b)+v_b \\right\\}$
                where 
                $H_s = \\dfrac{1}{2}\\left(1 + \\tanh{[\\kappa_v(z-z_\\mathrm{vc}]}\\right)$
            W_eqn_z_calibrated (sympy.Eq):
                weathering number calibrated using model parameters: 
                $W = \\dfrac{w_0}{k v_0}$
            vs_eqn_z_calibrated (sympy.Eq):
                surface-normal erosion rate calibrated using model parameters: 
                $v_0 =  v_r \\left\\{ (h-H_s[z,z_\\mathrm{vc},k_\\mathrm{v}])(1-v_b)+v_b \\right\\}$
                where 
                $H_s = \\dfrac{1}{2}\\left(1 + \\tanh{[\\kappa_v(z-z_\\mathrm{vc}]}\\right)$
        """
        self.physical_parameters = physical_parameters
        self.w0_eqn_z_calibrated \
            = eqns.w0_eqn_wr_z.subs(self.physical_parameters)
        self.v0_eqn_z_calibrated \
            = eqns.v0_eqn_vr_h_z.subs(self.physical_parameters)
        self.W_eqn_z_calibrated \
            = eqns.W_eqn.subs(self.physical_parameters)
        self.vs_eqn \
            = eqns.vs_eqn_w0_v0.subs({v_0:eqns.v0_eqn_vr_h_z.rhs})
        self.vs_eqn_z_calibrated = (
            self.vs_eqn
                .subs({w_0:self.w0_eqn_z_calibrated.rhs})
                .subs(self.physical_parameters)
        )
    
    def compute_vertical_profiles(self, n_pts: int=100,) -> None:
        """
        Compute dependence of various properties with height above channel base.

        Args:
            n_pts: number of sampling points along vertical profile  

        Attributes:
            n_pts (int): 
                record of number of sampling points
            z_array (NDArray):
                sample points along vertical profile 
            vs_eqn_z_calibrated (sympy.Eq): 
                weathering-mediated surface-normal erosion rate at steady state
                calibrated by model parameters supplied in dict
            w0_array (NDArray):
                surface weakness along vertical profile (at steady state)
            v0_array (NDArray):
                baseline (absent weathering-driven weakening) surface-normal 
                erosion rate along vertical profile  (at steady state)
            vs_array (NDArray):
                weathering-mediated surface-normal erosion rate along vertical 
                profile (at steady state)
            eta0_array (NDArray):
                surface weakness along vertical profile (at steady state)
            W_array (NDArray):
                weathering number along vertical profile
        """
        self.n_pts = n_pts
        self.z_array = np.linspace(0,1,self.n_pts)
        self.w0_array = np.array([
            np.float64(self.w0_eqn_z_calibrated.rhs.subs(z,z__)) 
            for z__ in self.z_array
        ])
        self.v0_array = np.array([
            np.float64(self.v0_eqn_z_calibrated.rhs.subs(z,z__)) 
            for z__ in self.z_array
        ])
        self.vs_array = np.array([
            np.float64(self.vs_eqn_z_calibrated.rhs.subs(z,z__)) 
            for z__ in self.z_array
        ])
        self.eta0_array = self.vs_array / self.v0_array
        self.W_array = np.array([
            np.float64(  (self.W_eqn_z_calibrated.rhs.subs(z,z__))
                            .subs({v_0:self.v0_array[idx]})
                            .subs({w_0:self.w0_array[idx]})  )
            for idx,z__ in enumerate(self.z_array)
        ])
        self.vs_eqn_z_calibrated \
            = self.vs_eqn.subs({w_0:self.w0_eqn_z_calibrated.rhs}) \
                         .subs(self.physical_parameters)
    
    def compute_cross_section(self) -> None:
        """
        Compute dependence of various properties with height above channel base.

        Attributes:
            vs_array (NDArray):
                weathering-mediated surface-normal erosion rate along vertical 
                profile (at steady state)                
            dzdy_array (NDArray):
                weathering-mediated vertical erosion rate along vertical profile 
                (at steady state)
            y_array (NDArray):
                horizontal positions of sample points along vertica profile
            ch_y_array (NDArray):
                horizontal positions of sample points along channel boundary
            ch_z_array (NDArray):
                vertical positions of sample points along channel boundary
        """
        self.vs_array = np.array([
            np.float64(self.vs_eqn_z_calibrated.rhs.subs(z,z__)) 
            for z__ in self.z_array
        ])
        self.dzdy_array \
            = 1/np.sqrt(((np.max(self.vs_array)+0.01)/self.vs_array)**2-1)
        self.y_array = integrate.cumulative_trapezoid(
            self.dzdy_array,
            x=self.z_array, 
            initial=0,
        )
        self.ch_y_array = np.concatenate((
            np.linspace(-np.max(self.z_array)*0.3, 0, self.z_array.shape[0]),
            self.y_array,
        ))
        self.ch_z_array = np.concatenate((
            np.zeros(self.y_array.shape[0]),
            self.z_array 
        ))