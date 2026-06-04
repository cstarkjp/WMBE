"""
---------------------------------------------------------------------

Module applying the weathering-mediated erosion model to bedrock channel geometry

---------------------------------------------------------------------

Requires Python packages/modules:
  -  :mod:`numpy`
  -  :mod:`scipy.integrate` (integrate)
  
Imports symbols from :mod:`.symbols` module

---------------------------------------------------------------------

"""


import numpy as np
from scipy import integrate
from .symbols import *



class ChannelWall:
    """
    Numerical solution of combined models of weathering-mediated erosion & bedrock channel
    
    Class that provides numerical solution of a combination model 1d weathering-mediated
    erosion and 1+1d bedrock channel cross-section (at the channel wall).
    
    Args:
        em (:class:`~.theory.WeatheringMediatedErosion`): 
                instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
        pdict (dict): model parameters dictionary
        
        
    """
    def __init__(self, em, pdict):
        """
        Initialize class instance.

        Attributes:
            pdict (:obj:`dict`) : 
                model parameters dictionary
            w0_eqn_z_calibrated (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                surface weakness calibrated using model parameters: 
                :math:`w_0 =  w_r H_s[z,z_\\mathrm{wc},k_\\mathrm{w}]`         
                where 
                :math:`H_s = \\dfrac{1}{2}\\left(1 + \\tanh{[\\kappa_w(z-z_\\mathrm{wc}]}\\right)`
            v0_eqn_z_calibrated (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                erosion rate calibrated using model parameters: 
                :math:`v_0 =  v_r \\left\{ (h-H_s[z,z_\\mathrm{vc},k_\\mathrm{v}])(1-v_b)+v_b \\right\}`
                where 
                :math:`H_s = \\dfrac{1}{2}\\left(1 + \\tanh{[\\kappa_v(z-z_\\mathrm{vc}]}\\right)`
            vs_eqn (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                surface-normal erosion rate (uncalibrated)
                :math:`v_0 =  v_r \\left\{ (h-H_s[z,z_\\mathrm{vc},k_\\mathrm{v}])(1-v_b)+v_b \\right\}`
                where 
                :math:`H_s = \\dfrac{1}{2}\\left(1 + \\tanh{[\\kappa_v(z-z_\\mathrm{vc}]}\\right)`
            W_eqn_z_calibrated  (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                weathering number calibrated using model parameters: 
                :math:`W = \\dfrac{w_0}{k v_0}`
            vs_eqn_z_calibrated (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                surface-normal erosion rate calibrated using model parameters: 
                :math:`v_0 =  v_r \\left\{ (h-H_s[z,z_\\mathrm{vc},k_\\mathrm{v}])(1-v_b)+v_b \\right\}`
                where 
                :math:`H_s = \\dfrac{1}{2}\\left(1 + \\tanh{[\\kappa_v(z-z_\\mathrm{vc}]}\\right)`
        """
        self.pdict = pdict
        self.w0_eqn_z_calibrated = em.w0_eqn_wr_z.subs(self.pdict)
        self.v0_eqn_z_calibrated = em.v0_eqn_vr_h_z.subs(self.pdict)
        self.W_eqn_z_calibrated  = em.W_eqn.subs(self.pdict)
        self.vs_eqn = em.vs_eqn_w0_v0.subs({v_0:em.v0_eqn_vr_h_z.rhs})
        self.vs_eqn_z_calibrated \
            = self.vs_eqn.subs({w_0:self.w0_eqn_z_calibrated.rhs}).subs(self.pdict)

    def compute_vertical_profiles(self, n_pts=100):
        """
        Compute dependence of various properties with height above channel base

        Args:
            n_pts (int): 
                number of sampling points along vertical profile  

        Attributes:
            n_pts (:obj:`int`): 
                record of number of sampling points
            z_array  (:class:`numpy.ndarray`) :
                sample points along vertical profile 
            vs_eqn_z_calibrated (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                weathering-mediated surface-normal erosion rate at steady state
                calibrated by model parameters supplied in :attr:`wme.pdict`
            w0_array (:class:`numpy.ndarray`) :
                surface weakness along vertical profile (at steady state)
            v0_array (:class:`numpy.ndarray`) :
                baseline (absent weathering-driven weakening) surface-normal erosion rate 
                along vertical profile  (at steady state)
            vs_array (:class:`numpy.ndarray`) :
                weathering-mediated surface-normal erosion rate along vertical profile 
                (at steady state)
            eta0_array (:class:`numpy.ndarray`) :
                surface weakness along vertical profile (at steady state)
            W_array (:class:`numpy.ndarray`) :
                weathering number along vertical profile

        """
        self.n_pts = n_pts
        self.z_array  = np.linspace(0,1,self.n_pts)
        self.w0_array = np.array([
            np.float64(self.w0_eqn_z_calibrated.rhs.subs(z,z__)) for z__ in self.z_array])
        self.v0_array = np.array([
            np.float64(self.v0_eqn_z_calibrated.rhs.subs(z,z__)) for z__ in self.z_array])
        self.vs_array = np.array([
            np.float64(self.vs_eqn_z_calibrated.rhs.subs(z,z__)) for z__ in self.z_array])
        self.eta0_array = self.vs_array/self.v0_array
        self.W_array = np.array([
            np.float64(  (self.W_eqn_z_calibrated.rhs.subs(z,z__))
                            .subs({v_0:self.v0_array[idx]})
                            .subs({w_0:self.w0_array[idx]})  )
                                for idx,z__ in enumerate(self.z_array)])
        self.vs_eqn_z_calibrated = self.vs_eqn.subs({w_0:self.w0_eqn_z_calibrated.rhs}) \
                                        .subs(self.pdict)
    
    
    def compute_cross_section(self):
        """
        Compute dependence of various properties with height above channel base
        
        Attributes:
            vs_array (:class:`numpy.ndarray`) :
                weathering-mediated surface-normal erosion rate along vertical profile 
                (at steady state)                
            dzdy_array (:class:`numpy.ndarray`) :
                weathering-mediated vertical erosion rate along vertical profile 
                (at steady state)
            y_array (:class:`numpy.ndarray`) :
                horizontal positions of sample points along vertica profile
            ch_y_array (:class:`numpy.ndarray`) :
                horizontal positions of sample points along channel boundary
            ch_z_array (:class:`numpy.ndarray`) :
                vertical positions of sample points along channel boundary

        """
        self.vs_array = np.array([
            np.float64(self.vs_eqn_z_calibrated.rhs.subs(z,z__)) for z__ in self.z_array])
        self.dzdy_array \
            = 1/np.sqrt(((np.max(self.vs_array)+0.01)/self.vs_array)**2-1)
        self.y_array \
            = integrate.cumtrapz(self.dzdy_array,x=self.z_array, initial=0)
        self.ch_y_array = np.concatenate(
            (np.linspace(-np.max(self.z_array)*0.3,0,self.z_array.shape[0]),self.y_array))
        self.ch_z_array \
            = np.concatenate(( np.zeros(self.y_array.shape[0]),self.z_array ))
    
        
        
        
        