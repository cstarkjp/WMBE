"""
---------------------------------------------------------------------

Module for :mod:`sympy` exposition of weathering-mediated erosion theory

---------------------------------------------------------------------

Requires Python packages/modules:
  -  :mod:`numpy`
  -  :mod:`sympy`

Imports symbols from :mod:`.symbols` module

---------------------------------------------------------------------

"""


import numpy as np
import sympy as sy
from sympy import Eq
from .symbols import *



class WeatheringMediatedErosion:
    """
    Weathering-mediated erosion theory in various forms
    
    """
    def __init__(self):
        """
        Initialize class instance.

        Attributes:
            W_eqn           (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                weathering number: 
                :math:`W = \\dfrac{w_0}{k v_0}`
            nus_eqn_W       (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                dimensionless steady-state erosion rate: 
                :math:`\\nu_s = \\dfrac{1}{2}(1+\\sqrt{1+4W})`
            nus_eqn_w0_v0   (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                dimensionless steady-state erosion rate:
                :math:`\\nu_s = \\frac{\\sqrt{1 + \\frac{4 w_{0}}{k v_{0}}}}{2} + \\frac{1}{2}`
            etas0_eqn_W     (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                steady-state surface weakness:
                :math:`\\eta_{s0} = \\dfrac{\\sqrt{4 W + 1}}{2} + \\frac{1}{2}`
            etas0_eqn_w0_v0 (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                steady-state surface weakness:
                :math:`\\eta_{s0} = \\dfrac{\\sqrt{1 + \\frac{4 w_{0}}{k v_{0}}}}{2} + \\frac{1}{2}`
            nus_eqn_etas0   (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                dimensionless steady-state erosion rate:
                :math:`\\nu_s =  \\eta_{s}(0)`
            vs_eqn_etas0_v0 (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                steady-state erosion rate:
                :math:`v_{s} = \\eta_{s}(0) v_{0}`
            v0_eqn_etas0_vs (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                erosion rate:
                :math:`v_{0} = \\dfrac{k v_{s}^{2}}{k v_{s} + w_{0}}`
            vs_eqn_w0_v0    (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                steady-state erosion rate:
                :math:`v_{s} = \\eta_{s}(0) v_{0}`
            v0_eqn_vs_w0    (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                erosion rate:
                :math:`v_{0} = \\dfrac{k v_{s}^{2}}{k v_{s} + w_{0}}`
            v0_eqn_vr_h_z   (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                erosion rate: 
                :math:`v_0 =  v_r \\left\{ (h-H_s[z,z_\\mathrm{vc},\\kappa_\\mathrm{v}])(1-v_b)+v_b \\right\}`
            w0_eqn_wr_z     (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                surface weakness:
                :math:`w_0 =  w_r H_s[z,z_\\mathrm{wc},\\kappa_w]`         
            
        """
        self.W_eqn           = Eq(W,w_0/(k*v_0))
        self.nus_eqn_W       = Eq(nu_s,(1+sy.sqrt(1+4*W))/2)
        self.nus_eqn_w0_v0   = self.nus_eqn_W.subs(W,self.W_eqn.rhs)
        self.etas0_eqn_W     = Eq(eta_s0, self.nus_eqn_W.rhs)
        self.etas0_eqn_w0_v0 = self.etas0_eqn_W.subs(W,self.W_eqn.rhs)
        self.nus_eqn_etas0   = Eq(nu_s, eta_s0)
        self.vs_eqn_etas0_v0 = Eq(v_s, self.nus_eqn_etas0.rhs*v_0)
        self.v0_eqn_etas0_vs = Eq(v_0, sy.solve(self.vs_eqn_etas0_v0,v_0)[0])
        self.vs_eqn_w0_v0    = Eq(v_s, self.nus_eqn_w0_v0.rhs*v_0)
        self.v0_eqn_vs_w0    = Eq(v_0, sy.solve(self.vs_eqn_w0_v0,v_0)[0])
        self.v0_eqn_vr_h_z   = Eq(v_0, 
                                  v_r*( (h-self.step(z,z_vc,kappa_v))*(1-v_b)+v_b ) 
                                  )
        self.w0_eqn_wr_z     = Eq(w_0, w_r*self.step(z,z_wc,kappa_w))


    @staticmethod
    def step(z,z0,kappa):
        """
        Step function at x0= :math:`z_0` and sharpness k= :math:`k`.
        
        :math:`H_s(z,z_0,\\kappa) = \\dfrac{1}{2}\\left(1 + \\tanh{[\\kappa(z-z_0)]}\\right)`

        Attributes:
            z (sympy.float)  : abscissa :math:`z`        
            z0 (float) : offset :math:`z_0`        
            kappa (float)  : step sharpness :math:`\\kappa`   
            
        Returns:
            float: step function :math:`H_s`               
        """
        return (1+sy.tanh((z-z0)*kappa))/2



        