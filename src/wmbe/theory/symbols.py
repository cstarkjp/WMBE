"""
Mathematical symbols for parameters. 

Make for pretty printing of parameter dict contents in Jupyter.

Attributes:
    W (Symbol): x
    w_0 (Symbol): x
    v_0 (Symbol): x
    k (Symbol): x
    nu_s (Symbol): x
    v_s (Symbol): x

    eta (Symbol): x
    eta_s (Symbol): x
    eta_s0 (Symbol): x

    v_r (Symbol): x
    h (Symbol): x
    z (Symbol): x
    w_r (Symbol): x
    v_b (Symbol): x

    kappa_v (Symbol): x
    kappa_w (Symbol): x
    z_vc (Symbol): x
    z_wc (Symbol): x
    nu_s_bar (Symbol): x
    Δχ (Symbol): x
    Δτ (Symbol): x
    x_domain_size (Symbol): x
    τ_domain_size (Symbol): x
"""

from sympy import Symbol, symbols

W: Symbol = Symbol("W", positive=True)
w_0: Symbol = Symbol("w_0", positive=True)
u_0: Symbol = Symbol("v_0", positive=True)
k: Symbol = Symbol("k", positive=True)
ν_s: Symbol = Symbol("\\nu_s", positive=True)
u_s: Symbol = Symbol("u_s", positive=True)

eta: Symbol = Symbol("\\eta", positive=True)
eta_s: Symbol = Symbol("\\eta_{s0}", positive=True)
eta_s0: Symbol = Symbol("\\eta_{s0}", positive=True)

v_r: Symbol = Symbol("v_r", positive=True)
h: Symbol = Symbol("h", positive=True)
z: Symbol = Symbol("z", positive=True)
w_r: Symbol = Symbol("w_r", positive=True)
v_b: Symbol = Symbol("v_b", positive=True)

kappa_v: Symbol = Symbol("\\kappa_v", positive=True)
kappa_w: Symbol = Symbol("\\kappa_w", positive=True)
z_vc: Symbol = Symbol("z_\\mathrm{vc}", positive=True)
z_wc: Symbol = Symbol("z_\\mathrm{wc}", positive=True)
nu_s_bar: Symbol = Symbol("\\overline{\\nu}_s", positive=True)
Δχ: Symbol = Symbol("\\Delta\\chi", positive=True)
Δτ: Symbol = Symbol("\\Delta\\tau", positive=True)
x_domain_size: Symbol = Symbol("\\chi_\\mathrm{domain}", positive=True)
τ_domain_size: Symbol = Symbol("\\tau_\\mathrm{domain}", positive=True)
