"""
Model Inoue & Li experimental data.

Inoue et al (2017`: https://doi.org/10.1016/j.geomorph.2017.02.018
Li et al (2016): https://doi.org/10.25103/jestr.093.10
"""
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from numpy.typing import NDArray
from pandas import DataFrame

warnings.filterwarnings("ignore")

__all__ = [
    "linear_model",
    "exponential_decay_model",
    "weakening_model",
    "WeatheringMediatedWeakness",
]

def linear_model(
        x: float|NDArray, 
        m: float, 
        c: float,
    ) -> float|NDArray:
    """
    Simple linear model of form: $y = m x + c$.
    """    
    # Args:
    #     x (float or NDArray) : coordinate
    #     m (float) : gradient
    #     c (float) : intercept

    # Returns:
    #     float or NDArray: y
    return m*x+c

def exponential_decay_model(
        x: float|NDArray, 
        m: float, 
        c: float,
    ) -> float|NDArray:
    r"""
    Shifted exponential decay model of form: $y = 1 + c \exp(-x/m)$.
    """    
    # Args:
    #     x (float or NDArray) : coordinate
    #     m (float) : e-folding scale
    #     c (float) : magnitude

    # Returns:
    #     float or NDArray: y
    return 1 + c*np.exp(-x/m)

def weakening_model(
        wetdryN_P: float, 
        k: float, 
        w0: float, 
        tau0: float,
    ) -> float:
    r"""
    Shifted exponential decay weathering model: 
        $w = 1 + w_0(\\tau+\\tau_0)\exp(-k\\chi)$
    """    
    # Args:
    #     wetdryN_P (numpy.ndarray) : pair $(\\tau,\chi)$
    #     k (float) : reciprocal e-folding scale $k$
    #     w0 (float) : magnitude $w_0$
    #     tau0 (float) : time offset $\\tau_0$

    # Returns:
    #     float: y
    tau = wetdryN_P[0]
    chi = wetdryN_P[1]
    return 1 + w0*(tau+tau0)*np.exp(-k*chi)

class WeatheringMediatedWeakness:
    """
    Rock weathering & erosion modeling.
    """
    def __init__(self) -> None:
        """
        """
        self.fits = dict()
        self.fits2d = dict()

    def fit_weakness_vs_time_linear_model(
            self, 
            data: DataFrame, 
            x_name: str, 
            y_name: str, 
            select: str|None=None,
        ) -> None:
        """
        Regress a 1d model against experimental data.
        
        Perform a linear regression fit of a linear model to the given
        experimental data, such as modeling the degree of rock weakness
        as a linear function of the number of wetting and drying cycles.
        """
        
        # Args:
        #     data_set (str) : which experimental dataset,
        #         as key to ddict element :class:`pandas.DataFrame`
        #     x_name (str) : abscissa :math:`x`
        #     y_name (str) : ordinate :math:`y`
        #     select (str) : which computation of weakness from rock strength
    
        # Attributes:
        #     fdict[data_set] or fdict[selection_name]  (:obj:`dict` element) : 
        #         model fit 
        #     w_s2_means (:class:`numpy.ndarray`) :
        #         mean values of weakness :math:`w`
        #     w_s2_stds (:class:`numpy.ndarray`) :
        #         standard deviations of weakness :math:`w`
        if select is not None:
            for selection in np.unique(data[select]):
                selection_name = f"{select}_{selection}"
                x = data.loc[data[select]==selection][x_name]
                y = data.loc[data[select]==selection][y_name]
                self.fits[selection_name] = curve_fit(linear_model, x, y,)
        else:
            x = data[x_name]
            y = data[y_name]
            self.fits["default"] = curve_fit(linear_model, x, y,)
            
        self.w_s2_means \
            = data.groupby("wetdryN").mean()["w_sigma2"]
        self.w_s2_stds \
            = data.groupby("wetdryN").std()[ "w_sigma2"]

    def fit_weakness_vs_time_and_depth_model(
            self, 
            data: DataFrame, 
            select: str
        ) -> None:
        """
        Regress a 2d model against experimental data.
        """
        # Args:
        #     data_set (:obj:`str`) : which experimental dataset,
        #         as key to ddict element :class:`pandas.DataFrame`
        #     select (:obj:`str`) : which computation of weakness from rock strength

        # Attributes:
        #     fdict[data_set] (:obj:`list`) : model surface fit as
        #         meshgrid X=wet/dry :math:`N`, Y=confining pressure :math:`P`
        #         and corresponding surface Z[X,Y] as list (X,Y,Z)
        #     sdict[data_set] (:obj:`list`) : model surface estimates at experimental values as
        #         meshgrid X=wet/dry :math:`N`, Y=confining pressure :math:`P`
        #         and corresponding surface Z[X,Y] as list (X,Y,Z)
        #     w_s2normed_means (:class:`numpy.ndarray`) : 
        #         mean values of model-normed weakness :math:`w`
        #     w_s2normed_stds (:class:`numpy.ndarray`) : 
        #         standard deviations of :math:`w`
    
        
        wdN_vec  = data.wetdryN
        P_vec    = data.P
        sig_vec  = data.sigmaC/180
        w_vec    = sig_vec**(-2)
        wdN_P_array = np.vstack((
            wdN_vec,
            P_vec,
        ))
        model_fit   = curve_fit(weakening_model, wdN_P_array, w_vec,)
        self.fits["default"] = model_fit
    
        n_pts = 30
        X = np.linspace(0, wdN_vec.max()*1.1, n_pts,)
        Y = np.linspace(0, P_vec.max()*1.1, n_pts,)
        X,Y = np.meshgrid(X, Y,)
        X_Y_array = np.vstack((
            X.reshape(n_pts**2), 
            Y.reshape(n_pts**2),
        ))
        Z = weakening_model(X_Y_array, *model_fit[0],).reshape(n_pts, n_pts,)
        self.fits[select+"_surface"] = (X, Y, Z,)
        
        P_vec = np.linspace(0,P_vec.max()*1.3)
        X,Y = np.meshgrid(np.unique(wdN_vec), P_vec)
        X_Y_array = np.vstack((
            X.reshape(X.shape[0]*Y.shape[1]), 
            Y.reshape(X.shape[0]*Y.shape[1]),
        ))
        self.fits2d["default"] = (
            X[0], 
            Y[:,0],
            weakening_model(
                X_Y_array, 
                *model_fit[0],
            ).reshape(X.shape[0], Y.shape[1],)
        )
        
        pd.options.mode.chained_assignment = None
        w_ref_vec = (self.fits2d["default"][2].T.copy())[:,0] - 1
        data["w_s2normed"] = 0.0
        for idx,wetdryN in enumerate(np.unique(data.wetdryN)):
            w = data.w_sigma2[data.wetdryN==wetdryN].copy()
            w_normed = (w-1)/w_ref_vec[idx]+1
            data.loc[(idx*3):(idx*3+3), "w_s2normed"] = w_normed

        self.w_s2normed_means \
            = data.groupby("P").mean()["w_s2normed"]
        self.w_s2normed_stds \
            = data.groupby("P").std()["w_s2normed"]