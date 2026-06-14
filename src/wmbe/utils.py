"""
Math utility functions.

Inoue et al (2017`: https://doi.org/10.1016/j.geomorph.2017.02.018
Li et al (2016): https://doi.org/10.25103/jestr.093.10
"""
import warnings
import numpy as np

from numpy.typing import NDArray

warnings.filterwarnings("ignore")

__all__ = [
    "linear_model",
    "exponential_decay_model",
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
    """
    Shifted exponential decay model of form: $y = 1 + c \\exp(-x/m)$.
    """    
    # Args:
    #     x (float or NDArray) : coordinate
    #     m (float) : e-folding scale
    #     c (float) : magnitude

    # Returns:
    #     float or NDArray: y
    return 1 + c*np.exp(-x/m)
