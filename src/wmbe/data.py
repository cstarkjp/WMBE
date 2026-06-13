"""
Import and model Inoue & Li experimental data.

Inoue et al (2017`: https://doi.org/10.1016/j.geomorph.2017.02.018
Li et al (2016): https://doi.org/10.25103/jestr.093.10
"""
import warnings
import numpy as np
import os
import pandas as pd
from scipy.optimize import curve_fit

from typing import Any, Callable
from collections.abc import Sequence
from numpy.typing import NDArray

warnings.filterwarnings("ignore")

def read_excel(
        self, 
        dir_name=("..", "data",), 
        file_name="Inoue_wetdryN_sigmaT",
        header=0, 
        skiprows=[1],
    ) -> None:
    """
    Read data from an Excel file into a :mod:`pandas` dataframe
    
    Attributes:
        data_set  (:obj:`str`)  : name of dataset in data dictionary
        dir_name  (:obj:`list`) : data source directory as platform-indepedent list
        file_name (:obj:`str`)  : data source filename
        header    (:obj:`int`)  : number of header (column title etc) rows
        skiprows  (:obj:`int`)  : number of rows to skip when creating dataframe

    """
    dir_name = os.path.join(*dir_name)
    if not os.path.exists(dir_name):
        print("Cannot find data directory")
        raise
    try:
        df = pd.read_excel(
            os.path.join(dir_name, file_name+".xlsx",),
            header=header, 
            skiprows=skiprows,
        )
    except OSError:  
        print("Cannot find data directory")
        raise
    except:  
        raise
    # self.data.update({data_set:df})
    self.df = df
    
class ExperimentalData:
    """
    Rock weathering experimental data import and modeling.
    """
    def __init__(self) -> None:
        """
        Initialize class instance.

        Attributes:
            data (:obj:`dict`) : data dictionary, to contain experimental results of 
                `Inoue et al (2017)`_ and `Li et al (2016)`_ variously on 
                rock tensile strength, rock compressive strength,
                erodibility, number of wetting/drying cycles, and confining pressure
            fits (:obj:`dict`) : 1d model dictionary, to contain regression fits of
                weathering model(s) to experimental data
            fits2d (:obj:`dict`) : 2d model dictionary, to contain regression fits of
                2d weathering model to experimental data
        """
        self.df = None
        
    def read_excel(
            self, 
            dir_name=("..", "data",), 
            file_name="Inoue_wetdryN_sigmaT",
            header=0, 
            skiprows=[1],
        ) -> None:
        """
        Read data from an Excel file into a :mod:`pandas` dataframe
        
        Attributes:
            data_set  (:obj:`str`)  : name of dataset in data dictionary
            dir_name  (:obj:`list`) : data source directory as platform-indepedent list
            file_name (:obj:`str`)  : data source filename
            header    (:obj:`int`)  : number of header (column title etc) rows
            skiprows  (:obj:`int`)  : number of rows to skip when creating dataframe
    
        """
        dir_name = os.path.join(*dir_name)
        if not os.path.exists(dir_name):
            print("Cannot find data directory")
            raise
        try:
            df = pd.read_excel(
                os.path.join(dir_name, file_name+".xlsx",),
                header=header, 
                skiprows=skiprows,
            )
        except OSError:  
            print("Cannot find data directory")
            raise
        except:  
            raise
        # self.data.update({data_set:df})
        self.df = df