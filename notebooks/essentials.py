import os
import numpy as np
import sympy as sy
sy.init_printing(pretty_print=True,wrap_line=True)
from sympy import Eq

import wmbe
from wmbe.data import ExperimentalData, read_excel
from wmbe.file import create_directories, export_plots, export_plot
from wmbe.model import WeatheringMediatedWeakness
from wmbe.solve1d import NumericalModel
from wmbe.channel1d import ChannelWallApplication
from wmbe.symbols import *
from wmbe.theory import Equations
from wmbe.viz import DataViz, SimViz, AppViz