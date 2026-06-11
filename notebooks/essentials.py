import os
import numpy as np
import sympy as sy
sy.init_printing(pretty_print=True,wrap_line=True)
from sympy import Eq

import wmbe
from wmbe.data import ExperimentalData
from wmbe.file import create_directories, export_plots, export_plot
from wmbe.solve1d import ErosionWeathering
from wmbe.symbols import *
from wmbe.theory import WeatheringMediatedErosion
from wmbe.viz import Viz