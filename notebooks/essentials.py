import os
import numpy as np
import sympy as sy
sy.init_printing(pretty_print=True,wrap_line=True)
from sympy import Eq

import wmbe
from wmbe.misc.file import (
    read_excel, create_directories, export_plots, export_plot,
)
from wmbe.numerical.model import WeatheringMediatedWeakness
from wmbe.numerical.solve1d import NumericalModel
from wmbe.applications.channel1d import ChannelWall
from wmbe.theory.symbols import *
from wmbe.theory.equations import Equations
from wmbe.viz.data import VizData
from wmbe.viz.simulations import VizSimulations
from wmbe.viz.applications import VizApplications
