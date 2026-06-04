"""
Initialize WMBE package.
"""

from wmbe import data
from wmbe import initialize
from wmbe import plot
from wmbe import save
from wmbe import solve1p1d
from wmbe import theory

__version__ = "2026.6.07"

__all__ = [
    "data",
    "initialize",
    "plot",
    "save",
    "solve",
    "solve1p1d",
    "symbols",
    "theory",
]