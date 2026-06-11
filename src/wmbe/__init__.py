"""
Initialize WMBE package.
"""

from wmbe import data
from wmbe import file
from wmbe import initialize
from wmbe import plot
from wmbe import solve1p1d
from wmbe import theory
from wmbe import viz

__version__ = "2026.6.11"

__all__ = [
    "data",
    "file",
    "initialize",
    "plot",
    "serialize",
    "solve",
    "solve1p1d",
    "symbols",
    "theory",
    "viz",
]