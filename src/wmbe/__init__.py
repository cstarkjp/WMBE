"""
Initialize WMBE package.
"""

from wmbe import data
from wmbe import file
from wmbe import initialize
from wmbe import channel1d
from wmbe import theory
from wmbe import viz

__version__ = "2026.6.13"

__all__ = [
    "channel1d",
    "data",
    "file",
    "initialize",
    "model",
    "serialize",
    "solve",
    "symbols",
    "theory",
    "viz",
]