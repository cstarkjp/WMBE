"""
Initialize WMBE package.
"""

from wmbe import file
from wmbe import initialize
from wmbe import channel1d
from wmbe import theory
from wmbe import viz_base

__version__ = "2026.6.13"

__all__ = [
    "channel1d",
    "file",
    "initialize",
    "model",
    "serialize",
    "solve",
    "symbols",
    "theory",
    "viz_base",
]