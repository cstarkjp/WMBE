"""
Initialize WMBE package.
"""

from wmbe import initialize
from wmbe import channel1d
from wmbe import file
from wmbe import theory
from wmbe import viz_base

__version__ = "2026.6.14"

__all__ = [
    "initialize",
    "channel1d",
    "file",
    "model",
    "serialize",
    "solve",
    "symbols",
    "theory",
    "viz_base",
]