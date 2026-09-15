"""Thermodynamic backends.

pyGCC is the only backend implemented. The abstract base exists so that the
numerical core never imports pyGCC directly, which keeps the core testable
against a stub and documents exactly what a backend must provide.
"""

from .base import ThermoBackend
from .pygcc_backend import PygccBackend

__all__ = ["ThermoBackend", "PygccBackend"]
