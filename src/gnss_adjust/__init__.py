"""Robust GNSS baseline network adjustment.

Synthetic baseline-level research prototype comparing classical weighted
least squares, DIA (Detection-Identification-Adaptation) outlier handling,
and robust M-estimation (Huber, Hampel). Not a raw RINEX/carrier-phase
processor.
"""

from .dia import DIAConfig, run_dia
from .robust import RobustConfig, hampel_weight, huber_weight, robust_adjust
from .simulation import NetworkConfig, NoiseModel, OutlierSpec, simulate_network
from .types import AdjustmentResult, BaselineObservation, DIAResult, RobustResult, Station
from .wls import solve_wls

__version__ = "0.1.0"

__all__ = [
    "AdjustmentResult",
    "BaselineObservation",
    "DIAConfig",
    "DIAResult",
    "NetworkConfig",
    "NoiseModel",
    "OutlierSpec",
    "RobustConfig",
    "RobustResult",
    "Station",
    "hampel_weight",
    "huber_weight",
    "robust_adjust",
    "run_dia",
    "simulate_network",
    "solve_wls",
]
