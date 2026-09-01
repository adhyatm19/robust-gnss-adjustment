"""Robust GNSS baseline network adjustment.

Synthetic baseline-level research prototype comparing classical weighted
least squares, DIA (Detection-Identification-Adaptation), specific-direction
outlier detection, and robust M-estimation (Huber, Hampel). Not a raw
RINEX/carrier-phase processor.
"""

from .dia import DIAConfig, run_dia
from .robust import RobustConfig, hampel_weight, huber_weight, robust_adjust
from .sd_outlier import (
    SDConfig,
    SDSnoopingResult,
    SpecificDirectionDetection,
    SpecificDirectionResult,
    angular_error_degrees,
    detect_specific_direction,
    run_sd_snooping,
    test_specific_direction,
)
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
    "SDConfig",
    "SDSnoopingResult",
    "Station",
    "SpecificDirectionDetection",
    "SpecificDirectionResult",
    "angular_error_degrees",
    "detect_specific_direction",
    "hampel_weight",
    "huber_weight",
    "robust_adjust",
    "run_dia",
    "run_sd_snooping",
    "simulate_network",
    "solve_wls",
    "test_specific_direction",
]
