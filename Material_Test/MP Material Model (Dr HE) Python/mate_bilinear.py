"""
Bilinear hysteretic uniaxial material model with kinematic hardening.

Python port of MateBilinear.m
    FEDEAS Lab, Release 2.2 (January 2001)
    Copyright (c) 1998, F. C. Filippou, UC Berkeley
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field, replace
from typing import List, Tuple

# Matches MATLAB's built-in scalar `eps` (2^-52).
_EPS = sys.float_info.epsilon


@dataclass
class BilinearMatData:
    """Material parameters for the bilinear model."""
    E: float    # initial elastic modulus
    fy: float   # yield stress
    Eh: float   # strain-hardening modulus


@dataclass
class BilinearHistory:
    """History variables at a converged (`past`) or trial (`pres`) state."""
    sig: float = 0.0    # stress
    Et: float = 0.0     # tangent modulus


@dataclass
class BilinearState:
    """
    Full material state handed to / returned from the model each step.

    `eps` mirrors the FEDEAS convention:
        eps[0] = total strain
        eps[1] = strain increment from last convergence
        eps[2] = strain increment from last iteration
        eps[3] = strain rate
    The model only consumes eps[0] and eps[1].
    """
    eps: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0, 0.0])
    past: BilinearHistory = field(default_factory=BilinearHistory)
    pres: BilinearHistory = field(default_factory=BilinearHistory)
    sig: float = 0.0
    Et: float = 0.0

    def commit(self) -> None:
        """Accept the trial state as converged (Past <- Pres)."""
        self.past = replace(self.pres)


def mate_bilinear(mat: BilinearMatData, state: BilinearState) -> BilinearState:
    """Advance one strain step. Reads `state.past` + `state.eps`, writes `state.pres`."""
    sig, Et = _state_determination(mat, state)
    state.sig = sig
    state.Et = Et
    state.pres = BilinearHistory(sig=sig, Et=Et)
    return state


def _state_determination(mat: BilinearMatData,
                         state: BilinearState) -> Tuple[float, float]:
    """State determination for the bilinear hysteretic stress-strain relation."""
    E, fy, Eh = mat.E, mat.fy, mat.Eh
    toler = _EPS

    epsy = fy / E
    sigP = state.past.sig            # stress  at last convergence
    Ep = state.past.Et               # modulus at last convergence
    epsC = state.eps[0]              # current strain
    deps = state.eps[1]              # strain increment from last convergence

    if deps == 0:
        return sigP, Ep

    if deps > 0:
        # min(hardening branch, elastic predictor)
        sig = min(fy + (epsC - epsy) * Eh, sigP + deps * E)
    else:
        # max(hardening branch, elastic predictor)
        sig = max(-fy + (epsC + epsy) * Eh, sigP + deps * E)

    # elastic if the elastic predictor governs, otherwise on the hardening branch
    Et = E if abs(sig - (sigP + deps * E)) <= toler else Eh
    return sig, Et
