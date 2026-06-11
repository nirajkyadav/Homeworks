"""
Menegotto-Pinto uniaxial material model with isotropic and kinematic
hardening, following Filippou et al. (EERC83-19).

Python port of MateMP.m
    FEDEAS Lab, Release 2.3 (March 2001)
    Copyright (c) 1998, F. C. Filippou, UC Berkeley
    MP routine contributed by P. Franchin & A. Lupoi (January 2001)

The model traces a smooth (Bauschinger) transition between two asymptotes,
one of slope E and one of slope b*E, whose curvature is governed by R.
Cross-references to the EERC83-19 equations (as numbered in the accompanying
course handout) are given inline below.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import List, Tuple


@dataclass
class MPMatData:
    """Material parameters for the Menegotto-Pinto model."""
    E: float                 # initial elastic modulus
    fy: float                # yield strength (initial)
    b: float = 0.01          # strain-hardening ratio (Es2 = b*E)
    R0: float = 20.0         # initial curvature of the elastic->plastic transition
    cR1: float = 0.925       # coefficient for the variation of R   (R1 in EERC83-19)
    cR2: float = 0.15        # coefficient for the variation of R   (R2 in EERC83-19)
    a1: float = 0.0          # isotropic-hardening coefficient
    a2: float = 0.0          # isotropic-hardening coefficient
    # Note: the original MateMP.m header lists cR1's default as 18.5. That is
    # almost certainly a typo (it would drive R negative for moderate xi).
    # 0.925 is the standard Filippou / OpenSees Steel02 value and is what
    # MaterialTester_MP.m uses.


@dataclass
class MPHistory:
    """History variables at a converged (`past`) or trial (`pres`) state."""
    sig: float = 0.0     # stress
    Et: float = 0.0      # tangent modulus
    epmin: float = 0.0   # max strain in compression
    epmax: float = 0.0   # max strain in tension
    epex: float = 0.0    # plastic excursion (extreme strain used for xi)
    ep0: float = 0.0     # strain at asymptote intersection (updated yield strain)
    s0: float = 0.0      # stress at asymptote intersection (updated yield stress)
    epr: float = 0.0     # strain at last reversal
    sr: float = 0.0      # stress at last reversal
    kon: int = 0         # 0 = virgin, 1 = loading (tension), 2 = unloading (compression)


@dataclass
class MPState:
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
    past: MPHistory = field(default_factory=MPHistory)
    pres: MPHistory = field(default_factory=MPHistory)
    sig: float = 0.0
    Et: float = 0.0

    def commit(self) -> None:
        """Accept the trial state as converged (Past <- Pres)."""
        self.past = replace(self.pres)


def mate_mp(mat: MPMatData, state: MPState) -> MPState:
    """Advance one strain step. Reads `state.past` + `state.eps`, writes `state.pres`."""
    b, R0, cR1, cR2 = mat.b, mat.R0, mat.cR1, mat.cR2
    a1, a2, E, fy = mat.a1, mat.a2, mat.E, mat.fy

    Es2 = b * E      # hardening modulus
    ey = fy / E      # initial yield strain (eps_y0)

    # retrieve history at the last converged state
    p = state.past
    sig, Et = p.sig, p.Et
    epmin, epmax = p.epmin, p.epmax
    epex, ep0, s0 = p.epex, p.ep0, p.s0
    epr, sr, kon = p.epr, p.sr, p.kon
    sigp = sig                       # stress at last converged state

    epm = max(abs(epmin), abs(epmax))
    epss = state.eps[0]              # total strain at current iteration
    depss = state.eps[1]             # total strain increment from last convergence

    if kon == 0:                                     # ---- virgin material ----
        if depss == 0:
            sig, Et = 0.0, E
        if depss > 0:                                # first loading in tension
            kon = 1
            ep0, s0, epex = epm, fy, epm
            sig, Et = _bauschinger(epex, ep0, ey, R0, cR1, cR2, epss, epr, b, s0, sr)
        if depss < 0:                                # first loading in compression
            kon = 2
            ep0, s0, epex = -epm, -fy, -epm
            sig, Et = _bauschinger(epex, ep0, ey, R0, cR1, cR2, epss, epr, b, s0, sr)
    else:                                            # ---- damaged material ----
        if (kon == 1 and depss > 0) or (kon == 2 and depss < 0):
            # continue loading in the current direction along the same asymptote
            sig, Et = _bauschinger(epex, ep0, ey, R0, cR1, cR2, epss, epr, b, s0, sr)

        elif kon == 1 and depss < 0:                 # reversal: tension -> compression
            kon = 2
            epr = epss - depss                       # last converged strain == reversal point
            sr = sigp
            if epr > epmax:
                epmax = epr                          # update tensile extreme (Eq. 7)
            epm = max(abs(epmin), abs(epmax))
            sst = max(fy * a1 * (epm / ey - a2), 0.0)            # isotropic shift (Eq. 15)
            ep0 = (sr + fy + sst - (E * epr + Es2 * ey)) / (Es2 - E)   # updated eps_y (Eq. 5)
            s0 = Es2 * (ep0 + ey) - fy - sst                          # updated sig_y (Eq. 6)
            epex = epmin
            sig, Et = _bauschinger(epex, ep0, ey, R0, cR1, cR2, epss, epr, b, s0, sr)

        elif kon == 2 and depss > 0:                 # reversal: compression -> tension
            kon = 1
            epr = epss - depss
            sr = sigp
            if epr < epmin:
                epmin = epr                          # update compressive extreme (Eq. 11)
            epm = max(abs(epmin), abs(epmax))
            sst = max(fy * a1 * (epm / ey - a2), 0.0)            # isotropic shift (Eq. 15)
            ep0 = (sr + Es2 * ey - (E * epr + fy + sst)) / (Es2 - E)   # updated eps_y (Eq. 9)
            s0 = fy + sst + Es2 * (ep0 - ey)                          # updated sig_y (Eq. 10)
            epex = epmax
            sig, Et = _bauschinger(epex, ep0, ey, R0, cR1, cR2, epss, epr, b, s0, sr)
        # else (depss == 0): no change; sig/Et retain their converged values

    # store the trial (Pres) state
    state.pres = MPHistory(sig=sig, Et=Et, epmin=epmin, epmax=epmax, epex=epex,
                           ep0=ep0, s0=s0, epr=epr, sr=sr, kon=kon)
    state.sig = sig
    state.Et = Et
    return state


def _bauschinger(epex: float, ep0: float, epy: float, R0: float, cR1: float,
                 cR2: float, epss: float, epr: float, b: float, s0: float,
                 sr: float) -> Tuple[float, float]:
    """
    Evaluate the normalized Menegotto-Pinto curve and its tangent.

        sigma* = b*eps* + (1-b)*eps* / (1 + |eps*|^R)^(1/R)        (Eq. 1)

    where eps* = (eps - epr) / (ep0 - epr) and the result is de-normalized
    onto the asymptote pair (epr, sr) -> (ep0, s0). `Et` is the closed-form
    derivative d(sigma)/d(eps).
    """
    xi = abs((epex - ep0) / epy)
    R = R0 * (1.0 - (cR1 * xi) / (cR2 + xi))        

    e_str = (epss - epr) / (ep0 - epr)
    abs_e_R = abs(e_str) ** R
    denom = (1.0 + abs_e_R) ** (1.0 / R)

    s_str = b * e_str + (1.0 - b) * e_str / denom
    sig = s_str * (s0 - sr) + sr

    dSdE = b + (1.0 - b) * (1.0 - abs_e_R / (1.0 + abs_e_R)) / denom
    Et = dSdE * (s0 - sr) / (ep0 - epr)

    return sig, Et
