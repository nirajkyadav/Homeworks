"""
Material tester for the Menegotto-Pinto uniaxial material model.

Python port of MaterialTester_MP.m
    Developed for UNR CEE 721 by H. Ebrahimian.
    UNITS: kips, inches.

Run as a script to drive the model through a (random by default) cyclic
strain history and plot the strain history and stress-strain response.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import matplotlib.pyplot as plt

from mate_mp import MPMatData, MPHistory, MPState, mate_mp


def mat_colon(a: float, step: float, b: float) -> np.ndarray:
    """
    Replicate MATLAB's `a:step:b`:
      * inclusive endpoint when it lands exactly on a multiple of `step`;
      * empty array when `step` points away from `b`.
    """
    if step == 0:
        return np.array([a]) if a == b else np.array([])
    n = int(np.floor((b - a) / step + 1e-10)) + 1   # small tol guards the endpoint
    if n <= 0:
        return np.array([])
    return a + step * np.arange(n)


def build_strain_history(peak_strain: np.ndarray, deps: float) -> np.ndarray:
    """Build a triangular strain path that visits each peak in sequence."""
    eps = mat_colon(0.0, deps * np.sign(peak_strain[0]), peak_strain[0])
    for i in range(1, len(peak_strain)):
        sgn = 1.0 if peak_strain[i - 1] < peak_strain[i] else -1.0
        seg = mat_colon(peak_strain[i - 1] + deps * sgn, deps * sgn, peak_strain[i])
        eps = np.concatenate([eps, seg])
    return eps


def run_model(mat: MPMatData, eps: np.ndarray):
    """Push a strain history through the model, returning (strain, stress) arrays."""
    ey = mat.fy / mat.E
    init = MPHistory(sig=0.0, Et=mat.E, epmin=-ey, epmax=ey, epex=0.0,
                     ep0=0.0, s0=0.0, epr=0.0, sr=0.0, kon=0)
    state = MPState(past=replace(init), pres=replace(init))

    n = len(eps)
    strain = np.empty(n)
    stress = np.empty(n)
    for i in range(n):
        state.eps[0] = float(eps[i])
        state.eps[1] = 0.0 if i == 0 else float(eps[i] - eps[i - 1])
        state = mate_mp(mat, state)
        strain[i] = state.eps[0]
        stress[i] = state.pres.sig
        state.commit()                      # State.Past = State.Pres
    return strain, stress


def main(seed: int | None = None, show: bool = True, save: bool = False):
    # --- material parameters ---
    mat = MPMatData(E=29000.0, fy=36.0, b=0.01, R0=20.0,
                    cR1=0.925, cR2=0.15, a1=0.0, a2=0.0)

    epsy = mat.fy / mat.E
    deps = 0.02 * epsy

    # --- define peak strains ---
    # Option A (manual, mirrors the .m): uncomment to use fixed targeted peaks.
    peak_strain = np.array([2, -2, 4, -4, 6, -6, 8, -8, 0]) * epsy

    # # Option B (default): random peak strains within +/- 10 * yield strain.
    # # The original used MATLAB's unseeded rand(); pass an int seed for reproducibility.
    # rng = np.random.default_rng(seed)
    # max_strain_limit = 10 * epsy
    # peak_strain = rng.random(30) * 2 * max_strain_limit - max_strain_limit

    eps = build_strain_history(peak_strain, deps)
    strain, stress = run_model(mat, eps)

    # --- plot input strain history ---
    fig1, ax1 = plt.subplots()
    ax1.plot(eps * 100, linewidth=1.0)
    ax1.grid(True)
    ax1.set_xlabel("Step")
    ax1.set_ylabel("Strain (%)")
    fig1.tight_layout()

    # --- plot stress-strain response ---
    fig2, ax2 = plt.subplots()
    ax2.plot(strain * 100, stress, linewidth=1.0)
    ax2.grid(True)
    ax2.set_xlabel("Strain (%)")
    ax2.set_ylabel("Stress (ksi)")
    fig2.tight_layout()

    if save:
        fig1.savefig("Strain.png", dpi=150, bbox_inches="tight")
        fig2.savefig("StressStrain.png", dpi=150, bbox_inches="tight")
    if show:
        plt.show()

    return strain, stress, eps


if __name__ == "__main__":
    main()
