"""
==============================================================================
CEE 721: Nonlinear Structural Analysis

Date: July 13, 2026
Author: Niraj Yadav

HW 4, Problem 1: Curve Tracing

Units: kips, inches, ksi

==============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

My = 320.0   # kip-in
Ky = 0.001   # 1/in

# def moment_curvature(kappa):
#     return (320/np.sqrt(0.001))*np.sqrt(kappa) # Parabolic

def moment_curvature(kappa):
    """
    Calculates Moment (M) for a given Curvature (kappa) using a 
    circular curve fitted in normalized coordinate space.
    
    DERIVATION:
       Let x = kappa / Ky  
       Let y = M / My     

       Equation: (x - 1)^2 + y^2 = 1^2

       y   = sqrt(1 - (x^2 - 2x + 1))
       y   = sqrt(2*x - x^2)
       
       M = My * y
    """

    x = kappa / Ky
    y_circular = np.sqrt(np.clip(2 * x - x**2, 0.0, None))       # np.clip prevents negative numbers inside sqrt
    
    M = My * y_circular   # = 320*((2000*kappa-1000_000*(kappa**2))**0.5)
    
    return np.where(kappa <= Ky, M, My)    # Circular

kappa = np.linspace(0.0, 0.0022, 800)
M_vals = moment_curvature(kappa)

fig, ax = plt.subplots(figsize=(6, 6))

ax.plot(kappa, M_vals, color = 'darkcyan', linewidth=4)

# Limits
ax.set_xlim(0, 0.0025)
ax.set_ylim(0, 400)

# Major grid spacing:
# 1 square horizontally = 0.001 (1/in)
# 1 square vertically   = 320 kip-in
ax.xaxis.set_major_locator(MultipleLocator(0.001))
ax.yaxis.set_major_locator(MultipleLocator(320))

ax.grid(True, linestyle='--', alpha=0.5)

# Make one x-grid interval equal to one y-grid interval
# (0.001 curvature ↔ 320 kip-in)
ax.set_aspect(Ky / My)

ax.set_xlabel(r'Curvature, $\kappa$ (1/in)')
ax.set_ylabel(r'Moment, $M$ (kip-in)')
ax.set_title('Moment-Curvature Relation')

plt.tight_layout()
plt.show()