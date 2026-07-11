"""
==============================================================================
CEE 721: Nonlinear Structural Analysis

Date: July 11, 2026
Author: Niraj Yadav

Problem 3: Comparison of the vertical deformation y(x), internal shear V(x), and internal moment M(x) for the three cases,

(a) Euler-Bernoulli beam theory
(b) pure shear beam theory
(c) Timoshenko beam theory

Units: kips, inches, ksi

==============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt

# Material & Section Properties
E = 29000.0             # ksi
nu = 0.3                # Poisson's ratio
G = E / (2 * (1 + nu))  # Shear modulus (ksi)

h = 20.04               # depth (in)
I = 3450.0              # Moment of Inertia (in^4)
As = 17.8               # Shear Area (in^2)
w0 = 2.0                # kips/in (Upward load)

L_ft_list = [2, 5, 20]  # Three values for span length in ft

fig, axs = plt.subplots(3, 3, figsize=(15, 12))

for idx, L_ft in enumerate(L_ft_list):

    L = L_ft * 12.0                         # Converted length to inches
    x = np.linspace(0, L, 200)              # Position Vectors: Division of L in 200 points (in inches)
    
    # Parameters calculation 
    hoL = h / L                             # Ratio of depth to span length
    phi = (12 * E * I) / (G * As * L**2)    # Shear Correction Term
    
    # (a) Euler-Bernoulli beam theory
    y_a = (w0 * L**2 / (16 * E * I)) * x**2 - (5 * w0 * L / (48 * E * I)) * x**3 + (w0 / (24 * E * I)) * x**4
    M_a = (w0 * L**2 / 8) - (5 / 8) * w0 * L * x + (w0 * x**2 / 2)
    V_a = - (5 / 8) * w0 * L + w0 * x
    
    # (b) Pure Shear beam theory
    y_b = (w0 * L / (2 * G * As)) * x - (w0 * x**2 / (2 * G * As))
    M_b = -(w0 * L / 2) * x + (w0 * x**2 / 2)
    V_b = -(w0 * L / 2) + w0 * x
    
    # (c) Timoshenko beam theory
    num = (w0 * L**2) / (2 * G * As) + (5 / 24) * (w0 * L**4) / (E * I)
    den = - (2 * E * I * L) / (G * As) - (2 * L**3) / 3
    a2 = num / den
    a1 = -2 * a2 * L - (w0 * L**2) / (2 * E * I)
    
    y_c = (-2 * E * I / (G * As)) * a2 * x + a1 * (x**2 / 2) - (w0 / (G * As)) * (x**2 / 2) + a2 * (x**3 / 3) + (w0 / (24 * E * I)) * x**4
    M_c = E * I * a1 + 2 * a2 * E * I * x + 0.5 * w0 * x**2
    V_c = 2 * E * I * a2 + w0 * x

    x_ft = x / 12.0                        # Convert x to feet for horizontal axis

    # Row 1: Deflection
    axs[0, idx].plot(x_ft, y_a, label='Euler-Bernoulli')
    axs[0, idx].plot(x_ft, y_b, label='Pure Shear')
    axs[0, idx].plot(x_ft, y_c, label='Timoshenko', linestyle='--')
    axs[0, idx].set_title(f'L = {L_ft} ft (h/L={hoL:.3f}, $\\phi$={phi:.2f})')
    axs[0, idx].set_ylabel('Deflection $y(x)$ [in]')
    axs[0, idx].grid(alpha=0.3)
    
    # Row 2: Shear Force
    axs[1, idx].plot(x_ft, V_a)
    axs[1, idx].plot(x_ft, V_b)
    axs[1, idx].plot(x_ft, V_c, linestyle='--')
    axs[1, idx].set_ylabel('Shear $V(x)$ [kips]')
    axs[1, idx].grid(alpha=0.3)
    
    # Row 3: Bending Moment
    axs[2, idx].plot(x_ft, M_a)
    axs[2, idx].plot(x_ft, M_b)
    axs[2, idx].plot(x_ft, M_c, linestyle='--')
    axs[2, idx].set_ylabel('Moment $M(x)$ [kip·in]')
    axs[2, idx].set_xlabel('Position $x$ [ft]')
    axs[2, idx].grid(alpha=0.3)

# Add legend to the top-left plot
axs[0, 0].legend()
plt.tight_layout()
plt.savefig("Q3_plots.pdf")
plt.show()