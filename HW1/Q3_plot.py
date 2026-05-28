import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ===========================================================================
# Material & Section Properties
# ===========================================================================
E = 29000.0             # ksi
nu = 0.3                # Poisson's ratio
G = E / (2 * (1 + nu))  # Shear modulus (ksi)

h = 20.04               # depth (in)
I = 3450.0              # Moment of Inertia (in^4)
A_s = 17.8              # Shear Area (in^2)

w_0 = 2.0               # kips/in

# Three beam lengths (Converted from ft to inches)
L_ft_list = [2, 5, 20]
L_in_list = [L * 12 for L in L_ft_list]     # [24, 60, 240] inches

# ===========================================================================
# Case a – Euler-Bernoulli
# ===========================================================================

def y_a(x, w_0, L, E, I):
    # Equation (1a)
    return (w_0 * L**2 / (16 * E * I)) * x**2 - (5 * w_0 * L / (48 * E * I)) * x**3 + (w_0 / (24 * E * I)) * x**4

def alpha_a(x, w_0, L, E, I):
    # Equation (2a)
    return (w_0 * L**2 / (8 * E * I)) * x - (5 * w_0 * L / (16 * E * I)) * x**2 + (w_0 / (6 * E * I)) * x**3

def gamma_a(x):
    # Equation (3a)
    return 0.0

def M_a(x, w_0, L):
    # Equation (4a)
    return (w_0 * L**2 / 8) - (5 / 8) * w_0 * L * x + (w_0 * x**2 / 2)

def V_a(x, w_0, L):
    # Equation (5a)
    return - (5 * w_0 * L / 8) + w_0 * x

# ===========================================================================
# Case b – Pure Shear
# ===========================================================================

def y_b(x, w_0, L, G, A_s):
    # Equation (1b)
    return (w_0 * L / (2 * G * A_s)) * x - (w_0 * x**2 / (2 * G * A_s))

def alpha_b(x):
    # Equation (2b)
    return 0.0

def gamma_b(x, w_0, L, G, A_s):
    # Equation (3b)
    return - (w_0 * L / (2 * G * A_s)) + (w_0 / (G * A_s)) * x

def M_b(x, w_0, L):
    # Equation (4b)
    return - (w_0 * L / 2) * x + (w_0 * x**2 / 2)

def V_b(x, w_0, L):
    # Equation (5b)
    return - (w_0 * L / 2) + w_0 * x

# ===========================================================================
# Case c – Timoshenko (general solution)
# ===========================================================================

def calculate_C2(w_0, L, E, I, G, A_s):
    # Equation for C2 
    numerator = (w_0 * L**2) / (2 * G * A_s) + (5 / 24) * (w_0 * L**4) / (E * I)
    denominator = - (2 * E * I * L) / (G * A_s) - (2 * L**3) / 3
    return numerator / denominator

def calculate_C1(C_2, w_0, L, E, I):
    # Equation for C1
    return -2 * C_2 * L - (w_0 * L**2) / (2 * E * I)

def y_c(x, w_0, E, I, G, A_s, C_1, C_2):
    # Equation (1)
    term1 = (-2 * E * I / (G * A_s)) * C_2 * x
    term2 = C_1 * (x**2 / 2)
    term3 = - (w_0 / (G * A_s)) * (x**2 / 2)
    term4 = C_2 * (x**3 / 3)
    term5 = (w_0 / (6 * E * I)) * (x**4 / 4)
    return term1 + term2 + term3 + term4 + term5

def alpha_c(x, w_0, E, I, C_1, C_2):
    # Equation (2)
    return C_1 * x + C_2 * (x**2) + (w_0 / (6 * E * I)) * (x**3)

def gamma_c(x, w_0, E, I, G, A_s, C_2):
    # Equation (3)
    return (2 * E * I / (G * A_s)) * C_2 + (w_0 / (G * A_s)) * x

def M_c(x, w_0, E, I, C_1, C_2):
    # Equation (4)
    return E * I * C_1 + (2 * C_2 * E * I) * x  + (w_0 * x**2) / 2

def V_c(x, w_0, E, I, C_2):
    # Equation (5)
    return 2 * E * I * C_2 + w_0 * x

# ===========================================================================
# Sweep parameters and evaluate at x = L/2
# ===========================================================================
# h/L sweep (vary L continuously to get a range of h/L values)
# Use 200 points from the smallest to largest L
L_sweep = np.array(L_in_list)
hL_sweep = h / L_sweep          # h/L

# phi sweep  phi = 12*E*I / (G*A_s*L^2)
phi_sweep = 12*E*I / (G*A_s*L_sweep**2)

# Storage: rows = [y, M, V], columns = [case_a, case_b, case_c]
y_a_vals, M_a_vals, V_a_vals = [], [], []
y_b_vals, M_b_vals, V_b_vals = [], [], []
y_c_vals, M_c_vals, V_c_vals = [], [], []

for L in L_sweep:
    x = L / 2          # evaluate at midspan

    # Case a
    y_a_vals.append(y_a(x, w_0, L, E, I))
    M_a_vals.append(M_a(x, w_0, L))
    V_a_vals.append(V_a(x, w_0, L))

    # Case b
    y_b_vals.append(y_b(x, w_0, L, G, A_s))
    M_b_vals.append(M_b(x, w_0, L))
    V_b_vals.append(V_b(x, w_0, L))

    # Case c
    C2 = calculate_C2(w_0, L, E, I, G, A_s)
    C1 = calculate_C1(C2, w_0, L, E, I)
    y_c_vals.append(y_c(x, w_0, E, I, G, A_s, C1, C2))
    M_c_vals.append(M_c(x, w_0, E, I, C1, C2))
    V_c_vals.append(V_c(x, w_0, E, I, C2))

y_a_vals = np.array(y_a_vals)
M_a_vals = np.array(M_a_vals)
V_a_vals = np.array(V_a_vals)

y_b_vals = np.array(y_b_vals)
M_b_vals = np.array(M_b_vals)
V_b_vals = np.array(V_b_vals)

y_c_vals = np.array(y_c_vals)
M_c_vals = np.array(M_c_vals)
V_c_vals = np.array(V_c_vals)


# ===========================================================================
# Plotting
# ===========================================================================

plt.figure(figsize=(14, 11))

# Row 1: Deflection (y)
plt.subplot(3, 2, 1)
plt.plot(hL_sweep, y_a_vals, label='(a) Euler-Bernoulli')
plt.plot(hL_sweep, y_b_vals, label='(b) Pure Shear')
plt.plot(hL_sweep, y_c_vals, label='(c) Timoshenko')
plt.xlabel(r'$h/L$')
plt.ylabel(r'$y(L/2) [in]$')
plt.title(r'Deflection Comparison')
plt.legend()

# Row 1: Deflection (y)
plt.subplot(3, 2, 2)
plt.plot(phi_sweep, y_a_vals, label='(a) Euler-Bernoulli')
plt.plot(phi_sweep, y_b_vals, label='(b) Pure Shear')
plt.plot(phi_sweep, y_c_vals, label='(c) Timoshenko')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$y(L/2) [in]$')
plt.title(r'Deflection Comparison')
plt.legend()

# Row 2: Shear Force (V)
plt.subplot(3, 2, 3)
plt.plot(hL_sweep, V_a_vals, label='(a) Euler-Bernoulli')
plt.plot(hL_sweep, V_b_vals, label='(b) Pure Shear')
plt.plot(hL_sweep, V_c_vals, label='(c) Timoshenko')
plt.xlabel(r'$h/L$')
plt.ylabel(r'$V(L/2) [kip]$')
plt.title(r'Shear Force Comparison')
plt.legend()

# Row 2: Shear Force (V)
plt.subplot(3, 2, 4)
plt.plot(phi_sweep, V_a_vals, label='(a) Euler-Bernoulli')
plt.plot(phi_sweep, V_b_vals, label='(b) Pure Shear')
plt.plot(phi_sweep, V_c_vals, label='(c) Timoshenko')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$V(L/2) [kip]$')
plt.title(r'Shear Force Comparison')
plt.legend()

# Row 3: Bending Moment (M)
plt.subplot(3, 2, 5)
plt.plot(hL_sweep, M_a_vals, label='(a) Euler-Bernoulli')
plt.plot(hL_sweep, M_b_vals, label='(b) Pure Shear')
plt.plot(hL_sweep, M_c_vals, label='(c) Timoshenko')
plt.xlabel(r'$h/L$')
plt.ylabel(r'$M(L/2) [kip \cdot in]$')
plt.title(r'Bending Moment Comparison')
plt.legend()

# Row 3: Bending Moment (M)
plt.subplot(3, 2, 6)
plt.plot(phi_sweep, M_a_vals, label='(a) Euler-Bernoulli')
plt.plot(phi_sweep, M_b_vals, label='(b) Pure Shear')
plt.plot(phi_sweep, M_c_vals, label='(c) Timoshenko')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$M(L/2) [kip \cdot in]$')
plt.title(r'Bending Moment Comparison')
plt.legend()

plt.tight_layout()
plt.show()