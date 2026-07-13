"""
==============================================================================
CEE 721: Nonlinear Structural Analysis

Updated Date: July 13, 2026
Author: Niraj Yadav

HW 4, Problem 3: Incremental Iterative Solution

Units: kips, inches, ksi

==============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from mate_mp import MPMatData, MPState, mate_mp

# the material model mate_mp turns ε into stress σ and tangent modulus Et = dσ/dε

# ==========================================
# PROBLEM PARAMETERS
# ==========================================
A = 1.0            # Cross-sectional area (in^2)
L = 10.0           # Length of the truss element (in)
P_yield = 60.0     # Yield force = sigma_y * A = 60 ksi * 1 in^2 = 60 kips

# Menegotto-Pinto Material Data
mat = MPMatData(
    E=29000.0, fy=60.0, b=0.05, 
    R0=20.0, cR1=0.925, cR2=0.15, 
    a1=0.0, a2=0.0
)

RELATIVE_TOLERANCE = 1e-4  # Relative Norm Unbalanced Force: 10^-2 to 10^-4
STEP_SIZE = 2.0            # 3.33% of P_yield (3-5% criteria)

# ==========================================
# FUNCTION : Incremental Iterative Solver
# ==========================================
def run_analysis(load_history, method='Full-NR'):
    state = MPState()
    u_current = 0.0         # current displacement
    u_committed = 0.0       # last converged displacement

    # Arrays for plotting
    history_u = [0.0]       # stores the displacement at the end of each successful load increment
    history_P = [0.0]       # stores the force at the end of each successful load increment
    iter_points_u = []      # stores the displacement at the end of each iteration
    iter_points_P = []      # stores the force at the end of each iteration
    all_steps_residuals = {}
    
    K_elastic = (A * mat.E) / L

    number_of_iterations = []

    for step, P_target in enumerate(load_history):     # Incremental Outer Loop 
        if step == 0: continue

        if method == 'Modified-NR':
            MAX_ITER = 2000
        else:
            MAX_ITER = 50
        
        no_of_iteration = 0
        current_step_residuals = []

        for idx in range(MAX_ITER):                    # Iterative Inner Loop
            no_of_iteration += 1
            # --------------------------------------
            # Structural State Determination
            # --------------------------------------
            epsilon = u_current / L                    
            state.eps[0] = epsilon        # Total Strain
            # state.eps[1] = epsilon - state.past.sig / mat.E if step==1 and idx==0 else epsilon - state.past.epr          # Strain increment from last convergence
            state.eps[1] = (u_current - u_committed) / L
            
            state = mate_mp(mat, state)               # Call material 
            
            R = state.sig * A
            psi = P_target - R                        # Unbalanced force residual

            # Initial Unbalanced force residual
            if idx == 0:
                psi_0 = abs(psi)
                if psi_0 < 1e-12:  # Prevent division by zero if step size is zero
                    psi_0 = 1.0

            rel_res = abs(psi) / psi_0   # Relative residual for the current iteration
            current_step_residuals.append(rel_res)
            
            # --------------------------------------
            # Iteration coordinates
            # --------------------------------------
            iter_points_u.append(u_current)
            iter_points_P.append(R)

            # --------------------------------------
            # Convergence Check
            # --------------------------------------
            if rel_res < RELATIVE_TOLERANCE:
                state.commit()
                u_committed = u_current  
                history_u.append(u_current)
                history_P.append(R) 
                number_of_iterations.append(no_of_iteration)
                break
            
            # --------------------------------------
            # Stiffness Update
            # --------------------------------------
            if method == 'Full-NR':
                K_T = (A * state.Et) / L
            elif method == 'Modified-NR':
                K_T = K_elastic
            elif method == 'Hybrid':
                K_T = K_elastic if idx == 0 else (A * state.Et) / L

            # --------------------------------------    
            # Displacement Update
            # --------------------------------------
            du = psi / K_T
            u_current += du
        else:
            print(f"Warning: Reached max iterations at load P = {P_target} using {method}")
        
        all_steps_residuals[step] = current_step_residuals

    print(f" Maximum Number of iterations: {max(number_of_iterations)}")
                    
    return np.array(history_u), np.array(history_P), np.array(iter_points_u), np.array(iter_points_P), all_steps_residuals, number_of_iterations

# ==========================================
# EXECUTION & PLOTTING
# ==========================================

# --------------------------------------------------
# PART A: Monotonic Pushover (Newton Raphson)
# --------------------------------------------------

P_monotonic = np.linspace(0, 200, int(200/STEP_SIZE) + 1)
u_a, P_a, iter_u_a, iter_P_a, history_dict_a, num_iter_a = run_analysis(P_monotonic, method='Full-NR')

MaxIter_a = max(num_iter_a)
TotalIter_a = sum(num_iter_a)

plt.figure(figsize=(10, 6))
plt.plot(u_a, P_a, 'b-', lw=1, label='Pushover Curve (Converged States)')
plt.scatter(iter_u_a, iter_P_a, color='red', s=5, alpha=0.7, label='Iteration Points')
plt.text(.79, 0.07, f'Max. Iterations: {MaxIter_a} \nTotal Iterations: {TotalIter_a}', transform=plt.gca().transAxes)
plt.title(f'Part (a): Monotonic Pushover Analysis (Newton-Raphson)')
plt.xlabel('Displacement, u (in)')
plt.ylabel('External Force, P (kips)')
plt.grid(True, alpha = 0.5, linestyle='--')
plt.legend()
# plt.show()

# --------------------------------------------------
# PART B: Modified Newton-Raphson
# --------------------------------------------------

u_b, P_b, iter_u_b, iter_P_b, history_dict_b, num_iter_b = run_analysis(P_monotonic, method='Modified-NR')

MaxIter_b = max(num_iter_b)
TotalIter_b = sum(num_iter_b)

plt.figure(figsize=(10, 6))
plt.plot(u_b, P_b, 'b-', lw=1, label='Pushover Curve (Converged States)')
plt.scatter(iter_u_b, iter_P_b, color='red', s=5, alpha=0.7, label='Iteration Points')
plt.text(.79, 0.07, f'Max. Iterations: {MaxIter_b} \nTotal Iterations: {TotalIter_b}', transform=plt.gca().transAxes)
plt.title(f'Part (b): Monotonic Pushover Analysis (Modified Newton-Raphson)')
plt.xlabel('Displacement, u (in)')
plt.ylabel('External Force, P (kips)')
plt.grid(True, alpha = 0.5, linestyle='--')
plt.legend()
# plt.show()

# --------------------------------------------------
# PART B: Convergence Comparison
# --------------------------------------------------
yield_step_idx = 31      # Step 31 corresponds to P = 62 kips (just after 60 kips yield point)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle(f'Part (b): Convergence comparison just after yield point')

res_a_yield = history_dict_a[yield_step_idx]
res_b_yield = history_dict_b[yield_step_idx]

# Left plot: Full Newton-Raphson
ax1.semilogy(range(1, len(res_a_yield) + 1), res_a_yield, 'b-o', lw=2, markersize=6)
ax1.axhline(RELATIVE_TOLERANCE, color='black', label=f'Relative Tolerance ({RELATIVE_TOLERANCE})')
ax1.set_title('Full Newton-Raphson')
ax1.set_xlabel('Iteration Number')
ax1.set_ylabel(r'Relative Residual Ratio $\left(\frac{|\psi|}{|\psi_0|}\right)$')
ax1.legend()

# Right plot: Modified Newton-Raphson
ax2.semilogy(range(1, len(res_b_yield) + 1), res_b_yield, 'r-o', lw=2, markersize=3)
ax2.axhline(RELATIVE_TOLERANCE, color='black', label=f'Relative Tolerance ({RELATIVE_TOLERANCE})')
ax2.set_title('Modified Newton-Raphson')
ax2.set_xlabel('Iteration Number')
ax2.set_ylabel(r'Relative Residual Ratio $\left(\frac{|\psi|}{|\psi_0|}\right)$')
ax2.legend()

plt.tight_layout()
# plt.show()

# --------------------------------------------------
# PART C: Cyclic Pushover Analysis 
# --------------------------------------------------

cyclic_targets = [0, 100, -100, 200, 0]
P_cyclic = [0.0]
for k in range(len(cyclic_targets)-1):
    start, end = cyclic_targets[k], cyclic_targets[k+1]
    steps = int(np.ceil(abs(end - start) / STEP_SIZE))
    P_cyclic.extend(np.linspace(start, end, steps + 1)[1:])

u_c, P_c, iter_u_c, iter_P_c, res_c, num_iter_c = run_analysis(P_cyclic, method='Hybrid')

MaxIter_c = max(num_iter_c)
TotalIter_c = sum(num_iter_c)

plt.figure(figsize=(10, 6))
plt.plot(u_c, P_c, 'b-', lw=1, label='Cyclic Response')
plt.scatter(iter_u_c, iter_P_c, color='red', s=1, alpha=0.7, label='Iteration Points')
plt.text(.79, 0.07, f'Max. Iterations: {MaxIter_c} \nTotal Iterations: {TotalIter_c}', transform=plt.gca().transAxes)
plt.axhline(0, color='black', linewidth=1, linestyle='--')
plt.axvline(0, color='black', linewidth=1, linestyle='--')
plt.title('Part (c): Quasi-Static Cyclic Pushover Analysis (Newton-Raphson)')
plt.xlabel('Displacement, u (in)')
plt.ylabel('External Force, P (kips)')
plt.grid(True, alpha = 0.5, linestyle='--')
plt.legend()
plt.show()
