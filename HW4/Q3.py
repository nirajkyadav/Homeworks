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

TOLERANCE = 1e-8
STEP_SIZE = 2.0       # 3.33% of P_yield (3-5% criteria)

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
    residuals_per_iteration = []
    
    K_elastic = (A * mat.E) / L

    number_of_iterations = []

    for step, P_target in enumerate(load_history):     # Incremental Outer Loop 
        if step == 0: continue

        if method == 'Modified-NR':
            MAX_ITER = 2000
        else:
            MAX_ITER = 50
        
        no_of_iteration = 0

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
            
            # --------------------------------------
            # Iteration coordinates
            # --------------------------------------
            iter_points_u.append(u_current)
            iter_points_P.append(R)
            residuals_per_iteration.append(abs(psi))
            
            # --------------------------------------
            # Convergence Check
            # --------------------------------------
            if abs(psi) < TOLERANCE:
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

    print(f" Maximum Number of iterations: {max(number_of_iterations)}")
                    
    return np.array(history_u), np.array(history_P), np.array(iter_points_u), np.array(iter_points_P), residuals_per_iteration, max(number_of_iterations)

# ==========================================
# EXECUTION & PLOTTING
# ==========================================

# --------------------------------------------------
# PART A: Monotonic Pushover (Newton Raphson)
# --------------------------------------------------

P_monotonic = np.linspace(0, 200, int(200/STEP_SIZE) + 1)
u_a, P_a, iter_u_a, iter_P_a, res_a, MaxIter_a = run_analysis(P_monotonic, method='Full-NR')

plt.figure(figsize=(10, 6))
plt.plot(u_a, P_a, 'b-', lw=1, label='Pushover Curve (Converged States)')
plt.scatter(iter_u_a, iter_P_a, color='red', s=5, alpha=0.7, label='Iteration Points')
plt.text(.79, 0.07, f'Max. Iterations: {MaxIter_a} \nTotal Iterations: {len(res_a)}', transform=plt.gca().transAxes)
plt.title(f'Part (a): Monotonic Pushover Analysis (Newton-Raphson)')
plt.xlabel('Displacement, u (in)')
plt.ylabel('External Force, P (kips)')
plt.grid(True, alpha = 0.5, linestyle='--')
plt.legend()
# plt.show()

# --------------------------------------------------
# PART B: Modified Newton-Raphson & Convergence Comparison
# --------------------------------------------------

u_b, P_b, iter_u_b, iter_P_b, res_b, MaxIter_b = run_analysis(P_monotonic, method='Modified-NR')

plt.figure(figsize=(10, 6))
plt.plot(u_b, P_b, 'b-', lw=1, label='Pushover Curve (Converged States)')
plt.scatter(iter_u_b, iter_P_b, color='red', s=5, alpha=0.7, label='Iteration Points')
plt.text(.79, 0.07, f'Max. Iterations: {MaxIter_b} \nTotal Iterations: {len(res_b)}', transform=plt.gca().transAxes)
plt.title(f'Part (b): Monotonic Pushover Analysis (Modified Newton-Raphson)')
plt.xlabel('Displacement, u (in)')
plt.ylabel('External Force, P (kips)')
plt.grid(True, alpha = 0.5, linestyle='--')
plt.legend()
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

u_c, P_c, iter_u_c, iter_P_c, res_c, MaxIter_c = run_analysis(P_cyclic, method='Hybrid')

plt.figure(figsize=(10, 6))
plt.plot(u_c, P_c, 'b-', lw=1, label='Cyclic Response')
plt.scatter(iter_u_c, iter_P_c, color='red', s=1, alpha=0.7, label='Iteration Points')
plt.text(.79, 0.07, f'Max. Iterations: {MaxIter_c} \nTotal Iterations: {len(res_c)}', transform=plt.gca().transAxes)
plt.axhline(0, color='black', linewidth=1, linestyle='--')
plt.axvline(0, color='black', linewidth=1, linestyle='--')
plt.title('Part (c): Quasi-Static Cyclic Pushover Analysis (Newton-Raphson)')
plt.xlabel('Displacement, u (in)')
plt.ylabel('External Force, P (kips)')
plt.grid(True, alpha = 0.5, linestyle='--')
plt.legend()
plt.show()
