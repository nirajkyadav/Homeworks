import numpy as np
import matplotlib.pyplot as plt
import time
import os

from mate_mp import MPMatData, MPState, MPHistory, mate_mp

# =============================================================================
# PROBLEM PARAMETERS
# =============================================================================
M = 360_000.0                     # lumped mass [kg]   (360 ton)
L = 0.5                           # truss length [m]
E = 200e9                         # elastic modulus [Pa]
T_0 = 0.15                        # target initial (elastic) period [s]
XI = 0.02                         # damping ratio [-]

omega_0 = 2.0 * np.pi / T_0       # natural circular frequency [rad/s]
K_0 = M * omega_0 ** 2            # required initial stiffness [N/m]
A = (K_0 * L) / E                 # Q1: truss area [m^2]
C = 2.0 * XI * M * omega_0        # Q2: damping coefficient [N.s/m]

print(f"Area = {A*1e6:.3f} mm^2, Damping Coefficient = {C/1e3:.3f} kN.s/m")

# Menegotto-Pinto Material Data
mat = MPMatData(
    E=E, 
    fy=0.2e9,       # Yield stress [Pa]
    b=0.05,
    R0=20.0, 
    cR1=0.925, 
    cR2=0.15, 
    a1=0.0, 
    a2=0.0
)

# =============================================================================
# MATERIAL STATE HELPER
# =============================================================================
def material_response(committed: MPHistory, eps_total: float, eps_incr: float):
    """Return (stress, tangent, trial_history) for a trial total strain."""
    trial = MPState(eps=[eps_total, eps_incr, 0.0, 0.0], past=committed)
    mate_mp(mat, trial)
    return trial.sig, trial.Et, trial.pres

# =============================================================================
# FUNCTION : Implicit Nonlinear Dynamic Analysis
# =============================================================================
def run_implicit_analysis(t_arr, acc_g, alpha, beta, max_iter=50):
    N = len(t_arr)
    dt = t_arr[1] - t_arr[0]
    
    u = np.zeros(N)                    # Displacement
    v = np.zeros(N)                    # Velocity
    a = np.zeros(N)                    # Acceleration
    P = -M * acc_g                     # Earthquake force vector

    committed = MPHistory()                                                           # virgin material
    sig0, Et0, committed = material_response(committed, 0.0, 0.0)
    R = sig0 * A                                                                      # = 0
    K_T = Et0 * A / L                                                                 # = K_0
    a[0] = (P[0] - C * v[0] - R) / M
    
    c1 = (1.0 / (beta * dt**2)) * M + (alpha / (beta * dt)) * C                       # See Notes for calculation
    c2 = (1.0 / (beta * dt)) * M + ((alpha / beta) - 1.0) * C                         # See Notes for calculation
    c3 = ((1.0 / (2.0 * beta)) - 1.0) * M + dt * ((alpha / (2.0 * beta)) - 1.0) * C   # See Notes for calculation
    
    for n in range(N - 1):                                                            # Time Step Loop
        p_eff = P[n+1] + c1 * u[n] + c2 * v[n] + c3 * a[n]                            # Effective External Load (History Terms)
        u_iter = u[n]
        R_iter = R
        K_T_iter = K_T
        trial_hist = committed                                                        # fallback if already converged
        eps = 1e-4
        
        for i in range(max_iter):                                                     # Newton Raphson Iteration Loop
            psi = p_eff - R_iter - c1 * u_iter                                        # unbalanced (residual) force
            K_T_dyn = K_T_iter + c1
            du  = psi / K_T_dyn                                                       # increment  

            energy = abs(psi * du)                                                    # energy increment 
            if i == 0:
                energy_ref = energy                                                   # reference 
            if energy < eps * max(energy_ref, 1e-30):                                 # relative energy test (1e-30 prevents div-by-zero)
                u_iter = u_iter + du                                                  
                break

            u_iter = u_iter + du

            # structural state determination at the trial displacement
            sig, Et, trial_hist = material_response(committed, u_iter / L, (u_iter - u[n]) / L)
            R_iter = sig * A
            K_T_iter = Et * A / L

        u[n + 1] = u_iter
        R = R_iter
        K_T = K_T_iter
        committed = trial_hist                  # commit converged material state
        
        # Update velocity and acceleration
        delta_u = u[n+1] - u[n]
        v[n+1] = ((alpha / (beta * dt)) * delta_u 
                    + (1.0 - (alpha / beta)) * v[n] 
                    + dt * (1.0 - (alpha / (2.0 * beta))) * a[n])
        a[n+1] = (P[n+1] - C * v[n+1] - R) / M
        
    return u, a

# =============================================================================
# FUNCTION: Explicit Nonlinear Dynamic Analysis
# =============================================================================
def run_explicit_analysis(t_arr, acc_g):
    N = len(t_arr)
    dt = t_arr[1] - t_arr[0]
    
    u = np.zeros(N)                      # Displacement
    v = np.zeros(N)                      # Velocity
    a = np.zeros(N)                      # Acceleration
    P = -M * acc_g                       # Earthquake force vector
    R = np.zeros(N)
    
    committed = MPHistory()
    sig0, _, committed = material_response(committed, 0.0, 0.0)
    R[0] = sig0 * A
    a[0] = (P[0] - C * v[0] - R[0]) / M                       # v0 = 0, R0 = 0
    u_minus_1 = u[0] - dt * v[0] + 0.5 * dt**2 * a[0]         # fictitious u[-1]

    c4 = (M / dt**2) + (C / (2.0 * dt))
    
    for n in range(N - 1):                               # Time Step Loop
        u_prev = u[n - 1] if n > 0 else u_minus_1
        p_eff = (P[n] - R[n]
                 + (2.0 * M / dt**2) * u[n]
                 - (M / dt**2 - C / (2.0 * dt)) * u_prev)
        u[n + 1] = p_eff / c4

        # state determination at the new displacement
        sig, _, committed = material_response(committed, u[n + 1] / L, (u[n + 1] - u[n]) / L)
        R[n + 1] = sig * A
        
        # Update velocity and acceleration
        v[n + 1] = (u[n + 1] - u_prev) / (2.0 * dt)
        a[n + 1] = (u[n + 1] - 2.0 * u[n] + u_prev) / dt**2
        
    return u, a

# =============================================================================
# EXECUTION & PLOTTING
# =============================================================================

# --------------------------------------------------
# Ground motion data, Resampled using MATLAB
# --------------------------------------------------
DATA_DIR = '/Users/niraj/Downloads/HE_Training/Homeworks/HW5/Records'
OUTPUT_DIR = DATA_DIR

# For Implicit Analysis
d_0_02 = np.loadtxt(os.path.join(DATA_DIR, 'Northridge Record.txt'), skiprows=1)
t_0_02, acc_0_02 = d_0_02[:, 0], d_0_02[:, 1]               # dt = 0.020 s
 
d_0_01 = np.loadtxt(os.path.join(DATA_DIR, 'Northridge_0_01.txt'), skiprows=1)
t_0_01, acc_0_01 = d_0_01[:, 0], d_0_01[:, 1]               # dt = 0.010 s
 
d_0_005 = np.loadtxt(os.path.join(DATA_DIR, 'Northridge_0_005.txt'), skiprows=1)
t_0_005, acc_0_005 = d_0_005[:, 0], d_0_005[:, 1]           # dt = 0.005 s
 
d_0_0025 = np.loadtxt(os.path.join(DATA_DIR, 'Northridge_0_0025.txt'), skiprows=1)
t_0_0025, acc_0_0025 = d_0_0025[:, 0], d_0_0025[:, 1]       # dt = 0.0025 s
 
# For Explicit Analysis
d_0_00005 = np.loadtxt(os.path.join(DATA_DIR, 'Northridge_0_00005.txt'), skiprows=1)
t_0_00005, acc_0_00005 = d_0_00005[:, 0], d_0_00005[:, 1]   # dt = 0.00005 s

d_0_0001 = np.loadtxt(os.path.join(DATA_DIR, 'Northridge_0_0001.txt'), skiprows=1)
t_0_0001, acc_0_0001 = d_0_0001[:, 0], d_0_0001[:, 1]       # dt = 0.0001 s

d_0_0002 = np.loadtxt(os.path.join(DATA_DIR, 'Northridge_0_0002.txt'), skiprows=1)
t_0_0002, acc_0_0002 = d_0_0002[:, 0], d_0_0002[:, 1]       # dt = 0.0002 s

d_0_05 = np.loadtxt(os.path.join(DATA_DIR, 'Northridge_0_05.txt'), skiprows=1)
t_0_05, acc_0_05 = d_0_05[:, 0], d_0_05[:, 1]               # dt = 0.05 s

d_0_07 = np.loadtxt(os.path.join(DATA_DIR, 'Northridge_0_07.txt'), skiprows=1)
t_0_07, acc_0_07 = d_0_07[:, 0], d_0_07[:, 1]               # dt = 0.07 s

# --------------------------------------------------
# PART 3A: Implicit Analysis with dt = 0.02s
# --------------------------------------------------
u_avg, a_avg = run_implicit_analysis(t_0_02, acc_0_02, alpha=1/2, beta=1/4) # Constant average, Unconditionally stable
u_lin, a_lin = run_implicit_analysis(t_0_02, acc_0_02, alpha=1/2, beta=1/6) # Linear acceleration, Conditionally stable

print("\nProblem 3a, Implicit Nonlinear Dynamic Analysis, dt = 0.02 s")
print(f"constant-avg : peak|u| = {np.abs(u_avg).max()*1e3:7.3f} mm, peak|a| = {np.abs(a_avg).max():7.3f} m/s^2")
print(f"linear-accel : peak|u| = {np.abs(u_lin).max()*1e3:7.3f} mm, peak|a| = {np.abs(a_lin).max():7.3f} m/s^2")
 
fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
ax[0].plot(t_0_02, u_avg * 1e3, 'b-',  lw=1.0, label=r'Constant avg. accel ($\alpha=1/2$, $\beta=1/4$)')
ax[0].plot(t_0_02, u_lin * 1e3, 'r--', lw=1.0, label=r'Linear accel ($\alpha=1/2$, $\beta=1/6$)')
ax[0].set_ylabel('Displacement [mm]'); ax[0].grid(alpha=0.3); ax[0].legend(loc='upper right')
ax[0].set_title(r'Problem 3a, Implicit Nonlinear Dynamic Analysis, $\Delta t = 0.02$ s')
ax[1].plot(t_0_02, a_avg, 'b-',  lw=1.0, label=r'Constant avg. accel ($\alpha=1/2$, $\beta=1/4$)')
ax[1].plot(t_0_02, a_lin, 'r--', lw=1.0, label=r'Linear accel ($\alpha=1/2$, $\beta=1/6$)')
ax[1].set_ylabel('Relative Acceleration [m/s$^2$]'); ax[1].set_xlabel('Time [s]')
ax[1].grid(alpha=0.3); ax[1].legend(loc='upper right')
fig.tight_layout()
fig.savefig("Problem_3a.pdf")
# plt.show()

# --------------------------------------------------
# PART 3B: Implicit Time-step convergence Analysis, constant-average acceleration method
# --------------------------------------------------
print("\nProblem 3b, Implicit Time-step convergence Analysis, constant-average acceleration method")
records_implicit = [(t_0_02,  acc_0_02,  0.020),
                    (t_0_01,  acc_0_01,  0.010),
                    (t_0_005, acc_0_005, 0.005),
                    (t_0_0025, acc_0_0025, 0.0025)
                    ]

colors = ['k', 'b', 'r', 'g']

fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
for (t, ag, dt), col in zip(records_implicit, colors):
    u, a = run_implicit_analysis(t, ag, alpha=0.5, beta=0.25)
    print(f"dt = {dt} s, peak|u| = {np.abs(u).max()*1e3:7.3f} mm, peak|a| = {np.abs(a).max():7.3f} m/s^2")
    ax[0].plot(t, u * 1e3, col + '-', lw=1.0, label=f'$\\Delta t$ = {dt:g} s')
    ax[1].plot(t, a,       col + '-', lw=1.0, label=f'$\\Delta t$ = {dt:g} s')

ax[0].set_ylabel('Displacement [mm]'); ax[0].grid(alpha=0.3); ax[0].legend(loc='upper right')
ax[0].set_title('Problem 3b, Implicit Time-step convergence Analysis using constant-average acceleration method')
ax[1].set_ylabel('Relative Acceleration [m/s$^2$]'); ax[1].set_xlabel('Time [s]')
ax[1].grid(alpha=0.3); ax[1].legend(loc='upper right')
fig.tight_layout()
fig.savefig("Problem_3b.pdf")
# plt.show()

# --------------------------------------------------
# PART 4A: Explicit Analysis based on CFL condition
# --------------------------------------------------
print("\nProblem 4a, Explicit Nonlinear Dynamic Analysis")
rho = 7850.0                     # Mass Density of Steel (kg/m^3)
c_wave  = np.sqrt(E / rho)       # wave speed (m/s)
dt_critical_cfl  = L / c_wave    # wave-propagation CFL (s)

print(f"Critical Time Step based on CFL condition = {dt_critical_cfl:.6f} s")

u_ex, a_ex = run_explicit_analysis(t_0_00005, acc_0_00005)
print(f"Explicit Analysis, peak|u| = {np.abs(u_ex).max()*1e3:7.3f} mm, peak|a| = {np.abs(a_ex).max():7.3f} m/s^2")

fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
ax[0].plot(t_0_00005, u_ex * 1e3, 'g-', lw=1.0); ax[0].set_ylabel('Displacement [mm]')
ax[0].grid(alpha=0.3); ax[0].set_title(r'Problem 4a, Explicit Nonlinear Dynamic Analysis, $\Delta t = 0.00005$ s')
ax[1].plot(t_0_00005, a_ex, 'g-', lw=1.0); ax[1].set_ylabel('Relative Acceleration [m/s$^2$]')
ax[1].set_xlabel('Time [s]'); ax[1].grid(alpha=0.3)
fig.tight_layout()
fig.savefig("Problem_4a.pdf")
# plt.show()

# --------------------------------------------------
# PART 4B: Explicit Time-step convergence Analysis
# --------------------------------------------------
print("\nProblem 4b, Explicit Time-step convergence Analysis")
records_explicit = [(t_0_00005,  acc_0_00005,  0.00005),
                    (t_0_0001,  acc_0_0001,  0.0001),
                    (t_0_0002, acc_0_0002, 0.0002),
                    (t_0_0025, acc_0_0025, 0.0025),
                    (t_0_02, acc_0_02, 0.02),
                    (t_0_05, acc_0_05, 0.05),
                    (t_0_07, acc_0_07, 0.07)
                    ]

colors = ['k', 'b', 'r', 'g','y','m','c']

fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
for (t, ag, dt), col in zip(records_explicit, colors):
    u, a = run_explicit_analysis(t, ag)
    print(f"dt = {dt} s, peak|u| = {np.abs(u).max()*1e3:7.3f} mm, peak|a| = {np.abs(a).max():7.3f} m/s^2")
    ax[0].plot(t, u * 1e3, col + '-', lw=1.0, label=f'$\\Delta t$ = {dt:g} s')
    ax[1].plot(t, a,       col + '-', lw=1.0, label=f'$\\Delta t$ = {dt:g} s')
 
ax[0].set_ylabel('Displacement [mm]'); ax[0].grid(alpha=0.3); ax[0].legend(loc='upper right')
ax[0].set_title('Problem 4b, Explicit Time-step convergence Analysis')
ax[1].set_ylabel('Relative Acceleration [m/s$^2$]'); ax[1].set_xlabel('Time [s]')
ax[1].grid(alpha=0.3); ax[1].legend(loc='upper right')
fig.tight_layout()
fig.savefig("Problem_4b.pdf")
# plt.show()

# --------------------------------------------------
# PART 5: Comparison of most efficient implicit and explicit method with total computation time
# --------------------------------------------------
print("\nProblem 5, Comparison of most efficient implicit and explicit method with total computation time")

# Most efficient implicit, let's take dt = 0.005s
t0 = time.perf_counter()
u_imp, a_imp = run_implicit_analysis(t_0_005, acc_0_005, alpha=0.5, beta=0.25)
time_imp = time.perf_counter() - t0

# Most efficient explicit, let's take dt = 0.0025s
t1 = time.perf_counter()
u_exp, a_exp = run_explicit_analysis(t_0_0025, acc_0_0025)
time_exp = time.perf_counter() - t1

print(f"Implicit (const-avg, dt=0.005): peak |u|={np.abs(u_imp).max()*1e3:7.3f} mm, peak |a|={np.abs(a_imp).max():6.3f} m/s^2, steps={len(t_0_005)}, Total computation time = {time_imp*1e3:7.1f} ms")
print(f"Explicit (central,   dt=0.0025): peak |u|={np.abs(u_exp).max()*1e3:7.3f} mm, peak |a|={np.abs(a_exp).max():6.3f} m/s^2, steps={len(t_0_0025)}, Total computation time = {time_exp*1e3:7.1f} ms")
print(f"speed ratio (implicit / explicit time) = {time_imp/time_exp:.2f}x")

fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
ax[0].plot(t_0_005, u_imp * 1e3, 'b-',  lw=1.0, label=r'Implicit const-avg, $\Delta t$=0.005 s')
ax[0].plot(t_0_0025, u_exp * 1e3, 'g--', lw=1.0, label=r'Explicit central, $\Delta t$=0.0025 s')
ax[0].set_ylabel('Displacement [mm]'); ax[0].grid(alpha=0.3)
ax[0].legend(loc='upper right'); ax[0].set_title('Problem 5, Comparison of most efficient implicit and explicit method')
ax[1].plot(t_0_005, a_imp, 'b-',  lw=1.0, label=r'Implicit const-avg, $\Delta t$=0.005 s')
ax[1].plot(t_0_0025, a_exp, 'g--', lw=1.0, label=r'Explicit central, $\Delta t$=0.0025 s' )
ax[1].set_ylabel('Relative Acceleration [m/s$^2$]'); ax[1].set_xlabel('Time [s]')
ax[1].grid(alpha=0.3); ax[1].legend(loc='upper right')
fig.tight_layout()
fig.savefig("Problem_5.pdf")
plt.show()
