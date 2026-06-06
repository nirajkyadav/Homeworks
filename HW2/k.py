import numpy as np

def get_float(prompt):
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("  Invalid input. Please enter a number.")

def print_matrix(name, M):
    print(f"\n{name}")
    print("-" * (11 * M.shape[1]))
    for row in M:
        print("  ".join(f"{v:+.4e}" for v in row))
    print("-" * (11 * M.shape[1]))

print("=" * 50)
print("  Beam Element Stiffness Matrix Calculator")
print("  Units: kip, in, ksi")
print("=" * 50)

E     = 29000.0        # ksi
A     = 11.8           # in^2
I     = 210.0          # in^4
L     = 20.0 * 12.0    # in
theta = 0.0            # degrees

theta_rad = np.radians(theta)
C = np.cos(theta_rad)
S = np.sin(theta_rad)

# --- Local stiffness matrix [k'] ---
f = E * I / L**3
r = A * L**2 / I

kp = f * np.array([
    [ r,    0,      0,     -r,    0,      0    ],
    [ 0,    12,     6*L,    0,   -12,     6*L  ],
    [ 0,    6*L,   4*L**2,  0,   -6*L,   2*L**2],
    [-r,    0,      0,      r,    0,      0    ],
    [ 0,   -12,    -6*L,    0,    12,    -6*L  ],
    [ 0,    6*L,   2*L**2,  0,   -6*L,   4*L**2]
])

# --- Transformation matrix [T] ---
T = np.array([
    [ C,  S,  0,  0,  0,  0],
    [-S,  C,  0,  0,  0,  0],
    [ 0,  0,  1,  0,  0,  0],
    [ 0,  0,  0,  C,  S,  0],
    [ 0,  0,  0, -S,  C,  0],
    [ 0,  0,  0,  0,  0,  1]
])

# --- Global stiffness matrix [k] = [T]^T [k'] [T] ---
kg = T.T @ kp @ T

print("\n" + "=" * 50)
print(f"  EI/L³ = {f:.4e} [kip/in]    AL²/I = {r:.4f}")
print("=" * 50)

print_matrix("[k'] — Local stiffness matrix (6×6)", kp)
print_matrix("[T]  — Transformation matrix (6×6)", T)
print_matrix("[k]  — Global stiffness matrix (6×6)", kg)