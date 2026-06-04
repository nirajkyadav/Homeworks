"""
Problem #3 + BONUS - Timoshenko cantilever
2-node LINEAR element (Problem #1)  vs  3-node QUADRATIC element (Problem #2)

Cantilever:  L=24 in, rect 3x6 in, E=29000 ksi, nu=0.3, tip load F=10 kips up.
Meshes: 2, 4, 8 uniform elements.   DOF/node: [w, alpha].
"""
import numpy as np
import matplotlib.pyplot as plt

# ===========================================================================
# Material & Section Properties
# ===========================================================================
E = 29000.0             # Modulus of elasticity (ksi)
nu = 0.3                # Poisson's ratio
G = E / (2 * (1 + nu))  # Shear modulus (ksi)

b = 3.0                 # width (in)
h = 6.0                 # depth (in)
L = 24.0                # length (in)

A = b * h               # Area (in^2)
I = (b * h**3) / 12.0   # Moment of Inertia (in^4)

ks = 5.0 / 6.0          # Shear correction factor (rectangular)
A_s = ks * A            # Shear Area (in^2)
GAs = G * A_s           # Shear Rigidity (kips)

F_tip = 10.0            # kips (up, at tip)

# ===========================================================================
# Element Stiffness Matrices
# ===========================================================================
def element_stiffness(le):                       # Problem #1, 2-node, [w1,a1,w2,a2]
    phi = 12.0*E*I/(GAs*le**2)
    c   = GAs*le/12.0
    Ke  = c*np.array([
        [ 12.0/le**2,  6.0/le,   -12.0/le**2,  6.0/le ],
        [  6.0/le,     4.0+phi,   -6.0/le,      2.0-phi],
        [-12.0/le**2, -6.0/le,    12.0/le**2,  -6.0/le ],
        [  6.0/le,     2.0-phi,   -6.0/le,      4.0+phi]])
    return Ke, phi

def Ke_quadratic(le):                            # Problem #2, 3-node, [w1,a1,w2,a2,w3,a3]
    g = GAs; eI = E*I
    return np.array([
     [7*g/(3*le),  g/2,                 -8*g/(3*le),  2*g/3,                g/(3*le),   -g/6],
     [g/2,         7*eI/(3*le)+2*g*le/15,-2*g/3,      -8*eI/(3*le)+g*le/15, g/6,         eI/(3*le)-g*le/30],
     [-8*g/(3*le), -2*g/3,               16*g/(3*le),  0.0,                -8*g/(3*le),  2*g/3],
     [2*g/3,       -8*eI/(3*le)+g*le/15, 0.0,          16*eI/(3*le)+8*g*le/15,-2*g/3,   -8*eI/(3*le)+g*le/15],
     [g/(3*le),    g/6,                 -8*g/(3*le), -2*g/3,                7*g/(3*le), -g/2],
     [-g/6,        eI/(3*le)-g*le/30,    2*g/3,       -8*eI/(3*le)+g*le/15, -g/2,        7*eI/(3*le)+2*g*le/15]])

# ===========================================================================
# LINEAR solver (2-node) - unchanged
# ===========================================================================
def solve(n_elements):
    le = L/n_elements
    n_nodes = n_elements+1
    ndof = 2*n_nodes

    K = np.zeros((ndof, ndof))
    Ke, phi = element_stiffness(le)
    Force_vector = np.zeros(ndof)
    Displacement_vector = np.zeros(ndof)
    Force_vector[2*(n_nodes-1)] = F_tip

    # (b) assemble
    for e in range(n_elements):
        dofs = [2*e, 2*e+1, 2*e+2, 2*e+3]
        for i in range(4):
            for j in range(4):
                K[dofs[i], dofs[j]] += Ke[i, j]

    # (c) solve
    free_dofs = np.arange(2, ndof)
    Displacement_vector[free_dofs] = np.linalg.solve(
        K[np.ix_(free_dofs, free_dofs)], Force_vector[free_dofs])

    x = np.linspace(0, L, n_nodes)
    w = Displacement_vector[0::2]
    alpha = Displacement_vector[1::2]

    # (d) element state determination -> BMD
    xm, M = [], []
    for e in range(n_elements):
        de = np.array([w[e], alpha[e], w[e+1], alpha[e+1]])
        V_i, M_i, V_j, M_j = Ke @ de
        xm += [x[e],  x[e+1]]
        M  += [-M_i,  M_j]
    return dict(le=le, phi=phi, x=x, w=w, alpha=alpha, xm=np.array(xm), M=np.array(M))

# ===========================================================================
# QUADRATIC solver (3-node) - meshed exactly like the linear one
#   n_elements quadratic elements -> n_nodes = 2*n_elements + 1 (shared ends)
#   element e uses nodes [2e, 2e+1, 2e+2]
# ===========================================================================
def solve_quad(n_elements):
    le = L/n_elements                    # full length of each 3-node element
    n_nodes = 2*n_elements + 1
    ndof = 2*n_nodes

    K = np.zeros((ndof, ndof))
    Ke = Ke_quadratic(le)
    Force_vector = np.zeros(ndof)
    Displacement_vector = np.zeros(ndof)
    Force_vector[2*(n_nodes-1)] = F_tip

    # (b) assemble
    for e in range(n_elements):
        n0 = 2*e                          # first global node of this element
        dofs = [2*n0, 2*n0+1, 2*n0+2, 2*n0+3, 2*n0+4, 2*n0+5]
        for i in range(6):
            for j in range(6):
                K[dofs[i], dofs[j]] += Ke[i, j]

    # (c) solve
    free_dofs = np.arange(2, ndof)
    Displacement_vector[free_dofs] = np.linalg.solve(
        K[np.ix_(free_dofs, free_dofs)], Force_vector[free_dofs])

    x = np.linspace(0, L, n_nodes)
    w = Displacement_vector[0::2]
    alpha = Displacement_vector[1::2]

    # (d) element state determination -> BMD  (use END-node moments: -left, +right)
    xm, M = [], []
    for e in range(n_elements):
        n0 = 2*e
        dofs = [2*n0, 2*n0+1, 2*n0+2, 2*n0+3, 2*n0+4, 2*n0+5]
        fe = Ke @ Displacement_vector[dofs]      # [V1,M1,V2,M2,V3,M3]
        xm += [x[2*e], x[2*e+2]]                 # left + right end nodes
        M  += [-fe[1], fe[5]]                    # -left, +right
    return dict(le=le, x=x, w=w, alpha=alpha, xm=np.array(xm), M=np.array(M))

# ===========================================================================
# Analytical Solution
# ===========================================================================
xx = np.linspace(0, L, 100)
Deflection_x = F_tip*xx**2*(3*L-xx)/(6*E*I) + F_tip*xx/GAs
Analytical_Tip_Deflection = F_tip*L**3/(3*E*I) + F_tip*L/GAs

# ===========================================================================
# Quadratic shape-function interpolation (your notes: xi in [-1,1], nodes -1,0,1)
#   N1 = (xi/2)(xi-1),  N2 = 1 - xi^2,  N3 = (xi/2)(xi+1)
#   w(xi) = N1*w1 + N2*w2 + N3*w3  -> evaluated across each 3-node element
# ===========================================================================
def quad_interp(x_nodes, vals, npts=40):
    xs, ys = [], []
    n_elem = (len(x_nodes) - 1) // 2
    for e in range(n_elem):
        i = 2*e                                  # left, mid, right node indices
        xL, xM, xR = x_nodes[i], x_nodes[i+1], x_nodes[i+2]
        half = (xR - xL) / 2.0                   # x = xM + xi*(le/2)
        xi = np.linspace(-1, 1, npts)
        N1 = 0.5*xi*(xi - 1)                     # node @ xi = -1
        N2 = 1 - xi**2                           # node @ xi =  0  (middle)
        N3 = 0.5*xi*(xi + 1)                     # node @ xi = +1
        xs.append(xM + xi*half)
        ys.append(N1*vals[i] + N2*vals[i+1] + N3*vals[i+2])
    return np.concatenate(xs), np.concatenate(ys)

# ===========================================================================
# Meshes - linear (2,4,8) + single quadratic element
# ===========================================================================
results      = {n: solve(n)      for n in (2, 4, 8)}
results_quad = {n: solve_quad(n) for n in (1,)}

print("-----------------------------------------------------------------------------")
print(f"G={G:.2f} ksi  I={I} in^4  As={A_s} in^2  GAs={GAs:.1f} kips")
print(f"Analytical tip w: Timoshenko={Analytical_Tip_Deflection:.6f} in\n")
for n in (2, 4, 8):
    rl = results[n]
    print(f"{n} elems | LINEAR tip={rl['w'][-1]:.6f} in ({100*rl['w'][-1]/Analytical_Tip_Deflection:5.1f}%) "
          f"M_fix={rl['M'][0]:7.2f} ")

print("-----------------------------------------------------------------------------")
for n in (1,):
    rq = results_quad[n]
    print(f"{n} elem  | QUAD   tip={rq['w'][-1]:.6f} in "
          f"({100*rq['w'][-1]/Analytical_Tip_Deflection:5.1f}%) M_fix={rq['M'][0]:7.2f}")

# ===========================================================================
# Plots
# ===========================================================================
colors = {2:"#1f77b4", 4:"#ff7f0e", 8:"#2ca02c"}
quad_color = "#d62728"   # separate color for the quadratic curve

# Deflection  (quadratic drawn as a smooth parabola via the shape functions)
plt.figure(figsize=(9,5.5))
for n in (2, 4, 8):
    plt.plot(results[n]["x"], results[n]["w"], "-o", color=colors[n], label=f"Linear, {n} elem")
for n in (1,):
    rq = results_quad[n]
    xs, ws = quad_interp(rq["x"], rq["w"])                       # smooth parabolic field
    plt.plot(xs, ws, "--", color=quad_color, lw=1.8, label=f"Quad, {n} elem")
    plt.plot(rq["x"], rq["w"], "s", color=quad_color, mfc="white")   # mark the 3 nodes
plt.plot(xx, Deflection_x, "k-", lw=2, label="Analytical (Timoshenko)")
plt.xlabel("x (in)"); plt.ylabel("deflection w (in)")
plt.title("Vertical deflection - linear vs quadratic")
plt.xticks(np.arange(0, L+1, 3)); plt.legend(fontsize=8, ncol=2)
plt.grid(alpha=0.3); plt.tight_layout()
plt.savefig("deflection_bonus_question.pdf")

# BMD  (M = EI*alpha' is LINEAR within a quad element -> straight line is exact)
plt.figure(figsize=(9,5.5))
for n in (2, 4, 8):
    plt.plot(results[n]["xm"], results[n]["M"], "-o", color=colors[n], label=f"Linear, {n} elem")
for n in (1,):
    plt.plot(results_quad[n]["xm"], results_quad[n]["M"], "--s", color=quad_color,
             mfc="white", label=f"Quad, {n} elem")
plt.plot([0, L], [F_tip*L, 0], "k-", lw=2, label="Analytical  M=F(L-x)")
plt.xlabel("x (in)"); plt.ylabel("bending moment (kip-in)")
plt.title("Bending Moment Diagram - linear vs quadratic")
plt.xticks(np.arange(0, L+1, 3)); plt.gca().invert_yaxis()
plt.legend(fontsize=8, ncol=2); plt.grid(alpha=0.3); plt.tight_layout()
plt.savefig("bmd_bonus_question.pdf")

plt.show()