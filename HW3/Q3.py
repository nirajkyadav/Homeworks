"""
Problem #3 Solution

Consider a cantilever beam with a total length of 24 in and a rectangular cross section 
as shown. We want to use MATLAB or Python to assemble the stiffness matrix for this 
structure using the Timoshenko element derived in Problem #1, solve the equations, 
and find the deformation of the beam and the bending moment diagram (BMD). 
Follow these steps: 

a. Consider three uniform mesh discretizations: 2 elements, 4 elements, and 8 
elements.  

b. Write a simple code to assemble the structural stiffness matrix based on the 
element stiffness matrix derived in Problem #1 (no need to make the code 
sophisticated or fancy). 

c. Solve the equilibrium equation to find the nodal displacement vectors. Plot the 
vertical deflection of the beam. 

d. Use element state determination to find the nodal forces. Use the nodal moments 
to plot the BMD of the beam.

e. Compare the derived vertical deflection and BMD for different mesh 
discretizations (plot them clearly on top of each other to be able to compare).  

f. Compare the results with analytical solution - use the analytically-driven solution 
for cantilever beam from the notes or the analytically-driven stiffness matrix for 
linear-elastic Timoshenko beam element to solve the structure. State your 
observations with reasons.  

Cantilever:  L = 24 in,  rectangular section b = 3 in x h = 6 in
Material:    E = 29000 ksi,  nu = 0.3
Load:        F = 10 kips, vertical (up), at the free tip
Meshes:      2, 4, 8 uniform elements
DOF / node:  [ w (transverse), alpha (section rotation) ] 

Element stiffness from Problem #1:

        K_e = (GAs*le/12) * [[ 12/le^2,  6/le,  -12/le^2,  6/le ],
                             [  6/le,   4+phi,  -6/le,    2-phi],
                             [-12/le^2, -6/le,  12/le^2, -6/le ],
                             [  6/le,   2-phi,  -6/le,    4+phi]]
        phi = 12*E*I/(GAs*le^2)
"""

import numpy as np
import matplotlib.pyplot as plt

# ===========================================================================
# Material & Section Properties
# ===========================================================================

# Units: inches (in) for length, kips for force, and KSI for modulus of elasticity

E = 29000.0             # Modulus of elasticity (ksi)
nu = 0.3                # Poisson's ratio
G = E / (2 * (1 + nu))  # Shear modulus (ksi)

b = 3.0                 # width (in)
h = 6.0                 # depth (in)
L = 24.0                # length (in)

A = b * h               # Area (in^2)
I = (b * h**3) / 12.0   # Moment of Inertia (in^4)

ks = 5.0 / 6.0          # Shear correction factor for rectangular section
A_s = ks * A            # Shear Area (in^2)
GAs = G * A_s           # Shear Rigidity (kips)

F_tip = 10.0                # kips (up, at x = L or tip)           

# ===========================================================================
# Element Stiffness Matrix (from Problem #1)
# ===========================================================================

def element_stiffness(le):
    phi = 12.0*E*I/(GAs*le**2)
    c   = GAs*le/12.0
    Ke  = c*np.array([
        [ 12.0/le**2,  6.0/le,   -12.0/le**2,  6.0/le ],
        [  6.0/le,     4.0+phi,   -6.0/le,      2.0-phi],
        [-12.0/le**2, -6.0/le,    12.0/le**2,  -6.0/le ],
        [  6.0/le,     2.0-phi,   -6.0/le,      4.0+phi]])
    return Ke, phi

# ===========================================================================
# Assemble Global Stiffness Matrix, Apply Boundary Conditions, Solve for Displacements and Forces
# ===========================================================================

def solve(n_elements):

    le = L/n_elements                       # length of each element
    n_nodes = n_elements+1                  # number of nodes
    ndof = 2*n_nodes                        # total number of DOFs

    K = np.zeros((ndof, ndof))              # Global Stiffeness Matrix
    Ke, phi = element_stiffness(le)         # Individual Element Stiffness Matrix

    Force_vector = np.zeros(ndof)           # Force Vector
    Displacement_vector = np.zeros(ndof)    # Displacement Vector

    Force_vector[2*(n_nodes-1)] = F_tip     # tip vertical load

    # ----------------------------------------------
    # (b) Assemble Global Stiffness Matrix
    # ----------------------------------------------

    for e in range(n_elements):                     
        dofs = [2*e, 2*e+1, 2*e+2, 2*e+3]   # global DOFs for the current element
        
        # K[np.ix_(dofs, dofs)] += Ke
        for i in range(4):
            for j in range(4):
                K[dofs[i], dofs[j]] += Ke[i, j]

    # ----------------------------------------------
    # (c) Solving [K_reduced]{d_free} = {F_free} to find Displacement Vector
    # ----------------------------------------------
    # Since the leftmost Node (Node 0) is fixed, w0 (DOF 0) and alpha0 (DOF 1) = 0 ; other nodes are free to have displacements
    free_dofs = np.arange(2, ndof)  
    K_reduced = K[np.ix_(free_dofs, free_dofs)]
    F_reduced = Force_vector[free_dofs]

    Displacement_vector[free_dofs] = np.linalg.solve(K_reduced, F_reduced)

    # ----------------------------------------------
    # (d) Element State Determination, Use {F} = [K]{d} to find the Force Vector
    # ----------------------------------------------

    x = np.linspace(0, L, n_nodes)         # x-coordinate of each node
    xm = []                                # x-positions for moment diagram
    M  = []                                # moment values

    w = Displacement_vector[0::2]          # w[0], w[1], ..., w[n_nodes-1]
    alpha = Displacement_vector[1::2]      # alpha[0], alpha[1], ..., alpha[n_nodes-1]

    for e in range(n_elements):
        node_i = e
        node_j = e + 1

        de = np.array([w[node_i], alpha[node_i], w[node_j], alpha[node_j]])     # displacement vector for this element

        V_i, M_i, V_j, M_j = Ke @ de                                            # element forces [V_i, M_i, V_j, M_j]

        xm.append(x[node_i])
        M.append(-M_i)                                                          # moment at start of element

        xm.append(x[node_j])
        M.append(M_j)                                                           # moment at end of element

    return dict(le=le, phi=phi, x=x, w=w, alpha=alpha, xm=np.array(xm), M=np.array(M))

# ===========================================================================
# Analytical Solution
# ===========================================================================

xx = np.linspace(0, L, 100)
Deflection_x = F_tip*xx**2*(3*L-xx)/(6*E*I)+F_tip*xx/GAs

Analytical_Tip_Deflection = F_tip*L**3/(3*E*I) + F_tip*L/GAs

# ===========================================================================
# Approximation using Meshes of size 2, 4, and 8
# ===========================================================================

results = {n: solve(n) for n in (2, 4, 8)}

print("-----------------------------------------------------------------------------")
print(f"G={G:.2f} ksi  I={I} in^4  As={A_s} in^2  GAs={GAs:.1f} kips")
print(f"Analytical tip w: Timoshenko={Analytical_Tip_Deflection:.6f} in\n")
for n, r in results.items():
    print(f"{n} elems: le={r['le']:5.2f}  phi_e={r['phi']:7.4f}  "
          f"tip w={r['w'][-1]:.6f} in ({100*r['w'][-1]/Analytical_Tip_Deflection:5.1f}% of Timoshenko)  "
          f"M_fixed={r['M'][0]:7.2f} kip-in")

# ===========================================================================
# Plot Results
# ===========================================================================

colors = {2:"#1f77b4", 4:"#ff7f0e", 8:"#2ca02c"}

# ----------------------------------------------
# Deflection
# ----------------------------------------------

plt.figure(figsize=(8,5))                         
for n, r in results.items():
    plt.plot(r["x"], r["w"], "-o", color=colors[n], label=f"{n} elements")
plt.plot(xx, Deflection_x, "k--", lw=1.5, label="Analytical (Timoshenko)")
plt.xlabel("x (in)"); plt.ylabel("deflection w (in)")
plt.title("Vertical deflection - mesh comparison")
plt.xticks(np.arange(0, L+1, 3))
plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
plt.savefig("deflection.pdf")

# ----------------------------------------------
# Bending Moment Diagram
# ----------------------------------------------

plt.figure(figsize=(8,5))                         
for n, r in results.items():
    plt.plot(r["xm"], r["M"], "-o", color=colors[n], label=f"{n} elements")
plt.plot([0, L], [F_tip*L, 0], "k--", lw=1.5, label="Analytical  M=F(L-x)")
plt.xlabel("x (in)"); plt.ylabel("bending moment (kip-in)")
plt.title("Bending Moment Diagram - mesh comparison")
plt.xticks(np.arange(0, L+1, 3))
plt.gca().invert_yaxis(); plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
plt.savefig("bmd.pdf")

plt.show()