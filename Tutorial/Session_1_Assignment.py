"""
================================================================================
OpenSees Tutorial - Session 1: 3D Frame
Date: 23 June 2026
Author: Niraj Yadav

Assignments:
For the 5-story frame, do the following:
1) Modal analysis
2) Static analysis under gravity loads
3) Visualization of mode shapes

Frame Details:

Baywidth_X = 5.0
Baywidth_Z = 6.0
Baywidth_Y = 4.0

No_of_Bays_X = 2
No_of_Bays_Z = 1
No_of_Bays_Y = 5
================================================================================
"""

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import opsvis as opsv
import os

ops.wipe()
ops.model('BasicBuilder', '-ndm', 3, '-ndf', 6)

#----------------------------------------------------------------------------------
# Geometry, Dimensions And Units (m, s, kN) 
# Global axes: X, Z, Y (vertical)
#----------------------------------------------------------------------------------

g = 9.81                                  # Acceleration due to gravity (m/s^2)
DL = 6.5                                  # Dead load (kN/m^2)
LL = 2.5                                  # Live load (kN/m^2)
gamma_concrete = 24.0                     # Unit weight of concrete (kN/m^3)

# Column properties
A = 0.2025                                # Area of column (m^2)
E = 3e7                                   # Young’s Modulus of concrete (kN/m^2)
G = 125e5                                 # Shear modulus of concrete (kN/m^2)
J = 0.007                                 # Torsional moment of inertia of cross section (m^4)
Iy = 0.45**4 / 12                         # Second moment of area about the local y-axis (m^4)
Iz = 0.45**4 / 12                         # Second moment of area about the local z-axis (m^4)

# Beam X properties
A_x = 0.16                                # Area of beam X (m^2)
E_x = 3e7                                 # Young’s Modulus (kN/m^2)
G_x = 125e5                               # Shear modulus (kN/m^2)
J_x = 0.004                               # Torsional moment of inertia of cross section (m^4)
Iy_x = 0.4**4 / 12                        # Second moment of area about the local y-axis (m^4)
Iz_x = 0.4**4 / 12                        # Second moment of area about the local z-axis (m^4)

# Beam Z properties
A_z = 0.18                                # Area of beam Z (m^2)
E_z = 3e7                                 # Young’s Modulus (kN/m^2)
G_z = 125e5                               # Shear modulus (kN/m^2)
J_z = 0.005                               # Torsional moment of inertia of cross section (m^4)
Iy_z = 0.45 * (0.4**3) / 12               # Second moment of area about the local y-axis (m^4)
Iz_z = 0.4 * (0.45**3) / 12               # Second moment of area about the local z-axis (m^4)

# --------------------------------------------------------------------------------
# Nodes, Tag: XYZ , X = 1, 2, 3 | Z = 1, 2 | Y = 0, 1, 2, 3, 4, 5
# --------------------------------------------------------------------------------

# Floor 0 nodes
ops.node(101, 0.0, 0.0, 0.0)
ops.node(201, 5.0, 0.0, 0.0)
ops.node(301, 10.0, 0.0, 0.0)
ops.node(102, 0.0, 0.0, 6.0)
ops.node(202, 5.0, 0.0, 6.0)
ops.node(302, 10.0, 0.0, 6.0)

# Floor 1 nodes
ops.node(111, 0.0, 4.0, 0.0)
ops.node(211, 5.0, 4.0, 0.0)
ops.node(311, 10.0, 4.0, 0.0)
ops.node(112, 0.0, 4.0, 6.0)
ops.node(212, 5.0, 4.0, 6.0)
ops.node(312, 10.0, 4.0, 6.0)

# Floor 2 nodes
ops.node(121, 0.0, 8.0, 0.0)
ops.node(221, 5.0, 8.0, 0.0)
ops.node(321, 10.0, 8.0, 0.0)
ops.node(122, 0.0, 8.0, 6.0)
ops.node(222, 5.0, 8.0, 6.0)
ops.node(322, 10.0, 8.0, 6.0)

# Floor 3 nodes
ops.node(131, 0.0, 12.0, 0.0)
ops.node(231, 5.0, 12.0, 0.0)
ops.node(331, 10.0, 12.0, 0.0)
ops.node(132, 0.0, 12.0, 6.0)
ops.node(232, 5.0, 12.0, 6.0)
ops.node(332, 10.0, 12.0, 6.0)

# Floor 4 nodes
ops.node(141, 0.0, 16.0, 0.0)
ops.node(241, 5.0, 16.0, 0.0)
ops.node(341, 10.0, 16.0, 0.0)
ops.node(142, 0.0, 16.0, 6.0)
ops.node(242, 5.0, 16.0, 6.0)
ops.node(342, 10.0, 16.0, 6.0)

# Floor 5 nodes
ops.node(151, 0.0, 20.0, 0.0)
ops.node(251, 5.0, 20.0, 0.0)
ops.node(351, 10.0, 20.0, 0.0)
ops.node(152, 0.0, 20.0, 6.0)
ops.node(252, 5.0, 20.0, 6.0)
ops.node(352, 10.0, 20.0, 6.0)

#----------------------------------------------------------------------------------
# Geometric Transformations
#----------------------------------------------------------------------------------

Col_TransfTag = 1
Beam_X_TransfTag = 2
Beam_Z_TransfTag = 3

# geomTransf(  'Linear',        transfTag,     *vecxz)
ops.geomTransf('Linear',    Col_TransfTag, *[0, 0, 1])  # Local x-z plane vector for columns
ops.geomTransf('Linear', Beam_X_TransfTag, *[0, 1, 0])  # Local x-z plane vector for beams in X-direction
ops.geomTransf('Linear', Beam_Z_TransfTag, *[0, 1, 0])  # Local x-z plane vector for beams in Z-direction

#----------------------------------------------------------------------------------
# Elements
# Element tag = start_node * 1000 + end_node
#----------------------------------------------------------------------------------

# Column Elements
for y in range(5):          # y-grid lines 0, 1, 2, 3, 4
    for x in [1, 2, 3]:     # x-grid lines 1, 2, 3
        for z in [1, 2]:    # z-grid lines 1, 2

            node_i = 100 * x + 10 * y + z        
            node_j = 100 * x + 10 * (y + 1) + z  
            eleTag = node_i * 1000 + node_j

            ops.element('elasticBeamColumn', eleTag, node_i, node_j, A, E, G, J, Iy, Iz, Col_TransfTag)

# X-Beam Elements
for y in range(1, 6):       # y-grid lines 1, 2, 3, 4, 5
    for x in [1, 2]:        # x-grid lines 1, 2
        for z in [1, 2]:    # z-grid lines 1, 2

            node_i = 100 * x + 10 * y + z
            node_j = 100 * (x + 1) + 10 * y + z  
            eleTag = node_i * 1000 + node_j

            ops.element('elasticBeamColumn', eleTag, node_i, node_j, A_x, E_x, G_x, J_x, Iy_x, Iz_x, Beam_X_TransfTag)

# Z-Beam Elements
for y in range(1, 6):       # y-grid lines 1, 2, 3, 4, 5
    for x in [1, 2, 3]:     # x-grid lines 1, 2, 3
        for z in [1]:       # z-grid line 1

            node_i = 100 * x + 10 * y + z
            node_j = 100 * x + 10 * y + (z + 1)
            eleTag = node_i * 1000 + node_j

            ops.element('elasticBeamColumn', eleTag, node_i, node_j, A_z, E_z, G_z, J_z, Iy_z, Iz_z, Beam_Z_TransfTag)

#----------------------------------------------------------------------------------
# Nodal Masses  
#----------------------------------------------------------------------------------

mass_floor = (DL + 0.25 * LL) / g   # Mass of floor (kilo-kg/m^2)
mass_concrete = gamma_concrete / g  # Mass of concrete (kilo-kg/m^3)

for y in range(1, 6):               # y-grid lines 1, 2, 3, 4, 5
    if y == 5:                      # top floor 
        trib_y = 4.0 / 2.0          # Tributary width in y-direction for top floor
    else:
        trib_y = 4.0                # Tributary width in y-direction for other floors

    for z in [1, 2]:                # z-grid line 1, 2
        trib_z = 6.0 / 2.0          # Tributary width in z-direction

        for x in [1, 2, 3]:         # x-grid lines 1, 2, 3
            
            if x == 1 or x == 3:    # corner nodes
                trib_x = 5.0 / 2.0  # Tributary width in x-direction for corner nodes
            elif x == 2:
                trib_x = 5.0        # Tributary width in x-direction for interior nodes
            
            mass_slab = trib_x * trib_z * mass_floor
            
            mass_col = 0.45*0.45 * trib_y * mass_concrete
            mass_beamX = 0.40*0.40 * trib_x * mass_concrete
            mass_beamZ = 0.40*0.45 * trib_z * mass_concrete
            
            total_node_mass = mass_slab + mass_col + mass_beamX + mass_beamZ
            
            nodeTag = 100 * x + 10 * y + z

            ops.mass(nodeTag, total_node_mass, total_node_mass, total_node_mass, 0.0, 0.0, 0.0)

# for node in ops.getNodeTags():
#     mass = ops.nodeMass(node)
#     print(f"Node {node}: Mass = {mass}")

#----------------------------------------------------------------------------------
# Boundary Conditions, The base nodes are fixed
#----------------------------------------------------------------------------------

ops.fix(101, 1, 1, 1, 1, 1, 1)
ops.fix(201, 1, 1, 1, 1, 1, 1)
ops.fix(301, 1, 1, 1, 1, 1, 1)
ops.fix(102, 1, 1, 1, 1, 1, 1)
ops.fix(202, 1, 1, 1, 1, 1, 1)
ops.fix(302, 1, 1, 1, 1, 1, 1)

#----------------------------------------------------------------------------------
# Nodal Loads  
#----------------------------------------------------------------------------------

ops.timeSeries('Constant', 100)
ops.pattern('Plain', 100, 100)

weight_floor = DL + LL              # Weight of floor (kN/m^2)
weight_concrete = gamma_concrete    # Weight of concrete (kN/m^3)

for y in range(1, 6):               # y-grid lines 1, 2, 3, 4, 5
    if y == 5:                      # top floor 
        trib_y = 4.0 / 2.0          # Tributary width in y-direction for top floor
    else:
        trib_y = 4.0                # Tributary width in y-direction for other floors

    for z in [1, 2]:                # z-grid line 1, 2
        trib_z = 6.0 / 2.0          # Tributary width in z-direction

        for x in [1, 2, 3]:         # x-grid lines 1, 2, 3
            
            if x == 1 or x == 3:    # corner nodes
                trib_x = 5.0 / 2.0  # Tributary width in x-direction for corner nodes
            elif x == 2:
                trib_x = 5.0        # Tributary width in x-direction for interior nodes
            
            weight_slab = trib_x * trib_z * weight_floor
            
            weight_col = 0.45*0.45 * trib_y * weight_concrete
            weight_beamX = 0.40*0.40 * trib_x * weight_concrete
            weight_beamZ = 0.40*0.45 * trib_z * weight_concrete
            
            total_node_weight = weight_slab + weight_col + weight_beamX + weight_beamZ
            
            nodeTag = 100 * x + 10 * y + z

            ops.load(nodeTag, 0.0, -total_node_weight, 0.0, 0.0, 0.0, 0.0)

# ----------------------------------------------------------------------------------
# Recorders  
# ----------------------------------------------------------------------------------

# Output folder
output_folder = os.path.join("Tutorial", "Session_1_Results")
os.makedirs(output_folder, exist_ok=True)

# Node recorders
ops.recorder('Node', '-file', os.path.join(output_folder, 'Acc312.out'), '-node', 312, '-dof', 1, 2, 3, 4, 5, 6, 'accel')
ops.recorder('Node', '-file', os.path.join(output_folder, 'Vel312.out'), '-node', 312, '-dof', 1, 2, 3, 4, 5, 6, 'vel')
ops.recorder('Node', '-file', os.path.join(output_folder, 'Disp312.out'), '-node', 312, '-dof', 1, 2, 3, 4, 5, 6, 'disp')
ops.recorder('Node', '-file', os.path.join(output_folder, 'R101.out'), '-node', 101, '-dof', 1, 2, 3, 4, 5, 6, 'reaction')
ops.recorder('Element', '-file', os.path.join(output_folder, 'F101111.out'), '-ele', 101111, 'force')

# --------------------------------------------------------------------------------
# Modal/Eigenvalue Analysis 
# --------------------------------------------------------------------------------

numModes = 6
lambdas = ops.eigen(numModes)      # Eigenvalues (rad^2/s^2) from the modal analysis

# Export modal properties
ops.modalProperties('-file', f'{output_folder}/modalProperties.out')

time_periods = [(2 * np.pi) / (lam**0.5) for lam in lambdas]         # T = 2π / √λ  
print("Initial Time Periods:", [f"{p:.5f}" for p in time_periods])

# --------------------------------------------------------------------------------
# Gravity Analysis 
# --------------------------------------------------------------------------------

ops.wipeAnalysis()

ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.algorithm('Linear')

NstepG = 1                                       # No of steps to perform gravity analysis
ops.integrator('LoadControl', 1.0 / NstepG)
ops.analysis('Static')

ops.analyze(NstepG)
ops.loadConst('-time', 0.0)                      # Hold gravity loads constant

# Plotting 
def Plotter():
    """Visualizes the 3D model geometry and applied nodal loads using opsvis."""

    opsv.plot_model(node_labels = 0, element_labels = 0)
    plt.title("3D Model")

    opsv.plot_load(sfac= 1.7, node_supports=True)
    plt.title("Applied Nodal Loads (kN)")

    for text in plt.gca().texts:
        try:
            value = float(text.get_text())
            text.set_text(f"{value:.2f}")
        except ValueError:
            pass

    plt.show()

Plotter()

R = 0.0           # Total Vertical Reaction
ops.reactions()   # Get reactions for all nodes
for node in [101,201,301,102,202,302]:
    R = R + ops.nodeReaction(node)[1]

print("Total Vertical Reaction =", R, "kN")

# --------------------------------------------------------------------------------
# Plotting Mode Shapes and Deformed Shape
# --------------------------------------------------------------------------------

def ModeShapesPlot():
    """Plots the deformed shape and all mode shapes for first numModes modes."""

    opsv.plot_defo()
    plt.title("Deformed Shape")

    for i in range(1, numModes+1):
        opsv.plot_mode_shape(i)
        plt.title(f"Mode {i}")

    plt.show()

ModeShapesPlot()
