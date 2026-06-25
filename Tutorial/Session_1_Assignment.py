import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import opsvis as opsv
import os

ops.wipe()
ops.model('BasicBuilder', '-ndm', 3, '-ndf', 6)

#----------------------------------------------------------------------------------
# Geometry, Dimensions And Units (m, s, kN) , Global axes X, Z, Y (vertical) 
#----------------------------------------------------------------------------------

# Baywidths
Baywidth_X = 5.0
Baywidth_Z = 6.0
Baywidth_Y = 4.0

# No. of Bays
No_of_Bays_X = 2
No_of_Bays_Z = 1
No_of_Bays_Y = 5

g = 9.81              
DL = 6.5              
LL = 2.5              
gamma_concrete = 24.0   

# Column properties
A = 0.2025
E = 3e7
G = 125e5
J = 0.007
Iy = 0.45**4 / 12
Iz = 0.45**4 / 12

# Beam X properties
A_x = 0.16
E_x = 3e7
G_x = 125e5
J_x = 0.004
Iy_x = 0.4**4 / 12
Iz_x = 0.4**4 / 12

# Beam Z properties
A_z = 0.18
E_z = 3e7
G_z = 125e5
J_z = 0.005
Iy_z = 0.45 * (0.4**3) / 12
Iz_z = 0.4 * (0.45**3) / 12

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

ops.geomTransf('Linear', Col_TransfTag, 0, 0, 1)  
ops.geomTransf('Linear', Beam_X_TransfTag, 0, 1, 0)  
ops.geomTransf('Linear', Beam_Z_TransfTag, 0, 1, 0)

#----------------------------------------------------------------------------------
# Elements
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

mass_floor = (DL + 0.25 * LL) / g 
mass_concrete = gamma_concrete / g 

for y in range(1, 6):               # y-grid lines 1, 2, 3, 4, 5
    if y == 5:                      # top 
        trib_y = 4.0 / 2.0
    else:
        trib_y = 4.0

    for z in [1, 2]:                # z-grid line 1, 2
        trib_z = 6.0 / 2.0

        for x in [1, 2, 3]:         # x-grid lines 1, 2, 3
            
            if x == 1 or x == 3:    # corner
                trib_x = 5.0 / 2.0  
            elif x == 2:
                trib_x = 5.0        
            
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
# Boundary Conditions 
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

weight_floor = DL + LL
weight_concrete = gamma_concrete

for y in range(1, 6):               # y-grid lines 1, 2, 3, 4, 5
    if y == 5:                      # top 
        trib_y = 4.0 / 2.0
    else:
        trib_y = 4.0

    for z in [1, 2]:                # z-grid line 1, 2
        trib_z = 6.0 / 2.0

        for x in [1, 2, 3]:         # x-grid lines 1, 2, 3
            
            if x == 1 or x == 3:    # corner
                trib_x = 5.0 / 2.0  
            elif x == 2:
                trib_x = 5.0        
            
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

output_folder = os.path.join("Tutorial", "Session_1_Results")
os.makedirs(output_folder, exist_ok=True)

ops.recorder('Node', '-file', os.path.join(output_folder, 'Acc312.out'), '-node', 312, '-dof', 1, 2, 3, 4, 5, 6, 'accel')
ops.recorder('Node', '-file', os.path.join(output_folder, 'Vel312.out'), '-node', 312, '-dof', 1, 2, 3, 4, 5, 6, 'vel')
ops.recorder('Node', '-file', os.path.join(output_folder, 'Disp312.out'), '-node', 312, '-dof', 1, 2, 3, 4, 5, 6, 'disp')
ops.recorder('Node', '-file', os.path.join(output_folder, 'R101.out'), '-node', 101, '-dof', 1, 2, 3, 4, 5, 6, 'reaction')
ops.recorder('Element', '-file', os.path.join(output_folder, 'F101111.out'), '-ele', 101111, 'force')

# --------------------------------------------------------------------------------
# Modal/Eigenvalue Analysis 
# --------------------------------------------------------------------------------

numModes = 6
lambdas = ops.eigen(numModes)  
ops.modalProperties('-file', f'{output_folder}/modalProperties.out')

periods = [(2 * np.pi) / (lam**0.5) for lam in lambdas]
print("Initial Time Periods:", [f"{p:.5f}" for p in periods])

def Plotter():
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

# --------------------------------------------------------------------------------
# Gravity Analysis 
# --------------------------------------------------------------------------------

ops.wipeAnalysis()

ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.algorithm('Linear')

NstepG = 1
ops.integrator('LoadControl', 1.0 / NstepG)
ops.analysis('Static')

ops.analyze(NstepG)
ops.loadConst('-time', 0.0)

Plotter()
ops.reactions()

R = 0.0
for node in [101,201,301,102,202,302]:
    R += ops.nodeReaction(node)[1]

print("Total Vertical Reaction =", R, "kN")

# --------------------------------------------------------------------------------
# Plotting Mode Shapes and Deformed Shape
# --------------------------------------------------------------------------------

def ModeShapesPlot():
    opsv.plot_defo()
    plt.title("Deformed Shape")

    for i in range(1, numModes+1):
        opsv.plot_mode_shape(i)
        plt.title(f"Mode {i}")

    plt.show()

# ModeShapesPlot()
