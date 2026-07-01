"""
================================================================================
OpenSees Tutorial - Session 3: Quasi-Static cyclic analysis of RC cantilever Beam
Date: 26 June 2026
Author: Niraj Yadav

Reference: Tcl code developed by Abdelrahamn using dispBeamColumn for Popov et al 1972

Source: https://www.youtube.com/watch?v=b64Izd2uq54&list=PL7ZAB0cr0KbOJqUeAw3sjjtmm_5FWpmi2

Units: in, kip, ksi
================================================================================
"""

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import opsvis as opsv
import os

ops.wipe()
ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)  # 2D model, 3 DOFs per node (Ux, Uy, Rz)

# ---------------------------------------------------------------------------------
# Material Definitions
# ---------------------------------------------------------------------------------
unconfined_concrete_tag = 1
confined_concrete_tag   = 2
steel_tag               = 3

# Concrete parameters (Units: kip, in, ksi)
fpc    = -5.03           # Unconfined concrete compressive strength [ksi]
fpcu   = 0.2 * fpc       # Concrete crushing strength [ksi]
epsc0  = -0.002          # Strain at peak unconfined strength
epsU   = -0.004          # Ultimate strain for unconfined concrete
fpcc   = -5.73           # Confined concrete compressive strength [ksi]
fpccu  = 0.2 * fpcc      # Confined concrete crushing strength [ksi]
epscc0 = -0.003          # Strain at peak confined strength
epscU  = -0.021          # Ultimate strain for confined concrete
ft     = 0.28            # Concrete tensile strength [ksi]
Ec     = 4300.0          # Initial elastic modulus of concrete [ksi]
Ets    = 0.1 * Ec        # Tension softening stiffness [ksi]

# Steel parameters
Fy     = 60.0            # Yield stress of reinforcing steel [ksi]
Es     = 29000.0         # Elastic modulus of steel [ksi]

# Define uniaxial material models
ops.uniaxialMaterial('Concrete02', unconfined_concrete_tag,  fpc,  epsc0,  fpcu,  epsU,   0.1,  ft,  Ets)     
ops.uniaxialMaterial('Concrete02', confined_concrete_tag,    fpcc, epscc0, fpccu, epscU,  0.1,  ft,  Ets)     
ops.uniaxialMaterial('Steel02', steel_tag,  Fy,  Es,  0.1, 20, 0.925, 0.15, 0.005, 1.0, 0.005, 1.0)           

# ---------------------------------------------------------------------------------
# Section & Fiber Geometry
# ---------------------------------------------------------------------------------
colWidth = 15.0          # Total column section width [in]
colDepth = 29.0          # Total column section depth [in]
dia      = 1.128         # Nominal diameter of longitudinal bars [in]
As       = 1.0           # Area of an individual steel bar [in^2]
c        = 2.5           # Clear concrete cover [in]
sv       = 3.0           # Vertical center-to-center spacing between bar layers [in]

# Coordinate tracking boundaries relative to section centroid
y1 = colDepth / 2.0
y2 = y1 - c - dia / 2.0  # Location of outer reinforcing layer
y3 = y2 - sv             # Location of inner reinforcing layer
z1 = colWidth / 2.0
z2 = z1 - c - dia / 2.0  # Lateral position of bars

ops.section('Fiber', 1)
# Concrete patches: patch('rect', matTag, numSubdivY, numSubdivZ, bottom_y, left_z, top_y, right_z)
ops.patch('rect', unconfined_concrete_tag,  2,  1,   -y1, -z1,  -y2,  z1)       # Bottom cover
ops.patch('rect', unconfined_concrete_tag,  2,  1,    y2, -z1,   y1,  z1)       # Top cover
ops.patch('rect', confined_concrete_tag, 15,  1,   -y2, -z1,   y2,  z1)         # Confined core

# Reinforcing steel layers: layer('straight', matTag, numBars, barArea, start_y, start_z, end_y, end_z)
ops.layer('straight', steel_tag,  4,  As,   -y2,  -z2,  -y2,   z2)              # Bottom layer
ops.layer('straight', steel_tag,  2,  As,   -y3,  -z2,  -y3,   z2)              # Lower-middle layer
ops.layer('straight', steel_tag,  2,  As,    y3,  -z2,   y3,   z2)              # Upper-middle layer
ops.layer('straight', steel_tag,  4,  As,    y2,  -z2,   y2,   z2)              # Top layer

# ---------------------------------------------------------------------------------
# Nodes, Boundary Conditions, and Elements
# ---------------------------------------------------------------------------------
L  = 70.0 / 3.0          # Length segment of column element [in]
Lr = 8.0                 # Rigid arm length extension [in]

ops.node(1,   0.0,        0.0)
ops.node(2,   L,          0.0)
ops.node(3,   2 * L,      0.0)
ops.node(4,   3 * L,      0.0)
ops.node(5,   3 * L + Lr, 0.0)  # Node where displacement loading is applied

ops.fix(1,  1, 1, 1)     # Base cantilever node fully fixed
ops.geomTransf('Linear', 1)

no_of_integration_points = 7
ops.beamIntegration('Lobatto', 1, 1, no_of_integration_points)

# element('dispBeamColumn', eleTag, *eleNodes, transfTag, integrationTag)
ops.element('dispBeamColumn', 1, 1, 2, 1, 1)
ops.element('dispBeamColumn', 2, 2, 3, 1, 1)
ops.element('dispBeamColumn', 3, 3, 4, 1, 1)

# rigidLink(type, rNodeTag, cNodeTag)
ops.rigidLink('beam', 5, 4) 

"""
# ops.rigidLink('beam', 4, 5)
This makes, Node 4 a master node and Node 5 the slave node. With the Transformation constraint handler, the DOFs for the slave node are eliminated from the global system matrix. So, DisplacementControl applied on respnode (Node 5), is internally transferred to the master node (Node 4). Thus, Node 4 moves exactly with step size (dK), but because of the rigid body rotation (θ) of the arm, Node 5 undergoes additional displacement (u_5 = u_4 + L_r.θ), yielding greater current displacement than target displacement. 
Solution: we make Node 5 master, and Node 4 the slave.
"""

# ---------------------------------------------------------------------------------
# Recorders & Output Directories
# ---------------------------------------------------------------------------------

output_folder = os.path.join("Tutorial", "Session_3_Results", "disp_Beam_Column")
os.makedirs(output_folder, exist_ok=True)

ops.recorder('Node', '-file', os.path.join(output_folder, 'node5disp.out'),     '-time', '-node', 5, '-dof', 2, 'disp')
ops.recorder('Node', '-file', os.path.join(output_folder, 'node1reaction.out'), '-time', '-node', 1, '-dof', 2, 'reaction')

# ---------------------------------------------------------------------------------
# Cyclic Loading
# ---------------------------------------------------------------------------------

peakDisp = [
    1.30E-02, -1.30E-02, 1.30E-02, -1.30E-02,
    2.31E-02, -2.31E-02, 2.31E-02, -2.31E-02,
    7.75E-02, -7.75E-02, 7.75E-02, -7.75E-02,
    1.22E-01, -1.22E-01, 1.22E-01, -1.22E-01,
    0.6, -0.472329472, 0.6, -0.424317174,
    1.0823332332, -1.117653368, 0.888066638, -1.065672816,
    1.495960246, -1.490955241, 1.495960246, -1.490955241,
    1.81, -1.77034177, 1.81, -1.77034177, 
    2.4, 
    -2.391856142, 2.38967539, -2.391856142,
    2.68,
    -2.869297869,
    3.75,
    -3.445659946,
    2.99,
    -3.503217503,
]

respnode = 5

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

ops.load(5, 0.0, 1.0, 0.0)  # Reference lateral load applied at Node 5 (DOF 2)

dK0 = 0.001                 # increment

# ---------------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------------

currentDisp = 0.0
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('BandGeneral')

# Step sizes
step_factors = [1.0, 0.5, 0.1, 0.05, 0.01, 0.001]

convergence_tests = [
    ('RelativeTotalNormDispIncr', 1.0e-4),
    # ('RelativeNormDispIncr', 1.0e-4),
    ('RelativeNormUnbalance', 1.0e-2),
    ('RelativeEnergyIncr', 1.0e-3) 
]

solution_algorithms = [
    ('Newton', None),                  # Standard Newton-Raphson
    ('NewtonLineSearch', 0.8),     
    ('Broyden', 8),                    
    ('BFGS', None),   
    ('ModifiedNewton', None),     
    ('ModifiedNewton', '-initial')   
]

for i in range(len(peakDisp)):
    cycleDisp = peakDisp[i] - currentDisp
    sign = 1 if cycleDisp > 0 else -1
    dK = dK0 * sign

    numiter = round(abs(cycleDisp / dK))
    print(f"Running {numiter} steps with dK={dK:.6f} to go from {currentDisp:.6f} to {peakDisp[i]:.6f}")

    # Initial analysis configurations
    ops.integrator('DisplacementControl', respnode, 2, dK)
    ops.test('RelativeNormDispIncr', 1.0e-4, 100)
    ops.algorithm('Newton')
    ops.analysis('Static')
    
    ok1 = ops.analyze(numiter)
    currentDisp = ops.nodeDisp(respnode, 2)

    # If the analysis fails
    if ok1 != 0:
        print(f"Block execution failed. Engaging fallback matrix...")
        
        while (peakDisp[i] - currentDisp) * sign > 0:
            step_converged = False
            
            # 1. Loop through step subdivisions (Full size down to 1/1000th size)
            for factor in step_factors:
                if step_converged:
                    break
                
                new_dK = dK * factor
                ops.integrator('DisplacementControl', respnode, 2, new_dK)
                
                # 2. Loop through different relative test criteria
                for test_name, tolerance in convergence_tests:
                    if step_converged:
                        break
                        
                    ops.test(test_name, tolerance, 100, 0)
                    
                    # 3. Loop through different solution algorithms
                    for algo_name, param in solution_algorithms:
                        if param is not None:
                            ops.algorithm(algo_name, param)
                        else:
                            ops.algorithm(algo_name)
                        
                        ok = ops.analyze(1)
                        
                        if ok == 0:
                            print(f"   [SUCCESS] Step advanced via: {test_name} | {algo_name} | Step Size Factor: {factor}")
                            step_converged = True
                            break # Exit algorithm loop
            
            # If still no convergence reached
            if not step_converged:
                print(f"\nCRITICAL ERROR: NO convergence at displacement {currentDisp:.6f}")
                print("Analysis terminated.")
                break
            
            currentDisp = ops.nodeDisp(respnode, 2)
            
        print(f"Reached target: {currentDisp:.6f}")
    else:
        print(f"Reached target: {currentDisp:.6f}")