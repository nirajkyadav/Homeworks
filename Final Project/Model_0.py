"""
==============================================================================
OpenSees Finite Element Modeling Project

Date: June 30, 2026
Author: Niraj Yadav

Section Update: July 6, 2026

Units: Newtons (N), Millimeters (mm), Seconds (s)

Baseline model (Model 0): 
- This model uses displacement-based elements. 
- Plastic hinge elements length is one half of the depth of the associated section. 
- Other elements are of comparable length. 
- 5 Gauss-Lobotto integration points per element. 

Materials: Concrete02 & Steel02
==============================================================================
"""

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import opsvis as opsv
import os

ops.wipe()
ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)

#---------------------------------------------------------------------------------
# Materials, Units: N, mm, s
#---------------------------------------------------------------------------------

unconfined_concrete_tag = 1  # Cover concrete (unconfined)
confined_concrete_tag   = 2  # Core concrete (confined)
steel_beam_tag          = 3  # D16 Rebar
steel_col_tag           = 4  # HD16 Rebar

# uniaxialMaterial(  'Concrete02',                  matTag,   fpc,  epsc0, fpcu,   epsU, lambda,  ft,    Ets)
ops.uniaxialMaterial('Concrete02', unconfined_concrete_tag, -39.0, -0.002, -9.6, -0.004,    0.3, 2.0, 3600.0)
ops.uniaxialMaterial('Concrete02',   confined_concrete_tag, -48.0, -0.003, -9.6, -0.010,    0.3, 2.0, 3600.0)

# uniaxialMaterial(  'Steel02',         matTag,    Fy,       E0,     b,              *params)
ops.uniaxialMaterial('Steel02', steel_beam_tag, 295.0, 210000.0, 0.004, *[20.0, 0.925, 0.15])
ops.uniaxialMaterial('Steel02',  steel_col_tag, 500.0, 200000.0, 0.004, *[20.0, 0.925, 0.15])

#---------------------------------------------------------------------------------
# Section with quad element
#---------------------------------------------------------------------------------

def area(diameter):
    return (np.pi * diameter ** 2) / 4.0

# Beam Section (Section A-A)
beamWidth = 229
beamDepth = 457
beamCover = 42

ops.section('Fiber', 1)

ops.patch('rect', confined_concrete_tag,   6, 6, *[-186.5, -72.5], *[186.5,  72.5])  # core

# patch('quad', matTag, numSubdivIJ, numSubdivJK, *crdsI, *crdsJ, *crdsK, *crdsL)
ops.patch('quad', unconfined_concrete_tag, 2, 6, *[186.5, -72.5],   *[228.5, -114.5], *[228.5, 114.5], *[186.5, 72.5])      # right cover 
ops.patch('quad', unconfined_concrete_tag, 2, 6, *[-228.5, -114.5], *[-186.5, -72.5], *[-186.5, 72.5], *[-228.5, 114.5])    # left cover 
ops.patch('quad', unconfined_concrete_tag, 6, 2, *[-228.5, -114.5], *[228.5, -114.5], *[186.5, -72.5], *[-186.5, -72.5])    # bottom cover 
ops.patch('quad', unconfined_concrete_tag, 6, 2, *[-186.5, 72.5],   *[186.5, 72.5],   *[228.5, 114.5], *[-228.5, 114.5])    # top cover 

ops.layer('straight', steel_beam_tag, 3, area(16), *[186.5,  72.5],  *[186.5,  -72.5])
ops.layer('straight', steel_beam_tag, 2, area(16), *[153.5,  72.5],  *[153.5,  -72.5])
ops.layer('straight', steel_beam_tag, 2, area(16), *[-186.5, 72.5],  *[-186.5, -72.5])

# Column Section (Section B-B and C-C)
colWidth = 305
colDepth = 406
colCover = 43

ops.section('Fiber', 2)

ops.patch('rect',   confined_concrete_tag, 6, 6, *[-160.0, -109.5], *[ 160.0, 109.5])  # core

ops.patch('quad', unconfined_concrete_tag, 2, 6, *[160.0, -109.5],  *[203.0, -152.5],  *[203.0, 152.5],  *[160.0, 109.5])    # right cover
ops.patch('quad', unconfined_concrete_tag, 2, 6, *[-203.0, -152.5], *[-160.0, -109.5], *[-160.0, 109.5], *[-203.0, 152.5])   # left cover
ops.patch('quad', unconfined_concrete_tag, 6, 2, *[-203.0, -152.5], *[203.0, -152.5],  *[160.0, -109.5], *[-160.0, -109.5])  # bottom cover
ops.patch('quad', unconfined_concrete_tag, 6, 2, *[-160.0, 109.5],  *[160.0, 109.5],   *[203.0, 152.5],  *[-203.0, 152.5])   # top cover

ops.layer('straight', steel_col_tag, 3, area(16), *[160.0, 109.5],  *[160.0, -109.5])
ops.layer('straight', steel_col_tag, 2, area(16), *[  0.0, 109.5],  *[  0.0, -109.5])
ops.layer('straight', steel_col_tag, 3, area(16), *[-160.0,109.5],  *[-160.0,-109.5])

#---------------------------------------------------------------------------------
# Section without quad element
#---------------------------------------------------------------------------------

# def area(diameter):
#     return (np.pi * diameter ** 2) / 4.0

# # Beam Section (Section A-A)
# beamWidth = 229
# beamDepth = 457
# beamCover = 42

# ops.section('Fiber', 1)

# ops.patch('rect', confined_concrete_tag,   10, 1, *[-186.5, -114.5], *[186.5,  114.5])  
# ops.patch('rect', unconfined_concrete_tag, 1,  1, *[-228.5, -114.5], *[-186.5, 114.5]) 
# ops.patch('rect', unconfined_concrete_tag, 1,  1, *[186.5,  -114.5], *[228.5,  114.5]) 

# ops.layer('straight', steel_beam_tag, 3, area(16), *[186.5, 72.5],  *[186.5, -72.5])
# ops.layer('straight', steel_beam_tag, 2, area(16), *[153.5, 72.5],  *[153.5, -72.5])
# ops.layer('straight', steel_beam_tag, 2, area(16), *[-186.5, 72.5], *[-186.5, -72.5])

# # Column Section (Section B-B and C-C)
# colWidth = 305
# colDepth = 406
# colCover = 43

# ops.section('Fiber', 2)

# ops.patch('rect',   confined_concrete_tag, 10, 1, *[-160.0, -152.5], *[ 160.0, 152.5])
# ops.patch('rect', unconfined_concrete_tag,  1, 1, *[-203.0, -152.5], *[-160.0, 152.5])
# ops.patch('rect', unconfined_concrete_tag,  1, 1, *[ 160.0, -152.5], *[ 203.0, 152.5])

# ops.layer('straight', steel_col_tag, 3, area(16), *[160.0, 109.5],  *[160.0, -109.5])
# ops.layer('straight', steel_col_tag, 2, area(16), *[  0.0, 109.5],  *[  0.0, -109.5])
# ops.layer('straight', steel_col_tag, 3, area(16), *[-160.0,109.5],  *[-160.0,-109.5])

#---------------------------------------------------------------------------------
# Nodes (Origin at Joint Center), 1_ for right, 2_ for up, 3_ for left, 4_ for down
#---------------------------------------------------------------------------------

ops.node(0, 0.0, 0.0)                           # Origin

# Nodes on Right
ops.node(10,  203.0, 0.0)                       # Right Rigid Link Node
ops.node(11,  203.0 + 228.0, 0.0)               # Given Plastic Hinge Node
ops.node(12,  203.0 + 228.0 + (417/2)*1.0, 0.0)
ops.node(13,  203.0 + 228.0 + (417/2)*2.0, 0.0) # Right gravity load point P
ops.node(14,  848.0 + ((2119-848)/6)*1.0, 0.0)
ops.node(15,  848.0 + ((2119-848)/6)*2.0, 0.0)
ops.node(16,  848.0 + ((2119-848)/6)*3.0, 0.0)
ops.node(17,  848.0 + ((2119-848)/6)*4.0, 0.0)
ops.node(18,  848.0 + ((2119-848)/6)*5.0, 0.0)
ops.node(19,  848.0 + ((2119-848)/6)*6.0, 0.0)  # Right roller support

# Nodes on Left
ops.node(30,  -203.0, 0.0)                       # Left Rigid Link Node
ops.node(31,  -203.0 - 228.0, 0.0)               # Given Plastic Hinge Node
ops.node(32,  -203.0 - 228.0 - (417/2)*1.0, 0.0)
ops.node(33,  -203.0 - 228.0 - (417/2)*2.0, 0.0) # Left gravity load point P
ops.node(34,  -848.0 - ((2119-848)/6)*1.0, 0.0)
ops.node(35,  -848.0 - ((2119-848)/6)*2.0, 0.0)
ops.node(36,  -848.0 - ((2119-848)/6)*3.0, 0.0)
ops.node(37,  -848.0 - ((2119-848)/6)*4.0, 0.0)
ops.node(38,  -848.0 - ((2119-848)/6)*5.0, 0.0)
ops.node(39,  -848.0 - ((2119-848)/6)*6.0, 0.0)  # Left roller support

# Nodes on Top
ops.node(20, 0.0,  228.0)                        # Top joint face
ops.node(21, 0.0,  228.0 + 203.0)                # Given Plastic Hinge Node
ops.node(22, 0.0,  431.0 + ((1236-431)/4)*1.0)     
ops.node(23, 0.0,  431.0 + ((1236-431)/4)*2.0)     
ops.node(24, 0.0,  431.0 + ((1236-431)/4)*3.0)     
ops.node(25, 0.0,  431.0 + ((1236-431)/4)*4.0)   # Top tip loading node

# Nodes on Bottom
ops.node(40, 0.0,  -228.0)                        # Bottom joint face
ops.node(41, 0.0,  -228.0 - 203.0)                # Given Plastic Hinge Node
ops.node(42, 0.0,  -431.0 - ((1236-431)/4)*1.0)     
ops.node(43, 0.0,  -431.0 - ((1236-431)/4)*2.0)     
ops.node(44, 0.0,  -431.0 - ((1236-431)/4)*3.0)     
ops.node(45, 0.0,  -431.0 - ((1236-431)/4)*4.0)   # Bottom Hinge Support

# Boundary Conditions
ops.fix(19, 0, 1, 0)           # Right roller support 
ops.fix(39, 0, 1, 0)           # Left roller support 
ops.fix(45, 1, 1, 0)           # Bottom hinge support

#---------------------------------------------------------------------------------
# Elements
#---------------------------------------------------------------------------------

# Rigid Links
ops.rigidLink('beam', 0, 10)
ops.rigidLink('beam', 0, 20)
ops.rigidLink('beam', 0, 30)
ops.rigidLink('beam', 0, 40)

ops.geomTransf('Linear', 1)

# beamIntegration(. 'Lobatto', tag, secTag, N)
ops.beamIntegration('Lobatto',   1,      1, 5)  # For Beam (Section 1)
ops.beamIntegration('Lobatto',   2,      2, 5)  # For Column (Section 2)

# element('dispBeamColumn', eleTag, *eleNodes, transfTag, integrationTag)

# Right Beam Elements
ops.element('dispBeamColumn', 1011, *[10,11],          1,              1)
ops.element('dispBeamColumn', 1112, *[11,12],          1,              1)
ops.element('dispBeamColumn', 1213, *[12,13],          1,              1)
ops.element('dispBeamColumn', 1314, *[13,14],          1,              1)
ops.element('dispBeamColumn', 1415, *[14,15],          1,              1)
ops.element('dispBeamColumn', 1516, *[15,16],          1,              1)
ops.element('dispBeamColumn', 1617, *[16,17],          1,              1)
ops.element('dispBeamColumn', 1718, *[17,18],          1,              1)
ops.element('dispBeamColumn', 1819, *[18,19],          1,              1)

# Left Beam Elements
ops.element('dispBeamColumn', 3031, *[31,30],          1,              1)
ops.element('dispBeamColumn', 3132, *[32,31],          1,              1)
ops.element('dispBeamColumn', 3233, *[33,32],          1,              1)
ops.element('dispBeamColumn', 3334, *[34,33],          1,              1)
ops.element('dispBeamColumn', 3435, *[35,34],          1,              1)
ops.element('dispBeamColumn', 3536, *[36,35],          1,              1)
ops.element('dispBeamColumn', 3637, *[37,36],          1,              1)
ops.element('dispBeamColumn', 3738, *[38,37],          1,              1)
ops.element('dispBeamColumn', 3839, *[39,38],          1,              1)

# Top Column Elements
ops.element('dispBeamColumn', 2021, *[20,21],          1,              2) 
ops.element('dispBeamColumn', 2122, *[21,22],          1,              2)
ops.element('dispBeamColumn', 2223, *[22,23],          1,              2)
ops.element('dispBeamColumn', 2324, *[23,24],          1,              2)
ops.element('dispBeamColumn', 2425, *[24,25],          1,              2)

# Bottom Column Elements
ops.element('dispBeamColumn', 4041, *[41,40],          1,              2)
ops.element('dispBeamColumn', 4142, *[42,41],          1,              2)
ops.element('dispBeamColumn', 4243, *[43,42],          1,              2)
ops.element('dispBeamColumn', 4344, *[44,43],          1,              2)
ops.element('dispBeamColumn', 4445, *[45,44],          1,              2)

#---------------------------------------------------------------------------------
# Recorders
#---------------------------------------------------------------------------------

output_folder = os.path.join("Final_Project", "Results", "Model_0")
os.makedirs(output_folder, exist_ok=True)

# Node displacement and reaction force tracking
ops.recorder('Node', '-file', os.path.join(output_folder, 'col_tip_disp.out'), '-time', '-node', 25, '-dof', 1, 'disp')
ops.recorder('Node', '-file', os.path.join(output_folder, 'base_reaction.out'),'-time', '-node', 45, '-dof', 1, 'reaction')

# Cross-section force and deformation for two left beam elements as required
ops.recorder('Element', '-file', os.path.join(output_folder, 'ele3031_sec1_force.out'), '-time', '-ele', 3031, 'section', 1, 'force')
ops.recorder('Element', '-file', os.path.join(output_folder, 'ele3031_sec5_force.out'), '-time', '-ele', 3031, 'section', 5, 'force')
ops.recorder('Element', '-file', os.path.join(output_folder, 'ele3031_sec1_def.out'),   '-time', '-ele', 3031, 'section', 1, 'deformation')
ops.recorder('Element', '-file', os.path.join(output_folder, 'ele3031_sec5_def.out'),   '-time', '-ele', 3031, 'section', 5, 'deformation')

ops.recorder('Element', '-file', os.path.join(output_folder, 'ele3132_sec1_force.out'), '-time', '-ele', 3132, 'section', 1, 'force')
ops.recorder('Element', '-file', os.path.join(output_folder, 'ele3132_sec5_force.out'), '-time', '-ele', 3132, 'section', 5, 'force')
ops.recorder('Element', '-file', os.path.join(output_folder, 'ele3132_sec1_def.out'),   '-time', '-ele', 3132, 'section', 1, 'deformation')
ops.recorder('Element', '-file', os.path.join(output_folder, 'ele3132_sec5_def.out'),   '-time', '-ele', 3132, 'section', 5, 'deformation')

# Steel rebar stress and strain for top layer of two left beam elements as required
ops.recorder('Element', '-file', os.path.join(output_folder, 'ele3031_sec1_topRebar.out'), '-time', '-ele', 3031, 'section', 1, 'fiber', 186.5, 72.5, steel_beam_tag, 'stressStrain')
ops.recorder('Element', '-file', os.path.join(output_folder, 'ele3031_sec5_topRebar.out'), '-time', '-ele', 3031, 'section', 5, 'fiber', 186.5, 72.5, steel_beam_tag, 'stressStrain')
ops.recorder('Element', '-file', os.path.join(output_folder, 'ele3132_sec1_topRebar.out'), '-time', '-ele', 3132, 'section', 1, 'fiber', 186.5, 72.5, steel_beam_tag, 'stressStrain')
ops.recorder('Element', '-file', os.path.join(output_folder, 'ele3132_sec5_topRebar.out'), '-time', '-ele', 3132, 'section', 5, 'fiber', 186.5, 72.5, steel_beam_tag, 'stressStrain')

#---------------------------------------------------------------------------------
# Gravity Loads
#---------------------------------------------------------------------------------

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

# Downward symmetric vertical loads (55 kN = 55000 N)
ops.load(13, 0.0, -55000.0, 0.0)
ops.load(33, 0.0, -55000.0, 0.0)

# Analysis settings for gravity load application
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('RelativeNormDispIncr', 1e-4, 50)
ops.algorithm('Newton')
ops.integrator('LoadControl', 0.1)
ops.analysis('Static')
ops.analyze(10)

ops.loadConst('-time', 0.0)           # Set the time to zero an hold the loads constant

#---------------------------------------------------------------------------------
# Cyclic Lateral Loading Phase (Controlled at Node 25, DOF 1)
#---------------------------------------------------------------------------------

# Array containing the sequence of targeted drift peak displacement values (mm)
peakDisp = [
    11.1, -11.1, 
    29.6, -29.6, 29.6, -29.6, 
    44.4, -44.4, 44.4, -44.4, 
    59.2, -59.2, 59.2, -59.2, 
    74.0, -74.0, 74.0, -74.0, 
    88.8, -88.8, 88.8, -88.8, 
    103.6, -103.6, 103.6, -103.6
]

respnode = 25    # Tip node used to monitor and execute displacement control
lateral_dof = 1  # Horizontal translation degree of freedom (UX)
dK0 = 0.01        # Initial Displacement increment step size (mm)

ops.timeSeries('Linear', 2)
ops.pattern('Plain', 2, 2)
ops.load(respnode, 1.0, 0.0, 0.0)  # Reference lateral unit load

#---------------------------------------------------------------------------------
# Analysis
#---------------------------------------------------------------------------------

currentDisp = 0.0

# Step sizes
step_factors = [1.0, 0.5, 0.1, 0.05, 0.01, 0.001]

convergence_tests = [
    # ('RelativeNormDispIncr', 1.0e-4),
    ('RelativeNormUnbalance', 1.0e-2),
    ('RelativeEnergyIncr', 1.0e-3) 
]

solution_algorithms = [
    # ('NewtonLineSearch', 0.8), 
    ('KrylovNewton', None),    
    ('Broyden', 8),                    
    ('BFGS', None),   
    ('ModifiedNewton', None),     
    ('ModifiedNewton', '-initial'),   
    ('Newton', None)                  # Standard Newton-Raphson
]

for i in range(len(peakDisp)):

    cycleDisp = peakDisp[i] - currentDisp
    sign = 1 if cycleDisp > 0 else -1
    dK = dK0 * sign

    numiter = round(abs(cycleDisp / dK))
    print(f"\n---> Running {numiter} steps with dK={dK:.4f} mm to advance from {currentDisp:.4f} mm to {peakDisp[i]:.2f} mm")

    # Initial analysis configurations
    ops.integrator('DisplacementControl', respnode, lateral_dof, dK)

    # ops.test('RelativeTotalNormDispIncr', 1.0e-4, 100)
    # ops.test('RelativeNormUnbalance', 1.0e-2, 100)
    ops.test('RelativeNormDispIncr', 1.0e-4, 100)

    # ops.algorithm('Newton')
    ops.algorithm('NewtonLineSearch', 0.8) 
    
    ops.analysis('Static')
    
    ok1 = ops.analyze(numiter)
    currentDisp = ops.nodeDisp(respnode, lateral_dof)

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
                ops.integrator('DisplacementControl', respnode, lateral_dof, new_dK)
                
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
                            print(f"   [SUCCESS] Step advanced via: {test_name} | {algo_name} | Step Factor: {factor}")
                            step_converged = True
                            break  # Exit algorithm loop
            
            # If still no convergence reached
            if not step_converged:
                print(f"\nCRITICAL ERROR: NO convergence at displacement {currentDisp:.6f}")
                print("Analysis terminated.")
                break
            
            currentDisp = ops.nodeDisp(respnode, lateral_dof)
            
        print(f"Reached peak target: {currentDisp:.4f} mm")
    else:
        print(f"Reached peak target: {currentDisp:.4f} mm")

