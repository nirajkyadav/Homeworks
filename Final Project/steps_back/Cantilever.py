"""
================================================================================
Date: 15 July 2026 
Author: Niraj Yadav

Units: N, mm, s

Cantilever with unsymmetrical rectangular section
================================================================================
"""

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import opsvis as opsv
import os

ops.wipe()
ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)  # 2D model, 3 DOFs per node (Ux, Uy, Rz)

#---------------------------------------------------------------------------------
# Materials
#---------------------------------------------------------------------------------

unconfined_concrete_tag = 1  # Cover concrete (unconfined)
confined_concrete_tag   = 2  # Core concrete (confined)
steel_beam_tag          = 3  # D16 Rebar

# uniaxialMaterial(  'Concrete02',                  matTag,   fpc,  epsc0, fpcu,   epsU, lambda,  ft,    Ets)
ops.uniaxialMaterial('Concrete02', unconfined_concrete_tag, -39.0, -0.002, -9.6, -0.004,    0.3, 2.0, 3600.0)
ops.uniaxialMaterial('Concrete02',   confined_concrete_tag, -48.0, -0.003, -9.6, -0.010,    0.3, 2.0, 3600.0)

# uniaxialMaterial(  'Steel02',         matTag,    Fy,       E0,     b,              *params)
ops.uniaxialMaterial('Steel02', steel_beam_tag, 295.0, 210000.0, 0.004, *[20.0, 0.925, 0.15])        

#---------------------------------------------------------------------------------
# Section
#---------------------------------------------------------------------------------

def area(diameter):
    return (np.pi * diameter ** 2) / 4.0

# Beam Section (Section A-A)
beamWidth = 229
beamDepth = 457
beamCover = 42

SecTag = 1

ops.section('Fiber', SecTag)

ops.patch('rect', confined_concrete_tag,   6, 6, *[-186.5, -72.5], *[186.5,  72.5])  # core

# patch('quad', matTag, numSubdivIJ, numSubdivJK, *crdsI, *crdsJ, *crdsK, *crdsL)
ops.patch('quad', unconfined_concrete_tag, 2, 6, *[186.5, -72.5],   *[228.5, -114.5], *[228.5, 114.5], *[186.5, 72.5])      # right cover 
ops.patch('quad', unconfined_concrete_tag, 2, 6, *[-228.5, -114.5], *[-186.5, -72.5], *[-186.5, 72.5], *[-228.5, 114.5])    # left cover 
ops.patch('quad', unconfined_concrete_tag, 6, 2, *[-228.5, -114.5], *[228.5, -114.5], *[186.5, -72.5], *[-186.5, -72.5])    # bottom cover 
ops.patch('quad', unconfined_concrete_tag, 6, 2, *[-186.5, 72.5],   *[186.5, 72.5],   *[228.5, 114.5], *[-228.5, 114.5])    # top cover 

ops.layer('straight', steel_beam_tag, 3, area(16), *[186.5,  72.5],  *[186.5,  -72.5])
ops.layer('straight', steel_beam_tag, 2, area(16), *[153.5,  72.5],  *[153.5,  -72.5])
ops.layer('straight', steel_beam_tag, 2, area(16), *[-186.5, 72.5],  *[-186.5, -72.5])

# ---------------------------------------------------------------------------------
# Nodes, Boundary Conditions, and Elements
# ---------------------------------------------------------------------------------

L_flex = 1916.0  # Flexible beam length (mm)
L_rigid = 203.0  # Rigid joint link (mm)

ops.node(1, 0.0, 0.0)
ops.node(4, L_flex, 0.0)
ops.node(5, L_flex + L_rigid, 0.0)


ops.fix(1,  1, 1, 1)     # Base cantilever node fully fixed
ops.geomTransf('Linear', 1)

no_of_integration_points = 5
ops.beamIntegration('Lobatto', 1, 1, no_of_integration_points)

ops.element('forceBeamColumn', 1, 1, 4, 1, 1, '-maxIter', 100, '-tol', 1e-7)
ops.rigidLink('beam', 5, 4)

# ---------------------------------------------------------------------------------
# Recorders & Output Directories
# ---------------------------------------------------------------------------------

output_folder = os.path.join("Final_Project", "steps_back", "Results")  # Output directory for results
os.makedirs(output_folder, exist_ok=True)

ops.recorder('Node', '-file', os.path.join(output_folder, 'node5disp.out'),     '-time', '-node', 5, '-dof', 2, 'disp')
ops.recorder('Node', '-file', os.path.join(output_folder, 'node1reaction.out'), '-time', '-node', 1, '-dof', 2, 'reaction')

# ---------------------------------------------------------------------------------
# Cyclic Loading
# ---------------------------------------------------------------------------------

peakDisp = [
    11.1, -11.1, 
    29.6, -29.6, 29.6, -29.6, 
    44.4, -44.4, 44.4, -44.4, 
    59.2, -59.2, 59.2, -59.2, 
    74.0, -74.0, 74.0, -74.0, 
    88.8, -88.8, 88.8, -88.8, 
    103.6, -103.6, 103.6, -103.6
]

respnode = 5

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

ops.load(5, 0.0, 1.0, 0.0)  # Reference lateral load applied at Node 5 (DOF 2)

dK0  = 0.1                  # Displacement increment

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
    ('KrylovNewton', None), 
    ('NewtonLineSearch', 0.8),     
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
    print(f"Running {numiter} steps with dK={dK:.6f} to go from {currentDisp:.6f} to {peakDisp[i]:.6f}")

    # Initial analysis configurations
    ops.integrator('DisplacementControl', respnode, 2, dK)
    ops.test('RelativeNormDispIncr', 1.0e-4, 100)
    ops.algorithm('KrylovNewton')
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
