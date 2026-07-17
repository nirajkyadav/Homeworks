"""
================================================================================
Date: 15 July 2026
Author: Niraj Yadav

Units: N, mm, s

Cyclic loading: Moment-curvature diagram of an unsymmetrical section

Beam Section from Final Project
================================================================================
"""

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import opsvis as opsv
import os

ops.wipe()
ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)  # 2D model with 3 DOFs per node (Ux, Uy, Rz)

# ---------------------------------------------------------------------------------
# Nodes
# Two coincident nodes at the origin; zeroLengthSection element connects them
# ---------------------------------------------------------------------------------

ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, 0.0)

ops.fix(1, 1, 1, 1)        # fixing the left node in all dofs
ops.fix(2, 0, 1, 0)        # fixing the right in translation along y

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

#---------------------------------------------------------------------------------
# Element
#---------------------------------------------------------------------------------

# element(  'zeroLengthSection', eleTag, *eleNodes, secTag)
ops.element('zeroLengthSection',      1,   *[1, 2], SecTag)

#---------------------------------------------------------------------------------
# Loads
#---------------------------------------------------------------------------------

ops.timeSeries('Linear', 1)                   # Linearly increasing time series (tag = 1)
ops.pattern('Plain', 2, 1, '-fact', 1.0)      # Load pattern (tag = 2) referencing time series 1
ops.load(2, 0.0, 0.0, 1.0)                    # Unit moment applied at node 2, DOF 3 (rotation)

#---------------------------------------------------------------------------------
# Recorder
#---------------------------------------------------------------------------------

output_folder = os.path.join("Final_Project", "steps_back", "Results")  # Output directory for results
os.makedirs(output_folder, exist_ok=True)                               # Create folder if it does not exist

# Record rotation (DOF 3) at node 2 at every step
ops.recorder('Node', '-file', os.path.join(output_folder, "Cyclic_mom_curv.out"), '-time', '-node', 2, '-dof', 3, 'disp')

#---------------------------------------------------------------------------------
# Analysis
#---------------------------------------------------------------------------------
ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGeneral')

phiy         = 4.259419e-06  # Yield curvature from Monotonic Moment Curvature Analysis
numinc       = 100           # Number of displacement-controlled increments

ops.algorithm('Newton')
ops.test('RelativeNormDispIncr', 1e-4, 100)
# ops.test('RelativeNormUnbalance', 1e-2, 100)
# ops.test('RelativeEnergyIncr', 1e-3, 100)

peaks = [0, 2, -2, 4, -4, 6, -6]  # Curvature amplitude targets [multiples of phiy]

curvature = []  # List to store curvature values (rotation at node 2, DOF 3) 
moment    = []  # List to store moment values (reaction at node 1, DOF 3) 

for i in range(len(peaks) - 1):
    """
    Step between each peaks for cyclic analysis
    cdK0 is positive during loading and negative during reversal.
    """
    start = peaks[i]      # Curvature amplitude multiple at start of this half-cycle
    end   = peaks[i + 1]  # Curvature amplitude multiple at end of this half-cycle

    # Curvature increment per step for this half-cycle
    cdK0 = (end - start) * phiy / numinc

    # Update the integrator with the new curvature increment for this half-cycle
    ops.integrator('DisplacementControl', 2, 3, cdK0)
    ops.analysis('Static')

    for j in range(numinc):
        """ Step through each curvature increment and collect M-φ data """
        ok = ops.analyze(1)  # Perform one analysis step (0 = success, <0 = failure)

        if ok != 0:
            print(f"Analysis failed at half-cycle {i + 1}, increment {j + 1}. Trying backup algorithms...")
            
            print("Trying NewtonLineSearch...")
            ops.algorithm('NewtonLineSearch')
            ok = ops.analyze(1)
            
            if ok != 0:
                print("NewtonLineSearch failed. Trying Broyden...")
                ops.algorithm('Broyden')
                ok = ops.analyze(1)
                
            if ok != 0:
                print("Broyden failed. Trying ModifiedNewton...")
                ops.algorithm('ModifiedNewton')
                ok = ops.analyze(1)
                
            if ok != 0:
                print("Broyden failed. Trying BFGS...")
                ops.algorithm('BFGS')
                ok = ops.analyze(1)
                
            if ok != 0:
                print("Broyden failed. Trying KrylovNewton...")
                ops.algorithm('KrylovNewton')
                ok = ops.analyze(1)
            
            ops.algorithm('Newton')
            
            if ok == 0:
                print("Backup algorithm successful! Resuming standard analysis.")
            else:
                print("All backup algorithms failed. Analysis terminating.")
                break # Exit the inner loop if it still fails

        if ok != 0:
            break # Exit the outer loop if analysis completely died

        phi = ops.nodeDisp(2, 3)    # Curvature: rotation (DOF 3) at free node 2 

        ops.reactions()             # Compute reactions at fixed nodes
        M = ops.nodeReaction(1, 3)  # Moment reaction at fixed node 1, DOF 3 

        curvature.append(phi * 1e3)
        moment.append(-M * 1e-6)      # Negative because the support provides resisting moment

# ---------------------------------------------------------------------------------
# Moment Curvature Diagram
# ---------------------------------------------------------------------------------

plt.figure(figsize=(8, 6))
plt.plot(curvature, moment)
plt.axhline(0, color='black', linewidth=0.8)
plt.axvline(0, color='black', linewidth=0.8)
plt.xlabel('Curvature (1/m)')
plt.ylabel('Moment (kN·m)')
plt.title('Moment–Curvature Diagram (Unsymmetrical Section)')
plt.grid(True, alpha=0.6)
plt.tight_layout()
plt.show()
