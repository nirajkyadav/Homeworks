"""
================================================================================
OpenSees Tutorial - Session 2: Moment Curvature Analysis
Date: 25 June 2026
Author: Niraj Yadav

Units: in, kip, ksi

Assignment:
1) Cyclic loading
   Develop moment-curvature diagram
   • Using different steel and concrete materials
   • Using different nonlinear analysis algorithms
   • Using different convergence tests and tolerances
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
# Material
# Steel01    : bilinear steel model
# Steel02    : steel model with Bauschinger effect
# Concrete01 : Kent-Scott-Park model (no tensile strength)
# Concrete02 : Kent-Scott-Park model (with linear tension softening)
#---------------------------------------------------------------------------------

# Material tags
steel_01    = 1   
steel_02    = 2   
concrete_01 = 3  
concrete_02 = 4   

# uniaxialMaterial(  'Steel01',      tag,     fy,       E,    b)                      b = strain-hardening ratio
ops.uniaxialMaterial('Steel01', steel_01,    60., 29000.0, 0.01)
ops.uniaxialMaterial('Steel02', steel_02,    60., 29000.0, 0.01, *[18, 0.925, 0.15])  # *[R0, cR1, cR2]: Bauschinger curve transition parameters

# uniaxialMaterial(  'Concrete01',         tag,  fpc,  epsc0, fpcu,   epsU)
ops.uniaxialMaterial('Concrete01', concrete_01, -6.0, -0.004, -5.0, -0.014)

# uniaxialMaterial(  'Concrete02',         tag,  fpc,  epsc0, fpcu,   epsU,  ft,  Ets, lambda)
ops.uniaxialMaterial('Concrete02', concrete_02, -6.0, -0.004, -5.0, -0.014, 0.1, 0.56,    280)

#---------------------------------------------------------------------------------
# Section
#---------------------------------------------------------------------------------

# Geometry
col_width = 16.0             # column width
col_depth = 16.0             # column depth
cover = 2.5                  # cover
As = 1.0                     # area of steel
y1 = col_depth/2.0           # distance from neutral axis to top/bottom face
z1 = col_width/2.0           # distance from neutral axis to left/right face

SecTag = 1

ops.section('Fiber', SecTag)

# patch(  'rect',      matTag, numSubdivY, numSubdivZ,  -y,  -z, +y, +z)
ops.patch('rect', concrete_02,          9,          1, -y1, -z1, y1, z1)

# layer(  'straight',   matTag, numFibers, As,     y_start,     z_start,       y_end,      z_end)
ops.layer('straight', steel_02,         4, As,  y1 - cover,  z1 - cover,  y1 - cover, cover - z1)
ops.layer('straight', steel_02,         4, As,  cover - y1,  z1 - cover,  cover - y1, cover - z1)

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

output_folder = os.path.join("Tutorial", "Session_2_Results")  # Output directory for results
os.makedirs(output_folder, exist_ok=True)                      # Create folder if it does not exist

# Record rotation (DOF 3) at node 2 at every step
ops.recorder('Node', '-file', os.path.join(output_folder, "Cyclic.out"), '-time', '-node', 2, '-dof', 3, 'disp')

#---------------------------------------------------------------------------------
# Analysis
#---------------------------------------------------------------------------------
ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGeneral')

phiy  = 0.0021 / (0.6 * (col_depth - cover))  # Yield curvature (1/in): εy / (0.6 * effective depth)
numinc = 100                                  # Number of displacement-controlled increments

ops.algorithm('Newton')
ops.test('RelativeNormDispIncr', 1e-4, 100)
# ops.test('RelativeNormUnbalance', 1e-2, 100)
# ops.test('RelativeEnergyIncr', 1e-3, 100)

peaks = [0, 2, -2, 4, -4, 6, -6]  # Curvature amplitude targets [multiples of phiy]

curvature = []  # List to store curvature values (rotation at node 2, DOF 3) [1/in]
moment    = []  # List to store moment values (reaction at node 1, DOF 3) [kip·in]

for i in range(len(peaks) - 1):
    """
    Step between each peaks for cyclic analysis
    cdK0 is positive during loading and negative during reversal.
    """
    start = peaks[i]      # Curvature amplitude multiple at start of this half-cycle
    end   = peaks[i + 1]  # Curvature amplitude multiple at end of this half-cycle

    # Displacement increment per step for this half-cycle [1/in]
    cdK0 = (end - start) * phiy / numinc

    # Update the integrator with the new displacement increment for this half-cycle
    ops.integrator('DisplacementControl', 2, 3, cdK0)
    ops.analysis('Static')

    for j in range(numinc):
        """ Step through each displacement increment and collect M-φ data """
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
            
            ops.algorithm('Newton')
            
            if ok == 0:
                print("Backup algorithm successful! Resuming standard analysis.")
            else:
                print("All backup algorithms failed. Analysis terminating.")
                break # Exit the inner loop if it still fails

        if ok != 0:
            break # Exit the outer loop if analysis completely died

        phi = ops.nodeDisp(2, 3)    # Curvature: rotation (DOF 3) at free node 2 [1/in]

        ops.reactions()             # Compute reactions at fixed nodes
        M = ops.nodeReaction(1, 3)  # Moment reaction at fixed node 1, DOF 3 [kip·in]

        curvature.append(phi)
        moment.append(-M)           # Negative because the support provides resisting moment

# ---------------------------------------------------------------------------------
# Moment Curvature Diagram
# ---------------------------------------------------------------------------------

plt.figure(figsize=(8, 6))
plt.plot(curvature, moment)
plt.axhline(0, color='black', linewidth=0.8)
plt.axvline(0, color='black', linewidth=0.8)
plt.xlabel('Curvature (1/in)')
plt.ylabel('Moment (kip·in)')
plt.title('Moment–Curvature Diagram (Steel 02-Concrete 02)')
plt.grid(True, alpha=0.6)
plt.tight_layout()
plt.show()
