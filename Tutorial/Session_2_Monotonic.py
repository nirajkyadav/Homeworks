"""
================================================================================
OpenSees Tutorial - Session 2: Moment Curvature Analysis
Date: 24 June 2026
Author: Niraj Yadav

Units: in, kip, ksi

Assignment:
1) Monotonic loading
   Develop moment-curvature diagram using Displacement control integrator
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
col_width = 16.0               # column width
col_depth = 16.0               # column depth
cover = 2.5                    # cover
As = 1.0                       # area of steel
y1 = col_depth/2.0             # distance from neutral axis to top/bottom face
z1 = col_width/2.0             # distance from neutral axis to left/right face

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
ops.recorder('Node', '-file', os.path.join(output_folder, "Monotonic_Disp_Controlled.out"), '-time', '-node', 2, '-dof', 3, 'disp')

#---------------------------------------------------------------------------------
# Analysis
#---------------------------------------------------------------------------------

ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGeneral')

phiy  = 0.0021 / (0.6 * (col_depth - cover))  # Yield curvature (1/in): εy / (0.6 * effective depth)
numinc = 1000                                 # Number of displacement-controlled increments
dK0   = (5 * phiy) / numinc                   # Curvature increment per step

# integrator(  'DisplacementControl', nodeTag, dof, incr)
ops.integrator('DisplacementControl',       2,   3,  dK0)     

ops.algorithm('Newton')
ops.test('RelativeNormDispIncr', 1e-4, 50)

ops.analysis('Static')

curvature = []   # List to store curvature values (rotation at node 2, DOF 3) [1/in]
moment    = []   # List to store moment values (reaction at node 1, DOF 3) [kip·in]

for i in range(numinc):
    """ Step through each displacement increment and collect M-φ data """
    ok = ops.analyze(1)    # (0 if successful, <0 if not successful)

    if ok != 0:
        print(f"Analysis failed at increment {i + 1}")
        break

    phi = ops.nodeDisp(2, 3)   # Curvature: rotation (DOF 3) at free node 2 [1/in]

    ops.reactions()            # Compute reactions at fixed nodes
    M = ops.nodeReaction(1, 3) # Moment reaction at fixed node 1, DOF 3 [kip·in]

    curvature.append(phi)
    moment.append(abs(M))      # Store absolute value (sign convention: positive moment)
 
# ---------------------------------------------------------------------------------
# Moment Curvature Diagram
# ---------------------------------------------------------------------------------

plt.figure(figsize=(8, 6))
plt.plot(curvature, moment)
plt.xlabel('Curvature (1/in)')
plt.ylabel('Moment (kip·in)')
plt.title('Moment–Curvature Diagram (Monotonic, Displacement Controlled)')
plt.grid(True, alpha=0.6)
plt.tight_layout()
plt.show()
