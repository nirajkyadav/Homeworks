"""
Moment Curvature Analysis (Cyclic)
Units (in,kip,ksi)
"""

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import opsvis as opsv
import os

ops.wipe()
ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)

# ---------------------------------------------------------------------------------
# Nodes
# ---------------------------------------------------------------------------------

ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, 0.0)

ops.fix(1, 1, 1, 1)
ops.fix(2, 0, 1, 0)

#---------------------------------------------------------------------------------
# Material
#---------------------------------------------------------------------------------

steel_01 = 1
steel_02 = 2
concrete_01 = 3
concrete_02 = 4

ops.uniaxialMaterial('Steel01', steel_01, 60., 29000.0, 0.01) 
ops.uniaxialMaterial('Steel02', steel_02, 60., 29000.0, 0.01, *[18,0.925,0.15]) 
ops.uniaxialMaterial('Concrete01', concrete_01, -6.0, -0.004, -5.0, -0.014)       
ops.uniaxialMaterial('Concrete02', concrete_02, -6.0, -0.004, -5.0, -0.014, 0.1, 0.56, 280) 

#---------------------------------------------------------------------------------
# Section
#---------------------------------------------------------------------------------

# Geometry
col_width = 16.0
col_depth = 16.0
cover = 2.5
As = 1.0
y1 = col_depth/2.0
z1 = col_width/2.0

SecTag = 1

ops.section('Fiber', SecTag)
ops.patch('rect', concrete_02, 9, 1, -y1, -z1, y1, z1)
ops.layer('straight', steel_02, 4, As, y1 - cover, z1 - cover, y1 - cover, cover - z1)
ops.layer('straight', steel_02, 4, As, cover - y1, z1 - cover, cover - y1, cover - z1)

#---------------------------------------------------------------------------------
# Element
#---------------------------------------------------------------------------------

ops.element('zeroLengthSection', 1, 1, 2, SecTag)

#---------------------------------------------------------------------------------
# Loads
#---------------------------------------------------------------------------------

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 2, 1, '-fact', 1.0)
ops.load(2, 0.0, 0.0, 1.0)

#---------------------------------------------------------------------------------
# Recorder
#---------------------------------------------------------------------------------

output_folder = os.path.join("Tutorial", "Session_2_Results")
os.makedirs(output_folder, exist_ok=True)

output_file = os.path.join(output_folder, "Cyclic.out")

ops.recorder('Node', '-file', output_file, '-time', '-node', 2, '-dof', 3, 'disp')

#---------------------------------------------------------------------------------
# Analysis
#---------------------------------------------------------------------------------

ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGeneral')

phiy = 0.0021/(0.6*(col_depth-cover))
numinc = 100

ops.algorithm('Newton')
# ops.test('RelativeNormDispIncr', 1e-4, 50)
ops.test('NormUnbalance', 1e-4, 50)

peaks = [0, 2, -2, 4, -4, 6, -6]
curvature = []
moment    = []

for i in range(len(peaks)-1):
    start = peaks[i]
    end = peaks[i+1]

    cdK0 = (end - start) * phiy / numinc
    ops.integrator('DisplacementControl', 2, 3, cdK0)
    ops.analysis('Static')

    for j in range(numinc):
        ok = ops.analyze(1)    # (0 if successful, <0 if not successful)

        if ok != 0:
            print(f"Analysis failed at increment {j + 1}")
            break

        phi = ops.nodeDisp(2, 3)

        ops.reactions()
        M = ops.nodeReaction(1, 3)
    
        curvature.append(phi)
        moment.append(-M)

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
