"""
================================================================================
Date: 15 July 2026
Author: Niraj Yadav

Units: N, mm, s

Monotonic loading: Moment-curvature diagram of an unsymmetrical section

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
# Recorders
#---------------------------------------------------------------------------------

output_folder = os.path.join("Final_Project", "steps_back", "Results")  # Output directory for results
os.makedirs(output_folder, exist_ok=True)                               # Create folder if it does not exist

# Record rotation (DOF 3) at node 2 at every step
ops.recorder('Node', '-file', os.path.join(output_folder, "Moment_Curvature_Monotonic.out"), '-time', '-node', 2, '-dof', 3, 'disp')

# Record extreme rebar layers to ensure tracking regardless of bending direction
ops.recorder('Element', '-file', os.path.join(output_folder, 'sec_topRebar.out'), '-time', '-ele', 1, 'section', 1, 'fiber', 186.5, 72.5, steel_beam_tag, 'stressStrain')
ops.recorder('Element', '-file', os.path.join(output_folder, 'sec_botRebar.out'), '-time', '-ele', 1, 'section', 1, 'fiber', -186.5, 72.5, steel_beam_tag, 'stressStrain')

#---------------------------------------------------------------------------------
# Analysis Setup
#---------------------------------------------------------------------------------

ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGeneral')

yield_strain = 295.0 / 210_000.0                               # Fy/E0 (0.00140476)
phiy_approx  = yield_strain / (0.6 * (beamDepth - beamCover))  # Empirical yield curvature approximation
numinc       = 1000                                            # Number of displacement-controlled increments
dK0          = (5 * phiy_approx) / numinc                      # Curvature increment per step

# integrator(  'DisplacementControl', nodeTag, dof, incr)
ops.integrator('DisplacementControl',       2,   3,  dK0)     

ops.algorithm('Newton')
ops.test('RelativeNormDispIncr', 1e-4, 50)
ops.analysis('Static')

curvature = []   # List to store curvature values (rotation at node 2, DOF 3) 
moment    = []   # List to store moment values (reaction at node 1, DOF 3)

yielded = False
phi_y_actual = None
M_y_actual = None

#---------------------------------------------------------------------------------
# Analysis Loop
#---------------------------------------------------------------------------------

for i in range(numinc):
    """ Step through each displacement increment and collect M-φ data """
    ok = ops.analyze(1)    # (0 if successful, <0 if not successful)

    if ok != 0:
        print(f"Analysis failed at increment {i + 1}")
        break

    phi = ops.nodeDisp(2, 3)   # Curvature: rotation (DOF 3) at free node 2 
    ops.reactions()            # Compute reactions at fixed nodes
    M = ops.nodeReaction(1, 3) # Moment reaction at fixed node 1, DOF 3 

    current_moment_kNm = abs(M * 1e-6)
    curvature.append(phi * 1e3)
    moment.append(current_moment_kNm)      

    # Programmatic real-time check for actual steel yielding
    if not yielded:
        try:
            # Extract strain values from both extreme rebar layers
            strain_top = abs(ops.eleResponse(1, 'section', 1, 'fiber', 186.5, 72.5, steel_beam_tag, 'stressStrain')[1])
            strain_bot = abs(ops.eleResponse(1, 'section', 1, 'fiber', -186.5, 72.5, steel_beam_tag, 'stressStrain')[1])
            
            if strain_top >= yield_strain or strain_bot >= yield_strain:
                yielded = True
                phi_y_actual = phi
                M_y_actual = current_moment_kNm
        except:
            # Fallback if eleResponse syntax varies across platform builds
            pass

if yielded and phi_y_actual is not None:
    print(f" Empirical phiy reference: {phiy_approx * 1e3:.6e} 1/m")
    print(f" Actual Yield Curvature:   {phi_y_actual * 1e3:.6e} 1/m")
    print(f" Actual Yield Moment:      {M_y_actual:.2f} kN·m")
else:
    print(" WARNING: Yield strain was not reached in the analysis range.")

# ---------------------------------------------------------------------------------
# Moment Curvature Diagram
# ---------------------------------------------------------------------------------

plt.figure(figsize=(8, 6))
plt.plot(curvature, moment)
plt.xlabel('Curvature (1/m)')
plt.ylabel('Moment (kN·m)')
plt.title('Moment–Curvature Diagram')
plt.grid(True, alpha=0.6)
plt.tight_layout()
plt.show()
