import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import opsvis as opsv
import os

ops.wipe()
ops.model('BasicBuilder', '-ndm', 2, '-ndf', 3)

#---------------------------------------------------------------------------------
# Materials
#---------------------------------------------------------------------------------

unconfined_concrete_tag = 1
confined_concrete_tag = 2
steel_tag = 3

fpc    = -5.03
fpcu   = 0.2 * fpc       
epsc0  = -0.002
epsU   = -0.004
fpcc   = -5.73
fpccu  = 0.2 * fpcc      
epscc0 = -0.003
epscU  = -0.021
ft     = 0.28
Ec     = 4300.0
Ets    = 0.1 * Ec       
Fy     = 60.0
Es     = 29000.0

# uniaxialMaterial('Concrete02', matTag, fpc, epsc0, fpcu, epsU, lambda, ft, Ets)
ops.uniaxialMaterial('Concrete02', unconfined_concrete_tag,  fpc,  epsc0,  fpcu,  epsU,   0.1,  ft,  Ets)     # Unconfined Concrete
ops.uniaxialMaterial('Concrete02', confined_concrete_tag,    fpcc, epscc0, fpccu, epscU,  0.1,  ft,  Ets)     # Confined Concrete

# uniaxialMaterial('Steel02', matTag, Fy, E, b, R0, cR1, cR2, a1, a2, a3, a4)
ops.uniaxialMaterial('Steel02', steel_tag,  Fy,  Es,  0.1, 20, 0.925, 0.15, 0.005, 1.0, 0.005, 1.0)           # Steel

#---------------------------------------------------------------------------------
# Section
#---------------------------------------------------------------------------------

colWidth = 15
colDepth = 29
dia      = 1.128
As       = 1
c        = 2.5
sv       = 3          # vertical spacing between bars

y1 = colDepth / 2.0
y2 = y1 - c - dia / 2
y3 = y2 - sv
z1 = colWidth / 2.0
z2 = z1 - c - dia / 2

ops.section('Fiber', 1)

ops.patch('rect', unconfined_concrete_tag,  2,  1,   -y1, -z1,  -y2,  z1)       # bottom cover
ops.patch('rect', unconfined_concrete_tag,  2,  1,    y2, -z1,   y1,  z1)       # top cover
ops.patch('rect', confined_concrete_tag, 15,  1,   -y2, -z1,   y2,  z1)         # core

ops.layer('straight', steel_tag,  4,  As,   -y2,  -z2,  -y2,   z2)              # bottom bars
ops.layer('straight', steel_tag,  2,  As,   -y3,  -z2,  -y3,   z2)              # lower-mid bars
ops.layer('straight', steel_tag,  2,  As,    y3,  -z2,   y3,   z2)              # upper-mid bars
ops.layer('straight', steel_tag,  4,  As,    y2,  -z2,   y2,   z2)              # top bars

#---------------------------------------------------------------------------------
# Nodes And Elements
#---------------------------------------------------------------------------------

L  = 70. / 3
Lr = 8.0

ops.node(1,   0.0,        0.0)
ops.node(2,   L,          0.0)
ops.node(3,   2 * L,      0.0)
ops.node(4,   3 * L,      0.0)
ops.node(5,   3 * L + Lr, 0.0)

ops.fix(1,  1, 1, 1)

ops.geomTransf('Linear', 1)

no_of_integration_points = 7
ops.beamIntegration('Lobatto', 1, 1, no_of_integration_points)

# element('dispBeamColumn', eleTag, *eleNodes, transfTag, integrationTag, '-cMass', '-mass', mass=0.0)
ops.element('dispBeamColumn', 1, 1, 2, 1, 1)
ops.element('dispBeamColumn', 2, 2, 3, 1, 1)
ops.element('dispBeamColumn', 3, 3, 4, 1, 1)

ops.rigidLink('beam', 5, 4)           # rigidLink(type, rNodeTag, cNodeTag)
"""
# ops.rigidLink('beam', 4, 5)
This makes, Node 4 a master node and Node 5 the slave node. With the Transformation constraint handler, the DOFs for the slave node are eliminated from the global system matrix. So, DisplacementControl applied on respnode (Node 5), is internally transferred to the master node (Node 4). Thus, Node 4 moves exactly with step size (dK), but because of the rigid body rotation (θ) of the arm, Node 5 undergoes additional displacement (u_5 = u_4 + L_r.θ), yielding greater current displacement than target displacement. 
Solution: we make Node 5 master, and Node 4 the slave.
"""

#---------------------------------------------------------------------------------
# Recorders
#---------------------------------------------------------------------------------

output_folder = os.path.join("Tutorial", "Session_3_Results", "disp_Beam_Column")
os.makedirs(output_folder, exist_ok=True)

ops.recorder('Node', '-file', os.path.join(output_folder, 'node5disp.out'),     '-time', '-node', 5, '-dof', 2, 'disp')
ops.recorder('Node', '-file', os.path.join(output_folder, 'node1reaction.out'), '-time', '-node', 1, '-dof', 2, 'reaction')

#---------------------------------------------------------------------------------
# Cyclic Loading
#---------------------------------------------------------------------------------

file_path = "peakdisp.txt"

with open(file_path, 'r') as fp:
    peakDisp = [float(line.strip()) for line in fp if line.strip()]

print(peakDisp)

respnode = 5

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

ops.load(5, 0.0, 1.0, 0.0)

# Assumption of yield curvature
epsy = 60.0 / 29000.0
dK0  = 0.025 * epsy                      # increment

#---------------------------------------------------------------------------------
# Analysis
#---------------------------------------------------------------------------------

currentDisp = 0.0
ok          = 0.0

ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('BandGeneral')

for i in range(len(peakDisp)):

    if ok == 0:
        cycleDisp = peakDisp[i] - currentDisp

        sign = 1 if cycleDisp > 0 else -1

        dK = dK0 * sign

        ops.integrator('DisplacementControl', respnode, 2, dK)
        # ops.test('RelativeNormDispIncr', 1.0e-5, 30)
        ops.test('NormUnbalance', 1e-4, 50)
        ops.algorithm('Newton')
        ops.analysis('Static')

        numiter = round(abs(cycleDisp / dK))

        print(f"---> Running {numiter} steps with dK={dK:.6f} to go from {currentDisp:.6f} to {peakDisp[i]:.6f}")

        ok1 = ops.analyze(numiter)

        currentDisp = ops.nodeDisp(respnode, 2)

        if ok1 != 0:
            ok = 0
            print(f"Try stuff, peak disp  = {peakDisp[i]}")
            print(f"Current disp          = {currentDisp}")
            print(f"Cyclic disp           = {cycleDisp}")
            counter = 1

            # while abs(currentDisp) < abs(peakDisp[i]) and ok == 0:
            while (peakDisp[i] - currentDisp) * sign > 0 and ok == 0:
                ok = 1
                
                while ok != 0:
                    if counter == 1:
                        # increase load stepsize
                        dK = dK0 * sign * 1.5
                        print(f"dK = {dK0}*{sign}*1.5 = {dK}")
                        counter = 2
                    elif counter == 2:
                        # increase load stepsize
                        dK = dK0 * sign * 2.00
                        print(f"dK = {dK0}*{sign}*2.0 = {dK}")
                        counter = 3
                    elif counter == 3:
                        # decrease load stepsize
                        dK = dK0 * sign * 0.5
                        print(f"dK = {dK0}*{sign}*0.5 = {dK}")
                        counter = 4
                    elif counter == 4:
                        # decrease load stepsize
                        dK = dK0 * sign * 0.1
                        print(f"dK = {dK0}*{sign}*0.1 = {dK}")
                        counter = 5
                    elif counter == 5:
                        # decrease load stepsize
                        dK = dK0 * sign * 0.05
                        print(f"dK = {dK0}*{sign}*0.05 = {dK}")
                        counter = 6
                    elif counter == 6:
                        # decrease load stepsize
                        dK = dK0 * sign * 0.01
                        print(f"dK = {dK0}*{sign}*0.01 = {dK}")
                        counter = 7
                    elif counter == 7:
                        # decrease load stepsize
                        dK = dK0 * sign * 0.001
                        print(f"dK = {dK0}*{sign}*0.001 = {dK}")
                        counter = 8

                    ops.integrator('DisplacementControl', respnode, 2, dK)
                    ops.analysis('Static')
                    ok = ops.analyze(1)

                    if ok != 0:
                        print("Test Relative Force")
                        ops.test('RelativeNormUnbalance', 1.0e-4, 50, 0)
                        ops.analysis('Static')
                        ok = ops.analyze(1)

                    if ok != 0:
                        print("Test Relative Displacement")
                        ops.test('RelativeNormDispIncr', 1.0e-4, 50, 0)
                        ops.analysis('Static')
                        ok = ops.analyze(1)

                    if ok != 0:
                        print("ModifiedNewton")
                        ops.algorithm('ModifiedNewton')
                        ops.analysis('Static')
                        ok = ops.analyze(1)

                    if ok != 0:
                        print("Newton Modified -initial")
                        ops.algorithm('ModifiedNewton', '-initial')
                        ops.analysis('Static')
                        ok = ops.analyze(1)

                    if ok != 0:
                        print("Test Relative Force")
                        ops.algorithm('Newton', '-initial')
                        ops.test('RelativeNormUnbalance', 1.0e-3, 100, 0)
                        ops.analysis('Static')
                        ok = ops.analyze(1)

                    if ok != 0:
                        print("Test Relative Displ")
                        ops.test('RelativeNormDispIncr', 1.0e-3, 100, 0)
                        ops.analysis('Static')
                        ok = ops.analyze(1)

                    if ok != 0:
                        print("Broyden")
                        ops.algorithm('Broyden', 8)
                        ops.analysis('Static')
                        ok = ops.analyze(1)

                    if ok != 0:
                        print("Newton Line Search")
                        ops.algorithm('NewtonLineSearch', 0.8)
                        ops.analysis('Static')
                        ok = ops.analyze(1)

                    if ok != 0:
                        print("BFGS")
                        ops.algorithm('BFGS')
                        ops.analysis('Static')
                        ok = ops.analyze(1)

                    if ok != 0 and counter == 8:
                        print("Continue to the next time step... ")
                        print(f"Check the results at {currentDisp}")
                        ops.test('RelativeEnergyIncr', 1.0e-3, 100, 5)
                        ops.analysis('Static')
                        ok = ops.analyze(1)

                # Reset counter and update current displacement after a successful sub-step
                counter = 0
                currentDisp = ops.nodeDisp(respnode, 2)
        
        else:
            # Executes if the main analysis step succeeded on the first try
            print(f"target disp  = {peakDisp[i]}")
            print(f"current disp = {ops.nodeDisp(respnode, 2)}")
            

    if ok < 0:
        print("Unsuccessful convergence!")

