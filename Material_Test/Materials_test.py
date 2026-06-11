# source /Users/niraj/Downloads/Projects/Project_Geo_Lab/venv_arm/bin/activate

import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np

ops.wipe()  

# ------------------------------------------------------------
# Material Properties
# ------------------------------------------------------------
steel_tag = 1

Fy_steel = 36.0                 # Yield stress (ksi)
E0_steel = 29000.0              # Initial modulus (ksi)
Bs = 0.01                       # strain-hardening ratio
params_steel = [20, 0.925, 0.15]  # [R0, cR1, cR2] controls transition curvature

ops.uniaxialMaterial("Steel02", steel_tag, Fy_steel, E0_steel, Bs, *params_steel) 

# ------------------------------------------------------------
# Testing Materials
# ------------------------------------------------------------
ops.testUniaxialMaterial(steel_tag)

# Yield strain and step increment for Option A
yield_strain = Fy_steel / E0_steel 
deps = 0.02 * yield_strain                                      # step increment = 2% of yield strain

# Target peak strain factors
peak_factors = np.array([0, 2, -2, 4, -4, 6, -6, 8, -8, 0])     # Multiples of yield strain
peaks = peak_factors * yield_strain

# Cyclic strain path generation
strain_segments = []
for i in range(len(peaks) - 1):
    start = peaks[i]
    end = peaks[i+1]
    
    num_steps = int(round(abs(end - start) / deps))              # Number of steps needed to maintain 'deps' step-size
    
    if i == 0:
        segment = np.linspace(start, end, num_steps + 1)
    else:
        segment = np.linspace(start, end, num_steps + 1)[1:]     # Cuts off the very first number of every segment
    strain_segments.append(segment)

strain_values_steel = np.concatenate(strain_segments)

stress_steel = []
strain_steel = []

# Stress values for each strain value
for eps in strain_values_steel:
    ops.setStrain(eps)
    stress = ops.getStress()
    strain = ops.getStrain()
    stress_steel.append(stress)
    strain_steel.append(strain)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

strain_percentage = np.array(strain_steel) * 100

# Figure 1: Strain History
plt.figure()
plt.plot(strain_percentage, color='tab:blue')       # x-axis is just the index (0, 1, 2, ...) of each value
plt.grid(True)
plt.title("Strain History")
plt.xlabel("Step")
plt.ylabel("Strain (%)")

# Figure 2: Hysteresis Stress-Strain Loop
plt.figure()
plt.plot(strain_percentage, stress_steel, color='tab:red')
plt.grid(True)
plt.title("Stress-Strain Hysteresis")
plt.xlabel("Strain (%)")
plt.ylabel("Stress (ksi)")

plt.show()