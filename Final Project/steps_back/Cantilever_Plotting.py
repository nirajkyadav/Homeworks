"""
================================================================================
Date: 15 July 2026
Author: Niraj Yadav

Units: N, mm, s

Cantilever with unsymmetrical rectangular section plotting
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import os

output_folder = os.path.join("Final_Project", "steps_back", "Results")

disp = np.loadtxt(os.path.join(output_folder, "node5disp.out"))
force = np.loadtxt(os.path.join(output_folder, "node1reaction.out"))

tip_disp  = disp[:, 1]
base_shear  = -1e-3 * force[:, 1] 

plt.figure(figsize=(9, 7))
plt.plot(tip_disp, base_shear, color='steelblue', linewidth=1.2)

plt.grid(True, linestyle='--', alpha=0.6)
plt.axhline(0, color='black', linewidth=1.0)
plt.axvline(0, color='black', linewidth=1.0)

plt.xlabel('Tip Deflection (mm)')
plt.ylabel('Applied Load P (kN)')
plt.title('Hysteresis Curve')

plt.tight_layout()
plt.show()