"""
==============================================================================
OpenSees Finite Element Modeling Project

Date: July 2, 2026
Author: Niraj Yadav

Units: Newtons (N), Millimeters (mm), Seconds (s)

Comparison of results for Model 0 vs Model 0 with doubled plastic hinge elements length
==============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import os

output_folder_Model_0 = os.path.join("Final_Project", "Results", "Model_0")
output_folder_Model_0_Doubled = os.path.join("Final_Project", "Results", "Model_0_Doubled")

#-----------------------------------------------------------------------------------------------
# Hysteresis Curve
#-----------------------------------------------------------------------------------------------

disp_data_Model_0 = np.loadtxt(os.path.join(output_folder_Model_0, 'col_tip_disp.out'))
force_data_Model_0 = np.loadtxt(os.path.join(output_folder_Model_0, 'base_reaction.out'))

disp_data_Model_0_Doubled = np.loadtxt(os.path.join(output_folder_Model_0_Doubled, 'col_tip_disp.out'))
force_data_Model_0_Doubled = np.loadtxt(os.path.join(output_folder_Model_0_Doubled, 'base_reaction.out'))

tip_disp_Model_0  = disp_data_Model_0[:, 1]
base_shear_Model_0  = -1e-3 * force_data_Model_0[:, 1] 
tip_disp_Model_0_Doubled  = disp_data_Model_0_Doubled[:, 1]
base_shear_Model_0_Doubled  = -1e-3 * force_data_Model_0_Doubled[:, 1] 

plt.figure(figsize=(9, 7))
plt.plot(tip_disp_Model_0, base_shear_Model_0, color='steelblue', linewidth=1.5, label='Model 0' )
plt.plot(tip_disp_Model_0_Doubled, base_shear_Model_0_Doubled, color='darkorange', linewidth=1.5, label='Model 0 With Doubled Plastic Hinge Elements Length')

# Set major grid spacing: every 50 mm on X, every 20 kN on Y.
plt.gca().xaxis.set_major_locator(MultipleLocator(50))
plt.gca().yaxis.set_major_locator(MultipleLocator(20))

plt.grid(True, linestyle='--', alpha=0.6)
plt.axhline(0, color='black', linewidth=1.0)
plt.axvline(0, color='black', linewidth=1.0)
plt.legend()
plt.xlabel('Tip Deflection (mm)')
plt.ylabel('Load P (kN)')
plt.title('Hysteresis Curve')

plt.tight_layout()
plt.show()  