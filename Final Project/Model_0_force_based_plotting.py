"""
==============================================================================
OpenSees Finite Element Modeling Project

Date: July 3, 2026
Author: Niraj Yadav

Units: Newtons (N), Millimeters (mm), Seconds (s)

Baseline model (Model 0) implemented with forceBeamColumn elements

Plotting of Results: 
1. Hysteresis Curve
2. Rebar Stress Strain Curve
3. Section Moment-Curvature Curve
==============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import os

output_folder = os.path.join("Final_Project", "Results", "Model_0_force_based")

#-----------------------------------------------------------------------------------------------
# Hysteresis Curve
#-----------------------------------------------------------------------------------------------

disp_data_disp_Beam_Column = np.loadtxt(os.path.join(output_folder, 'col_tip_disp.out'))
force_data_disp_Beam_Column = np.loadtxt(os.path.join(output_folder, 'base_reaction.out'))

tip_disp_disp_Beam_Column  = disp_data_disp_Beam_Column[:, 1]
base_shear_disp_Beam_Column  = -1e-3 * force_data_disp_Beam_Column[:, 1] 

plt.figure(figsize=(9, 7))
plt.plot(tip_disp_disp_Beam_Column, base_shear_disp_Beam_Column, color='steelblue', linewidth=1.5)

# Set major grid spacing: every 50 mm on X, every 20 kN on Y.
plt.gca().xaxis.set_major_locator(MultipleLocator(50))
plt.gca().yaxis.set_major_locator(MultipleLocator(20))

plt.grid(True, linestyle='--', alpha=0.6)
plt.axhline(0, color='black', linewidth=1.0)
plt.axvline(0, color='black', linewidth=1.0)

plt.xlabel('Tip Deflection (mm)')
plt.ylabel('Load P (kN)')
plt.title('Hysteresis Curve')

plt.tight_layout()
# plt.show()  

# Optional: export as a transparent PNG for report use
plt.axis('off')
plt.savefig(os.path.join(output_folder, 'Model_0_force_based_Hysteresis_Curve.png'), 
            transparent=True,       # Makes the background completely transparent
            bbox_inches='tight',    # Crops out excess empty margin space
            pad_inches=0,           # Removes padding
            dpi=400)
plt.show()

#-----------------------------------------------------------------------------------------------
# Top Rebar Stress Strain Curve
# Element 3031 : plastic-hinge element directly adjacent to the column joint (nodes 31 to 30)
# Element 3132 : next element further from the joint in the left beam (nodes 32 to 31)
#-----------------------------------------------------------------------------------------------

elements = [3033]
sections = [1, 2, 3, 4]
colors = {
    (3033, 1): 'steelblue',
    (3033, 2): 'darkorange',
    (3033, 3): 'firebrick',
    (3033, 4): 'forestgreen',
    (3033, 5): 'darkred'
}

for ele in elements:
    for sec in sections:

        data = np.loadtxt(os.path.join(output_folder, f'ele{ele}_sec{sec}_topRebar.out'))

        stress = data[:, 1]
        strain = data[:, 2]

        plt.figure(figsize=(9, 7))
        plt.plot(strain, stress, color=colors[(ele, sec)], linewidth=1.5)

        plt.grid(True, linestyle='--', alpha=0.6)
        plt.axhline(0, color='black', linewidth=1.0)
        plt.axvline(0, color='black', linewidth=1.0)

        plt.xlabel('Strain')
        plt.ylabel('Stress (MPa)')
        plt.title(f'Top Rebar Stress-Strain Curve\nElement {ele}, Section {sec}')

        plt.tight_layout()

# plt.show()

#-----------------------------------------------------------------------------------------------
# section moment-curvature
#-----------------------------------------------------------------------------------------------

for ele in elements:
    for sec in sections:

        def_data = np.loadtxt(os.path.join(output_folder, f"ele{ele}_sec{sec}_def.out"))
        force_data = np.loadtxt(os.path.join(output_folder, f"ele{ele}_sec{sec}_force.out"))

        curvature = 1e3 * def_data[:, 2]      # Curvature: 1/mm to 1/m
        moment = 1e-6 * force_data[:, 2]      # Moment: N.mm to kN.m

        plt.figure(figsize=(9, 7))
        plt.plot(curvature, moment, color=colors[(ele, sec)], linewidth=1.5)

        plt.grid(True, linestyle='--', alpha=0.6)
        plt.axhline(0, color='black')
        plt.axvline(0, color='black')

        plt.xlabel("Curvature (1/m)")
        plt.ylabel("Moment (kN-m)")
        plt.title(f"Moment-Curvature Response\nElement {ele} Section {sec}")

        plt.tight_layout()

plt.show()