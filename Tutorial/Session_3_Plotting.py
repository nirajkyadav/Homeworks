import numpy as np
import matplotlib.pyplot as plt
import os

output_folder = os.path.join("Tutorial", "Session_3_Results")

output_folder_disp_Beam_Column = os.path.join(output_folder, "disp_Beam_Column")
output_folder_force_Beam_Column = os.path.join(output_folder, "force_Beam_Column")

disp_data_disp_Beam_Column = np.loadtxt(os.path.join(output_folder_disp_Beam_Column, "node5disp.out"))
force_data_disp_Beam_Column = np.loadtxt(os.path.join(output_folder_disp_Beam_Column, "node1reaction.out"))

disp_data_force_Beam_Column = np.loadtxt(os.path.join(output_folder_force_Beam_Column, "node5disp.out"))
force_data_force_Beam_Column = np.loadtxt(os.path.join(output_folder_force_Beam_Column, "node1reaction.out"))

tip_disp_disp_Beam_Column  = disp_data_disp_Beam_Column[:, 1]
base_shear_disp_Beam_Column  = -force_data_disp_Beam_Column[:, 1] 

tip_disp_force_Beam_Column  = disp_data_force_Beam_Column[:, 1]
base_shear_force_Beam_Column  = -force_data_force_Beam_Column[:, 1] 

plt.figure(figsize=(9, 7))
plt.plot(tip_disp_disp_Beam_Column, base_shear_disp_Beam_Column, color='steelblue', linewidth=1.2, label = 'dispBeamColumn')
plt.plot(tip_disp_force_Beam_Column, base_shear_force_Beam_Column, color='firebrick', linewidth=1.2, label = 'forceBeamColumn')

plt.grid(True, linestyle='--', alpha=0.6)
plt.axhline(0, color='black', linewidth=1.0)
plt.axvline(0, color='black', linewidth=1.0)

plt.xlabel('Tip Deflection (in)')
plt.ylabel('Applied Load P (kips)')
plt.title('Hysteresis Curve')
plt.legend(loc='lower right', frameon=True, shadow=False)

plt.tight_layout()
plt.show()