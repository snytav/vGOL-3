import numpy as np  # import module(s)
import matplotlib.pyplot as plt  # import module(s)
# (blank line)
# Load and visualize the final structural state from your data file
data = np.genfromtxt("gol3_anomalous_conductivity.csv", delimiter=',', names=True)  # assignment/initialization
time = data['time'] * 1e12  # assignment/initialization
max_E = data['max_E']  # assignment/initialization
# (blank line)
plt.figure(figsize=(7, 4))  # assignment/initialization
plt.plot(time, max_E, color='tab:blue', lw=2)  # assignment/initialization
plt.title("Langmuir Microturbulence Saturation Profile")  # statement
plt.xlabel("Time (ps)")  # statement
plt.ylabel("Electric Field Amplitude $E_x$ (V/m)")  # statement
plt.grid(True, ls=":")  # assignment/initialization
plt.tight_layout()  # statement
plt.savefig("gol3_field_evolution.png", dpi=300)  # assignment/initialization
print("Saved turbulence evolution check to: gol3_field_evolution.png")  # print output
