import numpy as np  # import module(s)
import matplotlib.pyplot as plt  # import module(s)
# (blank line)
def evaluate_gol3_heat_suppression(csv_path):  # function definition
    print("=== PROCESSING GOL-3 TURBULENT HEAT SUPPRESSION ===")  # print output
# (blank line)
    # Load nonlinear telemetry matrix arrays
    data = np.genfromtxt(csv_path, delimiter=',', names=True)  # assignment/initialization
    time = data['time'] * 1e12 # Convert to picoseconds
    max_E = data['max_E']  # assignment/initialization
    q_parallel = data['avg_q_parallel']  # assignment/initialization
# (blank line)
    # Locate the peak saturation time of Langmuir turbulence
    saturation_idx = np.argmax(max_E)  # assignment/initialization
    t_sat = time[saturation_idx]  # assignment/initialization
# (blank line)
    # Evaluate anomalous heat flux suppression rates
    initial_flux = q_parallel[1] # Early linear stage flux
    min_flux_post_sat = np.min(q_parallel[saturation_idx:]) # Post-saturation compressed flux
# (blank line)
    measured_suppression = initial_flux / (min_flux_post_sat + 1e-30)  # assignment/initialization
# (blank line)
    print(f" Microturbulence Peak Saturation Point: {t_sat:.2f} ps")  # print output
    print(f" Initial Classical Core Heat Flux     : {initial_flux:.4e} W/m^2")  # print output
    print(f" Post-Turbulence Suppressed Heat Flux : {min_flux_post_sat:.4e} W/m^2")  # print output
    print(f" Measured Parallel Conductivity Drop  : {measured_suppression:.2f}-fold")  # print output
    print("====================================================")  # print output
# (blank line)
    # Render presentation graphics chart
    fig, ax1 = plt.subplots(figsize=(8, 5))  # assignment/initialization
# (blank line)
    color = 'tab:blue'  # assignment/initialization
    ax1.set_xlabel('Time (ps)')  # statement
    ax1.set_ylabel('Peak Langmuir Turbulence Field $E_x$ (V/m)', color=color)  # assignment/initialization
    ax1.plot(time, max_E, color=color, lw=2, label="Electric Field Amplitude")  # assignment/initialization
    ax1.tick_params(axis='y', labelcolor=color)  # assignment/initialization
    ax1.grid(True, ls=":")  # assignment/initialization
# (blank line)
    ax2 = ax1.twinx()    # assignment/initialization
    color = 'tab:red'  # assignment/initialization
    ax2.set_ylabel('In-Situ Parallel Heat Flux $q_\parallel$ ($W/m^2$)', color=color)  # assignment/initialization
    ax2.plot(time, q_parallel, color=color, lw=2, linestyle='--', label="Parallel Heat Flux")  # assignment/initialization
    ax2.tick_params(axis='y', labelcolor=color)  # assignment/initialization
# (blank line)
    plt.title("GOL-3 Simulation: Microturbulence Growth vs Axial Heat Flux Suppression")  # statement
    fig.tight_layout()  # statement
    plt.savefig("gol3_nonlinear_suppression.png", dpi=300)  # assignment/initialization
    plt.show()  # statement
# (blank line)
if __name__ == "__main__":  # control flow / block statement
    evaluate_gol3_heat_suppression("gol3_anomalous_conductivity.csv")  # statement
