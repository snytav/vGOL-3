import numpy as np  # import module(s)
import matplotlib.pyplot as plt  # import module(s)
# (blank line)
def extract_exact_growth_rate(csv_path):  # function definition
    print("=== AUTOMATED LINEAR PLATEAU DETECTION SYSTEM ===")  # print output
# (blank line)
    # Load the diagnostics file containing the instantaneous gamma states
    # Expected columns: step, time, E_k1_amp, instantaneous_gamma
    try:  # control flow / block statement
        data = np.genfromtxt(csv_path, delimiter=',', names=True)  # assignment/initialization
        time = data['time']  # assignment/initialization
        gamma = data['instantaneous_gamma']  # assignment/initialization
        e_amp = data['E_k1_amp']  # assignment/initialization
    except Exception as e:  # control flow / block statement
        print(f"Error loading file: {e}")  # print output
        print("Please ensure your C++ logger outputs 'time', 'E_k1_amp', and 'instantaneous_gamma'.")  # print output
        return  # return value
# (blank line)
    # 1. Setup a Sliding Window to Scan for the Plateau
    # We look for a region where the standard deviation of gamma is minimized (perfectly flat line)
    window_size = 20  # Number of data points in the sliding check window
    num_points = len(gamma)  # assignment/initialization
# (blank line)
    best_variance = float('inf')  # assignment/initialization
    plateau_start_idx = 0  # assignment/initialization
    plateau_end_idx = 0  # assignment/initialization
# (blank line)
    # Slide the window across the simulation timeline, avoiding the very edges
    for i in range(10, num_points - window_size - 10):  # control flow / block statement
        window_data = gamma[i : i + window_size]  # assignment/initialization
        window_amp = e_amp[i : i + window_size]  # assignment/initialization
# (blank line)
        # Calculate variance in this window
        variance = np.var(window_data)  # assignment/initialization
# (blank line)
        # Physics Guardrails: The plateau must have a positive growth rate,
        # and the wave amplitude must be growing (not noise floor/saturation)
        if variance < best_variance and np.mean(window_data) > 0 and (window_amp[-1] > window_amp[0]):  # control flow / block statement
            best_variance = variance  # assignment/initialization
            plateau_start_idx = i  # assignment/initialization
            plateau_end_idx = i + window_size  # assignment/initialization
# (blank line)
    # 2. Extract and Validate the Plateau Value
    if best_variance == float('inf'):  # control flow / block statement
        print("[ERROR]: No stable linear growth plateau region detected.")  # print output
        return  # return value
# (blank line)
    # Expand window slightly if adjacent points share the same flat profile
    while plateau_start_idx > 0 and abs(gamma[plateau_start_idx-1] - gamma[plateau_start_idx]) < 2 * np.sqrt(best_variance):  # control flow / block statement
        plateau_start_idx -= 1  # assignment/initialization
    while plateau_end_idx < num_points and abs(gamma[plateau_end_idx] - gamma[plateau_end_idx-1]) < 2 * np.sqrt(best_variance):  # control flow / block statement
        plateau_end_idx += 1  # assignment/initialization
# (blank line)
    plateau_time = time[plateau_start_idx:plateau_end_idx]  # assignment/initialization
    plateau_gamma = gamma[plateau_start_idx:plateau_end_idx]  # assignment/initialization
# (blank line)
    # Calculate exact physical quantities
    exact_mean_gamma = np.mean(plateau_gamma)  # assignment/initialization
    gamma_std_err = np.std(plateau_gamma) / np.sqrt(len(plateau_gamma))  # assignment/initialization
# (blank line)
    print(f" Detected Plateau Window : {plateau_time[0]*1e12:.2f} ps to {plateau_time[-1]*1e12:.2f} ps")  # print output
    print(f" Exact Mean Growth Rate  : {exact_mean_gamma:.4e} rad/s")  # print output
    print(f" Statistical Uncertainty : \u00b1 {gamma_std_err:.4e} rad/s")  # print output
    print("====================================================")  # print output
# (blank line)
    # 3. Generate Presentation-Grade Visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)  # assignment/initialization
# (blank line)
    # Plot 1: Wave Amplitude Growth
    ax1.semilogy(time * 1e12, e_amp, label="Fundamental Mode $|E_{k1}|$", color='black', lw=2)  # assignment/initialization
    ax1.axvspan(plateau_time[0]*1e12, plateau_time[-1]*1e12, color='green', alpha=0.15, label='Detected Linear Growth Stage')  # assignment/initialization
    ax1.set_ylabel("Wave Amplitude (V/m)")  # statement
    ax1.grid(True, which="both", ls=":")  # assignment/initialization
    ax1.legend()  # statement
    ax1.set_title("In-Situ Wave Instability Progression and Plateau Analysis")  # statement
# (blank line)
    # Plot 2: Instantaneous Gamma Plateau Verification
    ax2.plot(time * 1e12, gamma, label="Instantaneous $\Gamma(t)$", color='purple', lw=1.5)  # assignment/initialization
    ax2.plot(plateau_time * 1e12, plateau_gamma, color='crimson', lw=3, label='Verified Physics Plateau')  # assignment/initialization
    ax2.axhline(exact_mean_gamma, color='crimson', linestyle='--', alpha=0.7,   # assignment/initialization
               label=f'Mean $\Gamma$ = {exact_mean_gamma:.3e} rad/s')  # assignment/initialization
# (blank line)
    # Dynamic vertical axis limit focusing on the plateau area
    ax2.set_ylim(-0.5 * exact_mean_gamma, 2.0 * exact_mean_gamma)  # statement
    ax2.set_xlabel("Time (ps)")  # statement
    ax2.set_ylabel("Instantaneous Growth Rate (rad/s)")  # statement
    ax2.grid(True, ls=":")  # assignment/initialization
    ax2.legend()  # statement
# (blank line)
    plt.tight_layout()  # statement
    plt.savefig("plateau_detection_results.png", dpi=300)  # assignment/initialization
    print(" Diagnostic chart successfully exported to: plateau_detection_results.png")  # print output
# (blank line)
if __name__ == "__main__":  # control flow / block statement
    # Generate mock format or process your output path file string
    extract_exact_growth_rate("two_stream_diagnostics.csv")  # statement
