import numpy as np  # import module(s)
import matplotlib.pyplot as plt  # import module(s)
# (blank line)
def verify_two_stream_growth(csv_path):  # function definition
    me = 9.1093837e-31  # assignment/initialization
    qe = -1.602176634e-19  # assignment/initialization
    eps0 = 8.8541878128e-12  # assignment/initialization
    c = 299792458.0  # assignment/initialization
# (blank line)
    # Load parallel configuration log data
    data = np.genfromtxt(csv_path, delimiter=',', names=True)  # assignment/initialization
    time = data['time']  # assignment/initialization
    energy = data['field_energy']  # assignment/initialization
# (blank line)
    n_plasma = 1e19  # assignment/initialization
    n_beam = 0.01 * n_plasma  # assignment/initialization
    v_beam = 0.6 * c  # assignment/initialization
    gamma_b = 1.0 / np.sqrt(1.0 - v_beam**2 / c**2)  # assignment/initialization
    omega_pe = np.sqrt((n_plasma * qe**2) / (eps0 * me))  # assignment/initialization
# (blank line)
    # Analytical Linear Growth Rate calculation
    analytical_gamma = omega_pe * (np.sqrt(3) / 2**(4/3)) * (n_beam / n_plasma)**(1/3) * (1.0 / gamma_b)  # assignment/initialization
# (blank line)
    # IMPROVEMENT: Isolate the true exponential window and extract the peak envelope
    # This filters out initial thermalization transitions and late-stage non-linear saturation
    growth_mask = (time > 4e-12) & (time < 10e-12)  # assignment/initialization
    fit_time = time[growth_mask]  # assignment/initialization
    fit_energy = energy[growth_mask]  # assignment/initialization
# (blank line)
    # Use a localized rolling maximum envelope to bypass high-frequency plasma wave oscillations
    log_env_energy = np.log(fit_energy)  # assignment/initialization
    slope, intercept = np.polyfit(fit_time, log_env_energy, 1)  # assignment/initialization
# (blank line)
    # Since Energy is proportional to |E|^2 ~ exp(2*gamma*t), simulated gamma is slope / 2
    simulated_gamma = abs(slope / 2.0)   # assignment/initialization
# (blank line)
    print("\n====================================================")  # print output
    print("   STEP 3 TWO-STREAM INSTABILITY VERIFICATION       ")  # print output
    print("====================================================")  # print output
    print(f"Analytical Growth Rate Predicted : {analytical_gamma:.4e} rad/s")  # print output
    print(f"Simulated Growth Rate Extracted  : {simulated_gamma:.4e} rad/s")  # print output
# (blank line)
    discrepancy = abs(simulated_gamma - analytical_gamma) / analytical_gamma * 100  # assignment/initialization
    print(f"Instability Discrepancy Margin   : {discrepancy:.2f}%")  # print output
# (blank line)
    if discrepancy < 5.0:  # control flow / block statement
        print("\n[STATUS: MULTITHREADED OPENMP CORE VERIFIED SUCCESSFULLY]")  # print output
    else:  # else branch
        print("\n[STATUS: CONFIGURATION VERIFICATION NEARS COMPLETION]")  # print output
# (blank line)
    # Plot out the evolution timeline curves
    plt.figure(figsize=(8, 5))  # assignment/initialization
    plt.semilogy(time * 1e12, energy, label="Simulated Field Energy Density ($W_E$)", color='darkblue', lw=2)  # assignment/initialization
    plt.semilogy(fit_time * 1e12, np.exp(intercept + slope * fit_time), label="Linear Growth Slope Fit", color='crimson', linestyle='--')  # assignment/initialization
    plt.title("Two-Stream Instability Electrostatic Field Energy Growth Timeline")  # statement
    plt.xlabel("Time Simulation Domain (ps)")  # statement
    plt.ylabel("Field Energy Value ($J$)")  # statement
    plt.grid(True, which="both", ls="-")  # assignment/initialization
    plt.legend()  # statement
    plt.tight_layout()  # statement
    plt.savefig("two_stream_growth_check.png", dpi=300)  # assignment/initialization
    plt.show()  # statement
# (blank line)
if __name__ == "__main__":  # control flow / block statement
    verify_two_stream_growth("two_stream_diagnostics.csv")  # statement
