import numpy as np  # import module(s)
import matplotlib.pyplot as plt  # import module(s)
# (blank line)
def evaluate_pic_correctness(csv_path):  # function definition
    me = 9.1093837e-31  # assignment/initialization
    qe = -1.602176634e-19  # assignment/initialization
    c = 299792458.0  # assignment/initialization
    B0 = 2.5  # assignment/initialization
# (blank line)
    try:  # control flow / block statement
        data = np.genfromtxt(csv_path, delimiter=',', names=True)  # assignment/initialization
    except Exception as e:  # control flow / block statement
        print(f"Error loading {csv_path}: {e}")  # print output
        return  # return value
# (blank line)
    time = data['time']  # assignment/initialization
    ux = data['tracked_ux']  # assignment/initialization
    uy = data['tracked_uy']  # assignment/initialization
    field_energy = data['field_energy']  # assignment/initialization
# (blank line)
    u_sq = ux**2 + uy**2  # assignment/initialization
    gamma = np.sqrt(1.0 + u_sq / c**2)  # assignment/initialization
    vx = ux / gamma  # assignment/initialization
    vy = uy / gamma  # assignment/initialization
# (blank line)
    dt = time[1] - time[0] if len(time) > 1 else 5e-14  # assignment/initialization
    x_orbit = np.cumsum(vx) * dt  # assignment/initialization
    y_orbit = np.cumsum(vy) * dt  # assignment/initialization
# (blank line)
    measured_radius = (np.max(x_orbit) - np.min(x_orbit)) / 2.0  # assignment/initialization
    v_perp = np.sqrt(vx**2 + vy**2)  # assignment/initialization
    analytical_radius = (gamma * me * v_perp) / (abs(qe) * B0)  # assignment/initialization
# (blank line)
    print("\n" + "="*52)  # print output
    print("   AUTOMATED PHYSICS CORE HEALTH CHECK RESULTS      ")  # print output
    print("="*52)  # print output
    print(f"Analytical Larmor Radius Expected : {analytical_radius[-1] * 1e6:.4f} micrometers")  # print output
    print(f"Simulated Larmor Radius Measured  : {measured_radius * 1e6:.4f} micrometers")  # print output
# (blank line)
    error_margin = abs(measured_radius - analytical_radius[-1]) / analytical_radius[-1] * 100  # assignment/initialization
    print(f"Computed Orbit Tracker Discrepancy: {error_margin:.3f}%")  # print output
# (blank line)
    if error_margin < 5.0:  # control flow / block statement
        print("\n[STATUS: CORE PHYSICS VERIFIED SUCCESSFULLY]")  # print output
    else:  # else branch
        print("\n[STATUS: CORE PHYSICS ERROR ENCOUNTERED]")  # print output
# (blank line)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))  # assignment/initialization
# (blank line)
    ax1.plot(x_orbit * 1e6, y_orbit * 1e6, color='blue', lw=2)  # assignment/initialization
    ax1.set_title("Relativistic Electron Cyclotron Orbit Trajectory")  # statement
    ax1.set_xlabel("X Position Coordinate ($\\mu$m)")  # statement
    ax1.set_ylabel("Y Position Coordinate ($\\mu$m)")  # statement
    ax1.grid(True)  # statement
    ax1.axis('equal')  # statement
# (blank line)
    ax2.plot(ux, uy, color='darkred', linestyle='dashed')  # assignment/initialization
    ax2.set_title("Relativistic Momentum Vector Rotation Phase Space")  # statement
    ax2.set_xlabel("$u_x = \\gamma v_x$ (m/s)")  # assignment/initialization
    ax2.set_ylabel("$u_y = \\gamma v_y$ (m/s)")  # assignment/initialization
    ax2.grid(True)  # statement
    ax2.axis('equal')  # statement
# (blank line)
    ax3.plot(time * 1e9, gamma, color='green', lw=2)  # assignment/initialization
    ax3.set_title("Lorentz Energy Factor ($\\gamma$) Conservation History")  # statement
    ax3.set_xlabel("Time Timeline Execution Metrics (ns)")  # statement
    ax3.set_ylabel("Lorentz Gamma Factor ($\\gamma$)")  # statement
    ax3.grid(True)  # statement
# (blank line)
    ax4.plot(time * 1e9, field_energy, color='purple', lw=2)  # assignment/initialization
    ax4.set_title("Self-Consistent Grid Field Energy Profile")  # statement
    ax4.set_xlabel("Time Timeline Execution Metrics (ns)")  # statement
    ax4.set_ylabel("Total Electrostatic Field Energy (J)")  # statement
    ax4.grid(True)  # statement
# (blank line)
    plt.tight_layout()  # statement
    output_png = "pic_physics_verification.png"  # assignment/initialization
    plt.savefig(output_png, dpi=300)  # assignment/initialization
    print(f"\nVisualization exported successfully to: {output_png}")  # print output
# (blank line)
if __name__ == '__main__':  # control flow / block statement
    evaluate_pic_correctness('pic_diagnostics.csv')  # statement
