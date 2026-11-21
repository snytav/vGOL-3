#include <iostream>  // include iostream
#include <vector>  // include vector
#include <cmath>  // include cmath
#include <fstream>  // include fstream
#include <omp.h>  // include omp.h
#include <algorithm>  // include algorithm
#include <random>  // include random
#include <mpi.h>  // include mpi.h
// (blank line)
struct Particle { double x, ux, uy, uz, w; int id; int species; };   // type definition
constexpr double eps0 = 8.8541878128e-12, me = 9.1093837e-31, qe = -1.602176634e-19, c = 299792458.0;  // variable assignment/initialization
// (blank line)
int main(int argc, char* argv[]) {  // begin block
    MPI_Init(&argc, &argv);  // statement
    int mpi_rank, mpi_size;  // statement
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);  // statement
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);  // statement
// (blank line)
    const int Nx_global = 512, ppc_bg = 150, ppc_bm = 50, num_steps = 1500;  // variable assignment/initialization
    const double dx = 2e-5, L_global = Nx_global * dx, dt = 1e-14;  // variable assignment/initialization
// (blank line)
    const int Nx_local = Nx_global / mpi_size;  // variable assignment/initialization
    const double x_local_min = mpi_rank * Nx_local * dx;  // variable assignment/initialization
    const double x_local_max = (mpi_rank + 1) * Nx_local * dx;  // function declaration/call or variable initialization
// (blank line)
    std::vector<double> Ex(Nx_local + 4, 0.0), rho(Nx_local + 4, 0.0); // Expanded padding registers 
    std::vector<double> local_q_parallel(Nx_local + 1, 0.0);  // function declaration/call or variable initialization
    std::vector<Particle> local_particles;  // statement
// (blank line)
    const double n_plasma = 1e19, n_beam = 0.01 * n_plasma;  // variable assignment/initialization
    const double w_p = (n_plasma * L_global) / (Nx_global * ppc_bg);  // function declaration/call or variable initialization
    const double w_b = (n_beam * L_global) / (Nx_global * ppc_bm);  // function declaration/call or variable initialization
// (blank line)
    std::mt19937 gen(42 + mpi_rank);  // function declaration/call or variable initialization
    std::normal_distribution<double> dist_v(0.0, 0.01 * c);  // function declaration/call or variable initialization
// (blank line)
    for (int i = 0; i < Nx_local * ppc_bg; ++i) {  // begin block
        Particle p; p.x = x_local_min + ((i / ppc_bg) + 0.5) / ppc_bg * dx; p.w = w_p;  // variable assignment/initialization
        p.ux = dist_v(gen); p.uy = dist_v(gen); p.uz = dist_v(gen); p.species = 0; p.id = i;  // variable assignment/initialization
        local_particles.push_back(p);  // statement
    }  // end block
// (blank line)
    double v_beam = 0.6 * c, gamma_b = 1.0 / std::sqrt(1.0 - v_beam*v_beam/(c*c));  // function declaration/call or variable initialization
    for (int i = 0; i < Nx_local * ppc_bm; ++i) {  // begin block
        Particle p; p.x = x_local_min + ((i / ppc_bm) + 0.5) / ppc_bm * dx; p.w = w_b;  // variable assignment/initialization
        p.ux = gamma_b * v_beam; p.uy = 0.0; p.uz = 0.0; p.species = 1; p.id = i + 100000;  // variable assignment/initialization
// (blank line)
        // Apply perturbation seed with a hard global wrap safety catch
        p.x += 0.01 * dx * std::sin(2.0 * M_PI * p.x / L_global);  // variable assignment/initialization
        while (p.x >= L_global) p.x -= L_global;  // control flow statement
        while (p.x < 0.0) p.x += L_global;  // control flow statement
// (blank line)
        // Route only particles belonging to this local domain's footprint
        if (p.x >= x_local_min && p.x < x_local_max) {  // begin block
            local_particles.push_back(p);  // statement
        }  // end block
    }  // end block
// (blank line)
    double* rho_ptr = rho.data();  // function declaration/call or variable initialization
    int left_neighbor  = (mpi_rank == 0) ? mpi_size - 1 : mpi_rank - 1;  // function declaration/call or variable initialization
    int right_neighbor = (mpi_rank == mpi_size - 1) ? 0 : mpi_rank + 1;  // function declaration/call or variable initialization
// (blank line)
    std::ofstream diag_out;  // statement
    if(mpi_rank == 0) {  // begin block
        diag_out.open("gol3_anomalous_conductivity.csv");  // statement
        diag_out << "step,time,max_E,avg_q_parallel\n";  // statement
    }  // end block
// (blank line)
    for (int step = 1; step <= num_steps; ++step) {  // begin block
        std::fill(rho.begin(), rho.end(), 0.0);  // function declaration/call or variable initialization
// (blank line)
        // --- STEP A: Charge Deposition ---
        #pragma omp parallel for reduction(+:rho_ptr[:Nx_local+4]) schedule(static)  // preprocessor directive
        for (size_t i = 0; i < local_particles.size(); ++i) {  // begin block
            double x_norm = (local_particles[i].x - x_local_min) / dx;  // function declaration/call or variable initialization
            int j = static_cast<int>(std::floor(x_norm));   // function declaration/call or variable initialization
// (blank line)
            // Safety Clamp to shield memory footprint access boundaries
            j = std::max(0, std::min(j, Nx_local - 1));  // variable assignment/initialization
            double xi = x_norm - j;  // variable assignment/initialization
// (blank line)
            rho_ptr[j + 1] += (qe * local_particles[i].w * (1.0 - xi)) / dx;  // variable assignment/initialization
            rho_ptr[j + 2] += (qe * local_particles[i].w * xi) / dx;  // variable assignment/initialization
        }  // end block
// (blank line)
        // --- STEP B: Ghost Node Exchange ---
        MPI_Status status;   // statement
        double send_left_ghost = rho[1];   // variable assignment/initialization
        double recv_right_ghost = 0.0;  // variable assignment/initialization
        MPI_Sendrecv(&send_left_ghost, 1, MPI_DOUBLE, left_neighbor, 0, &recv_right_ghost, 1, MPI_DOUBLE, right_neighbor, 0, MPI_COMM_WORLD, &status);  // statement
        rho[Nx_local + 1] += recv_right_ghost;  // variable assignment/initialization
// (blank line)
        double local_sum = 0.0, global_sum = 0.0;  // variable assignment/initialization
        for(int j = 1; j <= Nx_local; ++j) local_sum += rho[j];  // control flow statement
        MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  // statement
        double ion_correction = -global_sum / Nx_global;  // variable assignment/initialization
        for(int j = 1; j <= Nx_local + 1; ++j) rho[j] += ion_correction;  // control flow statement
// (blank line)
        // --- STEP C: Poisson Matrix Integration ---
        std::fill(Ex.begin(), Ex.end(), 0.0);  // function declaration/call or variable initialization
        for (int j = 1; j <= Nx_local + 1; ++j) { Ex[j] = Ex[j - 1] + (rho[j] * dx) / eps0; }  // control flow statement
// (blank line)
        // --- STEP D: Boris-Vay Particle Push ---
        std::vector<Particle> keep_particles; std::vector<Particle> send_L_buffer, send_R_buffer;  // statement
        #pragma omp parallel  // preprocessor directive
        {  // begin block
            std::vector<Particle> t_keep, t_send_L, t_send_R;  // statement
            #pragma omp for schedule(static)  // preprocessor directive
            for (size_t i = 0; i < local_particles.size(); ++i) {  // begin block
                double x_norm = (local_particles[i].x - x_local_min) / dx;  // function declaration/call or variable initialization
                int j = static_cast<int>(std::floor(x_norm));   // function declaration/call or variable initialization
// (blank line)
                // SAFETY CRITICAL CLAMP: Restricts lookups to prevent thread segmentation violations
                j = std::max(0, std::min(j, Nx_local - 1));  // variable assignment/initialization
                double xi = x_norm - j;  // variable assignment/initialization
                double E = (1.0 - xi) * Ex[j + 1] + xi * Ex[j + 2];  // function declaration/call or variable initialization
// (blank line)
                double g = std::sqrt(1.0 + (local_particles[i].ux*local_particles[i].ux + local_particles[i].uy*local_particles[i].uy + local_particles[i].uz*local_particles[i].uz)/(c*c));  // function declaration/call or variable initialization
                local_particles[i].ux -= (qe * dt / me) * E;  // variable assignment/initialization
                local_particles[i].x  += (local_particles[i].ux / g) * dt;  // variable assignment/initialization
// (blank line)
                if (local_particles[i].x < x_local_min)       t_send_L.push_back(local_particles[i]);  // control flow statement
                else if (local_particles[i].x >= x_local_max) t_send_R.push_back(local_particles[i]);  // else branch
                else                                         t_keep.push_back(local_particles[i]);  // else branch
            }  // end block
            #pragma omp critical  // preprocessor directive
            {  // begin block
                keep_particles.insert(keep_particles.end(), t_keep.begin(), t_keep.end());  // statement
                send_L_buffer.insert(send_L_buffer.end(), t_send_L.begin(), t_send_L.end());  // statement
                send_R_buffer.insert(send_R_buffer.end(), t_send_R.begin(), t_send_R.end());  // statement
            }  // end block
        }  // end block
// (blank line)
        // --- STEP E: Symmetric MPI Particle Migration ---
        int send_size_L = send_L_buffer.size(), recv_size_R = 0;  // function declaration/call or variable initialization
        int send_size_R = send_R_buffer.size(), recv_size_L = 0;  // function declaration/call or variable initialization
        MPI_Sendrecv(&send_size_L, 1, MPI_INT, left_neighbor, 1, &recv_size_R, 1, MPI_INT, right_neighbor, 1, MPI_COMM_WORLD, &status);  // statement
        MPI_Sendrecv(&send_size_R, 1, MPI_INT, right_neighbor, 2, &recv_size_L, 1, MPI_INT, left_neighbor, 2, MPI_COMM_WORLD, &status);  // statement
// (blank line)
        std::vector<Particle> recv_L_buffer(recv_size_L), recv_R_buffer(recv_size_R);  // function declaration/call or variable initialization
        if (send_size_L > 0 || recv_size_R > 0) MPI_Sendrecv(send_L_buffer.data(), send_size_L*sizeof(Particle), MPI_BYTE, left_neighbor, 3, recv_R_buffer.data(), recv_size_R*sizeof(Particle), MPI_BYTE, right_neighbor, 3, MPI_COMM_WORLD, &status);  // control flow statement
        if (send_size_R > 0 || recv_size_L > 0) MPI_Sendrecv(send_R_buffer.data(), send_size_R*sizeof(Particle), MPI_BYTE, right_neighbor, 4, recv_L_buffer.data(), recv_size_L*sizeof(Particle), MPI_BYTE, left_neighbor, 4, MPI_COMM_WORLD, &status);  // control flow statement
// (blank line)
        for (auto& p : recv_R_buffer) { if (p.x >= L_global) p.x -= L_global; keep_particles.push_back(p); }  // control flow statement
        for (auto& p : recv_L_buffer) { if (p.x < 0.0) p.x += L_global; keep_particles.push_back(p); }  // control flow statement
        local_particles = std::move(keep_particles);  // variable assignment/initialization
// (blank line)
        // --- STEP F: In-Situ Microscopic Parallel Heat Flux Compilation ---
        if (step % 20 == 0) {  // begin block
            std::fill(local_q_parallel.begin(), local_q_parallel.end(), 0.0);  // function declaration/call or variable initialization
            double* q_ptr = local_q_parallel.data();  // function declaration/call or variable initialization
            #pragma omp parallel for reduction(+:q_ptr[:Nx_local+1]) schedule(static)  // preprocessor directive
            for (size_t i = 0; i < local_particles.size(); ++i) {  // begin block
                if (local_particles[i].species == 0) {  // begin block
                    double x_norm = (local_particles[i].x - x_local_min) / dx;  // function declaration/call or variable initialization
                    int j = static_cast<int>(std::floor(x_norm));  // function declaration/call or variable initialization
                    j = std::max(0, std::min(j, Nx_local - 1));  // variable assignment/initialization
// (blank line)
                    double g = std::sqrt(1.0 + (local_particles[i].ux*local_particles[i].ux)/(c*c));  // function declaration/call or variable initialization
                    double vx = local_particles[i].ux / g;  // variable assignment/initialization
                    q_ptr[j] += 0.5 * me * (vx * vx * vx) * local_particles[i].w / dx;  // variable assignment/initialization
                }  // end block
            }  // end block
// (blank line)
            double local_q_sum = 0.0, global_q_avg = 0.0, local_max_E = 0.0, global_max_E = 0.0;  // variable assignment/initialization
            for(int j = 1; j <= Nx_local; ++j) { local_q_sum += std::abs(local_q_parallel[j-1]); local_max_E = std::max(local_max_E, std::abs(Ex[j])); }  // control flow statement
            MPI_Reduce(&local_q_sum, &global_q_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  // statement
            MPI_Reduce(&local_max_E, &global_max_E, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);  // statement
// (blank line)
            if(mpi_rank == 0) {  // begin block
                diag_out << step << "," << step*dt << "," << global_max_E << "," << (global_q_avg / Nx_global) << "\n";  // statement
            }  // end block
        }  // end block
    }  // end block
    if(mpi_rank == 0) {  // begin block
        diag_out.close();  // statement
        std::cout << "\n====================================================\n";  // variable assignment/initialization
        std::cout << " GOL-3 NON-LINEAR SCALING RUN COMPLETION            \n";  // statement
        std::cout << "====================================================\n";  // variable assignment/initialization
        std::cout << " Status: Complete. Telemetry saved to gol3_anomalous_conductivity.csv\n";  // statement
        std::cout << "====================================================\n";  // variable assignment/initialization
    }  // end block
    MPI_Finalize(); return 0;  // statement
}  // end block
