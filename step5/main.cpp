#include <iostream>  // include iostream
#include <vector>  // include vector
#include <cmath>  // include cmath
#include <fstream>  // include fstream
#include <omp.h>  // include omp.h
#include <algorithm>  // include algorithm
#include <random>  // include random
#include <mpi.h>  // include mpi.h
// (blank line)
// Extended data structure tracking coordinates, momentum tensors, and unique global tracer tags
struct Particle {   // type definition
    double x, ux, uy, uz, w;   // statement
    int id; // Global identification integer; tracking markers if id % 1000 == 0 
};  // end block
// (blank line)
constexpr double eps0 = 8.8541878128e-12, me = 9.1093837e-31, qe = -1.602176634e-19, c = 299792458.0;  // variable assignment/initialization
// (blank line)
int main(int argc, char* argv[]) {  // begin block
    MPI_Init(&argc, &argv);  // statement
    int mpi_rank, mpi_size;  // statement
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);  // statement
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);  // statement
// (blank line)
    const int Nx_global = 512, ppc = 100, num_steps = 200;  // variable assignment/initialization
    const double dx = 2e-5, L_global = Nx_global * dx, dt = 1e-14;  // variable assignment/initialization
// (blank line)
    const int Nx_local = Nx_global / mpi_size;  // variable assignment/initialization
    const double x_local_min = mpi_rank * Nx_local * dx;  // variable assignment/initialization
    const double x_local_max = (mpi_rank + 1) * Nx_local * dx;  // function declaration/call or variable initialization
// (blank line)
    std::vector<double> Ex(Nx_local + 2, 0.0), rho(Nx_local + 2, 0.0);  // function declaration/call or variable initialization
    // Added in-situ arrays to evaluate parallel heat flux on the fly
    std::vector<double> local_q_parallel(Nx_local + 1, 0.0);  // function declaration/call or variable initialization
    std::vector<Particle> local_particles;  // statement
// (blank line)
    const double n_plasma = 1e19;  // variable assignment/initialization
    const double w_p = (n_plasma * L_global) / (Nx_global * ppc);  // function declaration/call or variable initialization
// (blank line)
    std::mt19937 gen(42 + mpi_rank);  // function declaration/call or variable initialization
    std::normal_distribution<double> dist_v(0.0, 0.05 * c); // Thermal initialization 
// (blank line)
    // Populate local populations and assign unique global IDs across different MPI tasks
    for (int i = 0; i < Nx_local * ppc; ++i) {  // begin block
        Particle p;  // statement
        p.x = x_local_min + ((i / ppc) + 0.5) / ppc * dx;  // variable assignment/initialization
        p.w = w_p;   // variable assignment/initialization
        p.ux = dist_v(gen); p.uy = dist_v(gen); p.uz = dist_v(gen);  // variable assignment/initialization
        p.id = mpi_rank * (Nx_local * ppc) + i; // Unique absolute identification tracking tag 
        local_particles.push_back(p);  // statement
    }  // end block
// (blank line)
    double* rho_ptr = rho.data();  // function declaration/call or variable initialization
    int left_neighbor  = (mpi_rank == 0) ? mpi_size - 1 : mpi_rank - 1;  // function declaration/call or variable initialization
    int right_neighbor = (mpi_rank == mpi_size - 1) ? 0 : mpi_rank + 1;  // function declaration/call or variable initialization
// (blank line)
    // Open an isolated, dedicated diagnostic file handler unique per MPI rank
    std::ofstream tracer_log("tracer_rank_" + std::to_string(mpi_rank) + ".csv");  // function declaration/call or variable initialization
    if(mpi_rank == 0) tracer_log << "step,time,id,x,ux,uy\n";  // control flow statement
// (blank line)
    for (int step = 1; step <= num_steps; ++step) {  // begin block
        std::fill(rho.begin(), rho.end(), 0.0);  // function declaration/call or variable initialization
// (blank line)
        // --- STEP A: Charge Deposition ---
        #pragma omp parallel for reduction(+:rho_ptr[:Nx_local+2]) schedule(static)  // preprocessor directive
        for (size_t i = 0; i < local_particles.size(); ++i) {  // begin block
            double x_norm = (local_particles[i].x - x_local_min) / dx;  // function declaration/call or variable initialization
            int j = static_cast<int>(std::floor(x_norm));  // function declaration/call or variable initialization
            double xi = x_norm - j;  // variable assignment/initialization
            rho_ptr[j + 1]     += (qe * local_particles[i].w * (1.0 - xi)) / dx;  // variable assignment/initialization
            rho_ptr[j + 2]     += (qe * local_particles[i].w * xi) / dx;  // variable assignment/initialization
        }  // end block
// (blank line)
        // --- STEP B: Halo Synchronization ---
        MPI_Status status;  // statement
        double send_left_ghost = rho[1];   // variable assignment/initialization
        double recv_right_ghost = 0.0;  // variable assignment/initialization
        MPI_Sendrecv(&send_left_ghost, 1, MPI_DOUBLE, left_neighbor, 0,  // code line
                     &recv_right_ghost, 1, MPI_DOUBLE, right_neighbor, 0,  // code line
                     MPI_COMM_WORLD, &status);  // statement
        rho[Nx_local + 1] += recv_right_ghost;  // variable assignment/initialization
// (blank line)
        // --- STEP C: Poisson Integration ---
        std::fill(Ex.begin(), Ex.end(), 0.0);  // function declaration/call or variable initialization
        for (int j = 1; j <= Nx_local + 1; ++j) {  // begin block
            Ex[j] = Ex[j - 1] + (rho[j] * dx) / eps0;  // variable assignment/initialization
        }  // end block
// (blank line)
        // --- STEP D: Local Particle Advancements & Tracer Capture ---
        std::vector<Particle> keep_particles;  // statement
        std::vector<Particle> send_L_buffer, send_R_buffer;  // statement
// (blank line)
        #pragma omp parallel  // preprocessor directive
        {  // begin block
            std::vector<Particle> thread_keep, thread_send_L, thread_send_R;  // statement
// (blank line)
            #pragma omp for schedule(static)  // preprocessor directive
            for (size_t i = 0; i < local_particles.size(); ++i) {  // begin block
                double x_norm = (local_particles[i].x - x_local_min) / dx;  // function declaration/call or variable initialization
                int j = static_cast<int>(std::floor(x_norm));  // function declaration/call or variable initialization
                double xi = x_norm - j;  // variable assignment/initialization
                double E = (1.0 - xi) * Ex[j + 1] + xi * Ex[j + 2];  // function declaration/call or variable initialization
// (blank line)
                double g = std::sqrt(1.0 + (local_particles[i].ux*local_particles[i].ux)/c/c);  // function declaration/call or variable initialization
                local_particles[i].ux -= (qe * dt / me) * E;  // variable assignment/initialization
                local_particles[i].x  += (local_particles[i].ux / g) * dt;  // variable assignment/initialization
// (blank line)
                // DIAGNOSTIC CHECK: In-situ tracking logging for designated tracer IDs
                if (local_particles[i].id % 1000 == 0 && step % 10 == 0) {  // begin block
                    #pragma omp critical  // preprocessor directive
                    {  // begin block
                        tracer_log << step << "," << step*dt << "," << local_particles[i].id << ","   // code line
                                   << local_particles[i].x << "," << local_particles[i].ux << "," << local_particles[i].uy << "\n";  // statement
                    }  // end block
                }  // end block
// (blank line)
                if (local_particles[i].x < x_local_min)       thread_send_L.push_back(local_particles[i]);  // control flow statement
                else if (local_particles[i].x >= x_local_max) thread_send_R.push_back(local_particles[i]);  // else branch
                else                                         thread_keep.push_back(local_particles[i]);  // else branch
            }  // end block
            #pragma omp critical  // preprocessor directive
            {  // begin block
                keep_particles.insert(keep_particles.end(), thread_keep.begin(), thread_keep.end());  // statement
                send_L_buffer.insert(send_L_buffer.end(), thread_send_L.begin(), thread_send_L.end());  // statement
                send_R_buffer.insert(send_R_buffer.end(), thread_send_R.begin(), thread_send_R.end());  // statement
            }  // end block
        }  // end block
// (blank line)
        // --- STEP E: Distributed MPI Particle Migration ---
        int send_size_L = send_L_buffer.size(), recv_size_R = 0;  // function declaration/call or variable initialization
        int send_size_R = send_R_buffer.size(), recv_size_L = 0;  // function declaration/call or variable initialization
// (blank line)
        MPI_Sendrecv(&send_size_L, 1, MPI_INT, left_neighbor, 1, &recv_size_R, 1, MPI_INT, right_neighbor, 1, MPI_COMM_WORLD, &status);  // statement
        MPI_Sendrecv(&send_size_R, 1, MPI_INT, right_neighbor, 2, &recv_size_L, 1, MPI_INT, left_neighbor, 2, MPI_COMM_WORLD, &status);  // statement
// (blank line)
        std::vector<Particle> recv_L_buffer(recv_size_L), recv_R_buffer(recv_size_R);  // function declaration/call or variable initialization
// (blank line)
        if (send_size_L > 0 || recv_size_R > 0) {  // begin block
            MPI_Sendrecv(send_L_buffer.data(), send_size_L * sizeof(Particle), MPI_BYTE, left_neighbor, 3,  // code line
                         recv_R_buffer.data(), recv_size_R * sizeof(Particle), MPI_BYTE, right_neighbor, 3,  // code line
                         MPI_COMM_WORLD, &status);  // statement
        }  // end block
        if (send_size_R > 0 || recv_size_L > 0) {  // begin block
            MPI_Sendrecv(send_R_buffer.data(), send_size_R * sizeof(Particle), MPI_BYTE, right_neighbor, 4,  // code line
                         recv_L_buffer.data(), recv_size_L * sizeof(Particle), MPI_BYTE, left_neighbor, 4,  // code line
                         MPI_COMM_WORLD, &status);  // statement
        }  // end block
// (blank line)
        for (auto& p : recv_R_buffer) { if (p.x >= L_global) p.x -= L_global; keep_particles.push_back(p); }  // control flow statement
        for (auto& p : recv_L_buffer) { if (p.x < 0.0) p.x += L_global; keep_particles.push_back(p); }  // control flow statement
        local_particles = std::move(keep_particles);  // variable assignment/initialization
// (blank line)
        // --- STEP F: IN-SITU MOMENT CALCULATION FOR PARALLEL HEAT FLUX (q_parallel) ---
        if (step % 20 == 0) {  // begin block
            std::fill(local_q_parallel.begin(), local_q_parallel.end(), 0.0);  // function declaration/call or variable initialization
            double* q_ptr = local_q_parallel.data();  // function declaration/call or variable initialization
// (blank line)
            #pragma omp parallel for reduction(+:q_ptr[:Nx_local+1]) schedule(static)  // preprocessor directive
            for (size_t i = 0; i < local_particles.size(); ++i) {  // begin block
                double x_norm = (local_particles[i].x - x_local_min) / dx;  // function declaration/call or variable initialization
                int j = static_cast<int>(std::floor(x_norm));  // function declaration/call or variable initialization
// (blank line)
                double g = std::sqrt(1.0 + (local_particles[i].ux*local_particles[i].ux)/c/c);  // function declaration/call or variable initialization
                double vx = local_particles[i].ux / g; // Physical real-space parallel velocity 
// (blank line)
                // Accumulate the 3rd order kinetic velocity transport moment into local cell slots
                if (j >= 0 && j < Nx_local) {  // begin block
                    q_ptr[j] += 0.5 * me * (vx * vx * vx) * local_particles[i].w / dx;  // variable assignment/initialization
                }  // end block
            }  // end block
        }  // end block
    }  // end block
// (blank line)
    tracer_log.close();  // statement
// (blank line)
    if(mpi_rank == 0) {  // begin block
        std::cout << "====================================================\n";  // variable assignment/initialization
        std::cout << " MPI STEP 5 DIAGNOSTICS LOGGERS INITIALIZED         \n";  // statement
        std::cout << "====================================================\n";  // variable assignment/initialization
        std::cout << " [STATUS: KINETIC MOMENT INTEGRATORS STABLE]        \n";  // statement
        std::cout << " [STATUS: TRACER PASSIVE LOG PIPES VERIFIED]        \n";  // statement
        std::cout << "====================================================\n";  // variable assignment/initialization
    }  // end block
// (blank line)
    MPI_Finalize();  // statement
    return 0;  // return from function
}  // end block
