#!/usr/bin/env python3
import csv, math, os, subprocess, sys, shutil  # import module(s)
os.chdir(os.path.dirname(__file__) or '.')  # statement
# Build and run with 4 ranks
subprocess.check_call(['mpicxx','-O2','-std=c++17','-fopenmp','-Wall','-Wextra','-pedantic','main.cpp','-o','hybrid_mpi_core_verify'])  # assignment/initialization
env=os.environ.copy(); env['OMP_NUM_THREADS']='1'  # assignment/initialization
r=subprocess.run(['mpirun','--oversubscribe','-np','4','./hybrid_mpi_core_verify'], env=env, capture_output=True, text=True)  # assignment/initialization
print(r.stdout)  # print output
print(r.stderr, file=sys.stderr)  # print output
sys.exit(r.returncode)  # statement
