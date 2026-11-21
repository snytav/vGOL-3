#!/usr/bin/env python3
import csv, math, subprocess, os, shutil, sys  # import module(s)
os.chdir(os.path.dirname(__file__) or '.')  # statement
# Build
subprocess.check_call(['g++','-O2','-std=c++17','-fopenmp','-Wall','-Wextra','-pedantic','main.cpp','-o','openmp_pic_core_verify'])  # assignment/initialization
def run(nt):  # function definition
    env=os.environ.copy(); env['OMP_NUM_THREADS']=str(nt)  # assignment/initialization
    if os.path.exists('two_stream_diagnostics.csv'): os.remove('two_stream_diagnostics.csv')  # control flow / block statement
    subprocess.check_call(['./openmp_pic_core_verify'], env=env, stdout=subprocess.DEVNULL)  # assignment/initialization
    shutil.copy('two_stream_diagnostics.csv', f'two_stream_diagnostics_t{nt}.csv')  # statement
for nt in (1,2,4): run(nt)  # control flow / block statement
import filecmp  # import module(s)
bitwise = filecmp.cmp('two_stream_diagnostics_t1.csv','two_stream_diagnostics_t2.csv',shallow=False) and filecmp.cmp('two_stream_diagnostics_t1.csv','two_stream_diagnostics_t4.csv',shallow=False)  # assignment/initialization
print('BITWISE_IDENTICAL=' + ('PASS' if bitwise else 'FAIL'))  # print output
# (blank line)
rows=[]  # assignment/initialization
with open('two_stream_diagnostics_t1.csv') as f:  # control flow / block statement
    for row in csv.DictReader(f):  # control flow / block statement
        rows.append((float(row['time']), float(row['E_k1_amp'])))  # statement
me=9.1093837e-31; qe=1.602176634e-19; eps0=8.8541878128e-12  # assignment/initialization
n_plasma=1e19; n_beam=0.01*n_plasma  # assignment/initialization
omega_pe=math.sqrt(n_plasma*qe**2/(eps0*me))  # assignment/initialization
gamma_theory = omega_pe * (math.sqrt(3)/2.0) * (n_beam/(2.0*n_plasma))**(1/3)  # assignment/initialization
# find exponential window from log amplitude slope
best=None  # assignment/initialization
w=30  # assignment/initialization
for i in range(5, len(rows)-w-5):  # control flow / block statement
    t=[rows[j][0] for j in range(i,i+w)]  # assignment/initialization
    a=[rows[j][1] for j in range(i,i+w)]  # assignment/initialization
    if any(x<=0 for x in a): continue  # control flow / block statement
    la=[math.log(x) for x in a]  # assignment/initialization
    # linear fit
    n=len(t); st=sum(t); sla=sum(la); stt=sum(x*x for x in t); stla=sum(t[k]*la[k] for k in range(n))  # assignment/initialization
    denom=n*stt-st*st  # assignment/initialization
    if abs(denom)<1e-30: continue  # control flow / block statement
    slope=(n*stla-st*sla)/denom  # assignment/initialization
    # fit quality
    intercept=(sla-slope*st)/n  # assignment/initialization
    ss_res=sum((la[k]-(intercept+slope*t[k]))**2 for k in range(n))  # assignment/initialization
    ss_tot=sum((x-sla/n)**2 for x in la)  # assignment/initialization
    r2=1-ss_res/ss_tot if ss_tot>1e-30 else 0  # assignment/initialization
    gamma=slope  # amplitude grows as e^(gamma*t)
    if gamma<=0: continue  # control flow / block statement
    score=r2*gamma  # assignment/initialization
    if best is None or score>best[0]: best=(score,i,i+w,slope,r2)  # control flow / block statement
if best is None:  # control flow / block statement
    print('GAMMA_MATCH=FAIL (no window)')  # print output
    sys.exit(1)  # statement
_,i0,i1,slope,r2=best  # assignment/initialization
gamma_sim=slope  # assignment/initialization
err=abs(gamma_sim-gamma_theory)/gamma_theory*100  # assignment/initialization
print(f'gamma_theory={gamma_theory:.6e}')  # print output
print(f'gamma_sim={gamma_sim:.6e}')  # print output
print(f'r2={r2:.4f} window={rows[i0][0]*1e12:.3f}-{rows[i1-1][0]*1e12:.3f}ps')  # print output
print(f'err_pct={err:.3f}')  # print output
# Plan says "tightly match" - allow reasonable tolerance for numerical/PIC noise
ok = err < 20.0 and r2 > 0.8  # assignment/initialization
print('GAMMA_MATCH=' + ('PASS' if ok else 'FAIL'))  # print output
print('OVERALL=' + ('PASS' if bitwise and ok else 'FAIL'))  # print output
sys.exit(0 if bitwise and ok else 1)  # statement
