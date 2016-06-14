
# coding: utf-8

# In[ ]:
import numpy as np
import sympy
import fipy as fp
import matplotlib.pyplot as plt
import os
import sys

c, rho_s, c_alpha, c_beta = sympy.symbols("c_var rho_s c_alpha c_beta")
f_0 = rho_s * (c - c_alpha)**2 * (c_beta - c)**2

sympy.diff(f_0, c, 2)
nx = int(sys.argv[1])
dx = float(sys.argv[2])

mesh = fp.PeriodicGrid2D(nx=nx, ny=nx, dx=dx, dy=dx)

c_alpha = 0.3
c_beta = 0.7
kappa = 2.0
M = 5.0
c_0 = 0.5
epsilon = 0.01
rho_s = 5.0

c_var = fp.CellVariable(mesh=mesh, name=r"$c$", hasOld=True)

# array of sample c-values: used in f versus c plot
vals = np.linspace(-.1, 1.1, 1000)

c_var = fp.CellVariable(mesh=mesh, name=r"$c$", hasOld=True)

x , y = np.array(mesh.x), np.array(mesh.y)

c_var[:] = c_0 + epsilon * (np.cos(0.105 * x) * np.cos(0.11 * y) +                             (np.cos(0.13 * x) * np.cos(0.087 * y))**2 +                             + np.cos(0.025 * x - 0.15 * y) * np.cos(0.07 * x - 0.02 * y))

out = sympy.diff(f_0, c, 2)
exec "f_0_var = " + repr(out)
#f_0_var = -A + 3*B*(c_var - c_m)**2 + 3*c_alpha*(c_var - c_alpha)**2 + 3*c_beta*(c_var - c_beta)**2
f_0_var

# bulk free energy density
def f_0(c):
    return rho_s*((c - c_alpha)**2)*((c_beta-c)**2)
def f_0_var(c_var):
    return 2*rho_s*((c_alpha - c_var)**2 + 4*(c_alpha - c_var)*(c_beta - c_var) + (c_beta - c_var)**2)
# free energy
def f(c):
    return (f_0(c)+ .5*kappa*(c.grad.mag)**2)

# plot free energy density versus c
def plotf_c():
    plt.figure(1)
    plt.xlabel('c')
    plt.ylabel('f_0')
    plt.plot(vals, f_0(vals))
    plt.show()
    
# save elapsed time and free energy at each data point
f_data = []
time_data = []
file_name = "1a{0}x{1}.txt".format(nx, dx)

def save_data(f, time):
    f_data.append(f.value)
    time_data.append(time)
    np.savetxt(file_name, zip(time_data, f_data))
    
eqn = fp.TransientTerm(coeff=1.) == fp.DiffusionTerm(M * f_0_var(c_var)) - fp.DiffusionTerm((M, kappa))
eqn

elapsed = 0.0
steps = 0
dt = 0.01
total_sweeps = 2
tolerance = 1e-1
total_steps = 4000
duration = 3000.0
viewers = []

c_var[:] = c_0 + epsilon * (np.cos(0.105 * x) * np.cos(0.11 * y) +                             (np.cos(0.13 * x) * np.cos(0.087 * y))**2 +                             + np.cos(0.025 * x - 0.15 * y) * np.cos(0.07 * x - 0.02 * y))
c_var.updateOld()
from fipy.solvers.pysparse import LinearLUSolver as Solver
solver = Solver()
print "Starting Solver."
while steps <= total_steps:
    res0 = eqn.sweep(c_var, dt=dt, solver=solver)

    for sweeps in range(total_sweeps):
        res = eqn.sweep(c_var, dt=dt, solver=solver)

    if res < res0 * tolerance:
        # record the volume integral of the free energy 
        save_data(f(c_var).cellVolumeAverage*mesh.numberOfCells, elapsed)
        # pickle the data on c as a function of space at this particular time
        file_dir = "1apickles{0}x{1}".format(nx, dx)
        if not os.path.isdir(file_dir):
            os.makedirs(file_dir)
        fp.dump.write({'time' : steps, 'var': c_var}, '1apickles{0}x{1}/1a_{2}.pkl'.format(nx, dx, steps))
        
        if (steps%(total_steps/10)==0):
            print steps
            print elapsed
            
        steps += 1
        elapsed += dt
        dt *= 1.1
        c_var.updateOld()
    else:
        dt *= 0.8
        c_var[:] = c_var.old

viewer.plot()
print 'elapsed_time:', elapsed

