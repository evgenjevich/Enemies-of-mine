
import sympy
import fipy as fp
import numpy as np
import sys
import os
import resource
import json
from fipy.solvers.pysparse import LinearLUSolver as Solver
import profile

def run():
    memory1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss #get initial memory peak

    #load in the json parameter file here
    jsonfile = sys.argv[1]

    if jsonfile:
        with open(jsonfile, 'rb') as ff:
            params = json.load(ff)

    else:
        params = dict()

    print 'my params:', params

    #extract the parameters
    N = params.get('N', 20)
    # dx = params.get('dx', 1)  
    total_steps = params.get('steps', 2)
    sumatra_label = params.get('sumatra_label', '')

    c, rho_s, c_alpha, c_beta = sympy.symbols("c_var rho_s c_alpha c_beta")
    f_0 = rho_s * (c - c_alpha)**2 * (c_beta - c)**2

    mesh = fp.PeriodicGrid2D(nx=N, ny=N, dx=.5, dy=.5)

    c_alpha = 0.3
    c_beta = 0.7
    kappa = 2.0
    M = 5.0
    c_0 = 0.5
    epsilon = 0.01
    rho_s = 5.0
    filepath = os.path.join('Data', sumatra_label)

    c_var = fp.CellVariable(mesh=mesh, name=r"$c$", hasOld=True)

    # array of sample c-values: used in f versus c plot
    vals = np.linspace(-.1, 1.1, 1000)

    c_var = fp.CellVariable(mesh=mesh, name=r"$c$", hasOld=True)

    x , y = np.array(mesh.x), np.array(mesh.y)

    out = sympy.diff(f_0, c, 2)

    exec "f_0_var = " + repr(out)
    #f_0_var = -A + 3*B*(c_var - c_m)**2 + 3*c_alpha*(c_var - c_alpha)**2 + 3*c_beta*(c_var - c_beta)**2

    def f_0(c):
        return rho_s*((c - c_alpha)**2)*((c_beta-c)**2)
    def f_0_var(c_var):
        return 2*rho_s*((c_alpha - c_var)**2 + 4*(c_alpha - c_var)*(c_beta - c_var) + (c_beta - c_var)**2)
    # free energy
    def f(c):
        return (f_0(c)+ .5*kappa*(c.grad.mag)**2)
    f_data = []
    time_data = []

    def save_data(f, time):
        f_data.append(f.value)
        time_data.append(time)
        np.savetxt(os.path.join(filepath, '1a.txt'), zip(time_data, f_data))

    eqn = fp.TransientTerm(coeff=1.) == fp.DiffusionTerm(M * f_0_var(c_var)) - fp.DiffusionTerm((M, kappa))

    elapsed = 0.0
    steps = 0
    dt = 0.01
    total_sweeps = 2
    tolerance = 1e-1
    # duration = 1000.0

    c_var[:] = c_0 + epsilon * (np.cos(0.105 * x) * np.cos(0.11 * y) + \
                                (np.cos(0.13 * x) * np.cos(0.087 * y))**2 + \
                                + np.cos(0.025 * x - 0.15 * y) * np.cos(0.07 * x - 0.02 * y))
    c_var.updateOld()
    solver = Solver()

    while steps < total_steps:
        res0 = eqn.sweep(c_var, dt=dt, solver=solver)

        for sweeps in range(total_sweeps):
            res = eqn.sweep(c_var, dt=dt, solver=solver)

        if (res < (res0 * tolerance)):
            steps += 1
    #       elapsed += dt
            dt *= 1.1
            c_var.updateOld()

            if (steps%(total_steps/10.0)==0):
                # record the volume integral of the free energy 
                save_data(f_0_var(c_var).cellVolumeAverage*mesh.numberOfCells, elapsed)
                # pickle the data on c as a function of space at this particular time
                fp.dump.write({'time' : steps, 'var': c_var}, os.path.join(filepath,'1a{0}.pkl'.format(steps)))

        else:
            dt *= 0.8
            c_var[:] = c_var.old

    print ' '

    #memory stuff saves
    filepath = os.path.join('Data', sumatra_label)
    #Keep track os how much memory was used and dump into a txt file
    memory2 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss #final memory peak
    memory_diff = (memory2 - memory1,)
    filename2 = 'memory_usage.txt'
    np.savetxt(os.path.join(filepath, filename2), memory_diff )

profile.run('run(); print')
