"""
Dedalus script for 3D Boussinesq Rayleigh-Benard convection within a rotating frame.

This script uses a Fourier basis in both the x and y directions with periodic boundary
conditions and a Chebyshev basis in the z direction with fixed flux boundary conditions.

Equations are non-dimensionalised by the viscous timescale (t in units of viscous time).
"""

import numpy as np
from mpi4py import MPI
import time

from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.tools import post
import pathlib

import logging
logger = logging.getLogger(__name__)

import run_param_file as rpf   # Imports a parameter file "run_param_file.py"

save_direc = "sim_data/3D/"
pathlib.Path(save_direc).mkdir(parents=True, exist_ok=True)


# Model Parameters
Lx, Ly, Lz = rpf.Lx, rpf.Ly, rpf.Lz
Nx, Ny, Nz = rpf.Nx, rpf.Ny, rpf.Nz
Pr = rpf.Pr
Ra = rpf.Ra
Ek = rpf.Ek
phi = rpf.phi
tau = rpf.tau

# Create bases and domain
x_basis = de.Fourier('x',Nx,interval=(-Lx/2.,Lx/2.),dealias=2)   # Fourier basis in x
y_basis = de.Fourier('y',Nx,interval=(-Ly/2.,Ly/2.),dealias=2)   # Fourier basis in y
z_basis = de.Chebyshev('z',Nz,interval=(0,Lz),dealias=2) # Chebyshev basis in z
domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64,mesh=[4,4])  # Defining domain
z = domain.grid(1, scales=1)                                   # accessing the z values

# Extra fields for time averaging and mean flow
average_fields = ['u_avgt', 'v_avgt', 'w_avgt', 'stress_uw', 'stress_vw']

# 3D Boussinesq hydrodynamics
problem = de.IVP(domain,variables=['T','p','u','v','w','Tz','uz','vz','wz', *average_fields])

# Defining model parameters
problem.parameters['Lx'] = Lx
problem.parameters['Ly'] = Ly
problem.parameters['Lz'] = Lz
problem.parameters['Ra'] = Ra
problem.parameters['Pr'] = Pr
problem.parameters['Ek'] = Ek
problem.parameters['phi'] = phi
problem.parameters['tau'] = tau 
problem.parameters['X'] = Ra/Pr


# Defining d/dz of T, u, v and w for reducing our equations to first order
problem.add_equation("dz(u) - uz = 0")
problem.add_equation("dz(v) - vz = 0")
problem.add_equation("dz(w) - wz = 0")
problem.add_equation("dz(T) - Tz = 0")

# mass continuity
problem.add_equation("dx(u) + dy(v) + wz = 0")
# x-component of the momentum equation
problem.add_equation("dt(u) + dx(p) - (dx(dx(u)) + dy(dy(u)) + dz(uz)) - 2 * (1 / Ek) * v * sin(phi) = - (u * dx(u) + v * dy(u) + w * uz)")
# y-component of the momentum equation
problem.add_equation("dt(v) + dy(p) - (dx(dx(v)) + dy(dy(v)) + dz(vz)) + 2 * (1 / Ek) * (u * sin(phi) - w * cos(phi)) = - (u * dx(v) + v * dy(v) + w * vz)")
# z-component of the momentum equation
problem.add_equation("dt(w) + dz(p) - (dx(dx(w)) + dy(dy(w)) + dz(wz)) - 2 * (1 / Ek) * v * cos(phi) - X * T = -(u * dx(w) + v * dy(w) + w * wz)")
# Temperature equation
problem.add_equation("Pr * dt(T) - (dx(dx(T)) + dy(dy(T)) + dz(Tz)) = - Pr * (u * dx(T) + v * dy(T) + w * Tz)")

# Time averaging
problem.add_equation("dt(u_avgt) = u / tau") # time averaged u
problem.add_equation("dt(v_avgt) = v / tau") # time averaged v
problem.add_equation("dt(w_avgt) = w / tau") # time averaged w
problem.add_equation("dt(stress_uw) = ((u - u_avgt) * (w - w_avgt)) / tau") # time averaged reynolds stress
problem.add_equation("dt(stress_vw) = ((v - v_avgt) * (w - w_avgt)) / tau") # time averaged reynolds stress

problem.add_bc("left(u) = 0")           # no-slip boundary
problem.add_bc("right(u) = 0")          # no-slip boundary
problem.add_bc("left(v) = 0")           # no-slip boundary
problem.add_bc("right(v) = 0")          # no-slip boundary
problem.add_bc("left(w) = 0")            # Impermeable bottom boundary
problem.add_bc("right(w) = 0",condition="(nx != 0) or (ny != 0)")   # Impermeable top boundary
problem.add_bc("right(p) = 0",condition="(nx == 0) and (ny == 0)")   # Required for equations to be well-posed - see https://bit.ly/2nPVWIg for a related discussion
problem.add_bc("right(T) = 0")           # Fixed temperature at upper boundary
problem.add_bc("left(Tz) = -1")           # Fixed flux at bottom boundary, F = F_cond


# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
logger.info('Solver built')

if not pathlib.Path('restart.h5').exists():
    # Initial conditions
    z = domain.grid(1)
    T = solver.state['T']
    Tz = solver.state['Tz']

    # Random perturbations, initialized globally for same results in parallel
    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand = np.random.RandomState(seed=42)
    noise = rand.standard_normal(gshape)[slices]

    # Linear background + perturbations damped at walls
    zb, zt = z_basis.interval
    pert =  1e-5 * noise * (zt - z) * (z - zb)
    T['g'] = 1-pert
    T.differentiate('z', out=Tz)

    # Initial timestep
    dt = rpf.initial_timestep

    fh_mode = 'overwrite'
else:
    #Restart
    write, last_dt = solver.load_state('restart.h5', -1)

    dt = last_dt
    stop_sim_time

# Integration parameters --- Note if these are all set to np.inf, simulation will perpetually run.
solver.stop_sim_time = rpf.end_sim_time
solver.stop_wall_time = rpf.end_wall_time
solver.stop_iteration = rpf.end_iterations

# CFL criterion
CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=0.5,max_change=1.5, min_change=0.5, max_dt=rpf.max_dt, threshold=0.05)
CFL.add_velocities(('u', 'v', 'w'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("sqrt(u**2 + v**2 + w**2)/Ra", name='Re')

# Saving snapshots
snapshots = solver.evaluator.add_file_handler(save_direc + 'snapshots', sim_dt=rpf.snapshot_freq, mode=fh_mode)
snapshots.add_system(solver.state)

# Analysis tasks
analysis = solver.evaluator.add_file_handler(save_direc + 'analysis', sim_dt=rpf.analysis_freq, mode=fh_mode)

# Mean Re
analysis.add_task("integ( integ( integ( sqrt(u**2 + v**2 + w**2) , 'x')/Lx , 'y')/Ly , 'z')/Lz", layout='g', name='Re')

# Flux decomposition - Internal energy equation
analysis.add_task("integ( integ( T*w,'x')*Pr/Lx , 'y')/Ly", layout='g', name='L_conv')
analysis.add_task("integ( integ( (-1)*Tz, 'x')/Lx , 'y')/Ly", layout='g', name='L_cond')

# Mean KE
analysis.add_task("integ( integ( integ(0.5*(u**2 + v**2 + w**2) , 'x')/Lx , 'y')/Ly , 'z')/Lz", layout='g', name='KE')

# Averages & Reynolds Stresses
# Must be called manually so that averages can be reset at each interval
averages = solver.evaluator.add_file_handler(save_direc + 'averages', sim_dt = np.inf, mode='overwrite')

# Warning: dirty hack
# Make dedalus think this handler has been evaluated at time t=0
# This stops the handlers being evaluated at the start:
averages.last_wall_div = averages.last_sim_div = averages.last_iter_div = 0

# Averages to save
averages.add_task("integ(u_avgt, 'y') / Ly", layout = 'g', name = 'u_avgt')
averages.add_task("integ(v_avgt, 'y') / Ly", layout = 'g', name = 'v_avgt')
averages.add_task("integ(w_avgt, 'y') / Ly", layout = 'g', name = 'w_avgt')
averages.add_task("integ(stress_uw, 'y') / Ly", layout='g', name='stress_uw')
averages.add_task("integ(stress_vw, 'y') / Ly", layout='g', name='stress_vw')
averages.add_task("dz(integ(stress_uw, 'y')) / Ly", layout='g', name='dz_stress_uw')
averages.add_task("dz(integ(stress_vw, 'y')) / Ly", layout='g', name='dz_stress_vw')

# Creating a parameter file
run_parameters = solver.evaluator.add_file_handler(save_direc + 'run_parameters', wall_dt=1e20, max_writes=1,mode=fh_mode)
run_parameters.add_task(Lx, name="Lx")
run_parameters.add_task(Ly, name="Ly")
run_parameters.add_task(Lz, name="Lz")
run_parameters.add_task(Ra, name="Ra")
run_parameters.add_task(Pr, name="Pr")
run_parameters.add_task(Ek, name="Ek")
run_parameters.add_task(phi, name="phi")
run_parameters.add_task(Nx, name="Nx")
run_parameters.add_task(Ny, name="Ny")
run_parameters.add_task(Nz, name="Nz")

run_parameters.add_task(rpf.snapshot_freq, name="snap_freq")
run_parameters.add_task(rpf.analysis_freq, name="ana_freq")
run_parameters.add_task(rpf.max_dt, name="max_dt")
run_parameters.add_task(tau, name = 'tau')

# Function in order to reset averages
def reset_averages():
    for var in average_fields:
        field = solver.state[var]
        nx, ny, nz = np.array(field['g']).shape
        field['g'] = np.zeros((nx,ny,nz))
    print('Reset time-averaged fields at time t = {}'.format(solver.sim_time))

# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()

    next_average = (np.ceil(solver.sim_time / tau) + 1) * tau

    while solver.ok:
        dt = CFL.compute_dt()
        dt = solver.step(dt)

        if (solver.iteration) == 1:
            # Prints various parameters to terminal upon starting the simulation
            logger.info('Parameter values imported form run_param_file.py:')
            logger.info('Lx = {}, Ly = {}, Lz = {}; (Resolution of {},{},{})'.format(Lx, Ly, Lz, Nx, Ny, Nz))
            logger.info('Ra = {}, Pr = {}, Ek = {}, phi = {}'.format(Ra, Pr, Ek, phi))
            logger.info('Snapshot files outputted every {}'.format(rpf.snapshot_freq))
            logger.info('Analysis files outputted every {}'.format(rpf.analysis_freq))
            if rpf.end_sim_time != np.inf:
                logger.info('Simulation finishes at sim_time = {}'.format(rpf.end_sim_time))
            elif rpf.end_wall_time != np.inf:
                logger.info('Simulation finishes at wall_time = {}'.format(rpf.end_wall_time))
            elif rpf.end_iterations != np.inf:
                logger.info('Simulation finishes at iteration {}'.format(rpf.end_iterations))
            else:
                logger.info('No clear end point defined. Simulation may run perpetually.')

        # When we reach the end of a time-averaging interval, save to analysis file
        if solver.sim_time >= next_average:
            next_average += tau
            solver.evaluate_handlers_now(dt, handlers=[averages])
            print(''.join([' '] * 200), end='\r')
            print('Calculated time-averaged fields at time t={}'.format(solver.sim_time))

        # When we reach the start of a time-averaging interval, zero the averages to restart them
        if solver.sim_time % tau < dt:
            print(''.join([' '] * 200), end='\r')
            reset_averages()

        if (solver.iteration-1) % 100 == 0:
            # Prints progress information include maximum Reynolds number every 100 iterations
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('Max Re = %f' %flow.max('Re'))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    # Prints concluding information upon reaching the end of the simulation.
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))
