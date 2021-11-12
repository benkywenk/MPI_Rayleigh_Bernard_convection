#Plotting tool for 2.5D simulations
from dedalus import public as de
import matplotlib.pyplot as plt
import matplotlib.animation as ani

import h5py
import numpy as np
import sys
import pathlib
import os

import run_param_file as rpf

direc = "sim_data/2halfD/"
save_direc = "figures/2halfD/"
run_name = rpf.run_name

plot_fluxes = False
plot_KE = False
mean_flow = True

if os.path.exists(save_direc) == False:
    pathlib.Path(save_direc).mkdir(parents=True)

# Importing run parameters
with h5py.File(direc + "run_parameters/run_parameters_"+ run_name +".h5", mode = 'r') as file:
    Pr = file['tasks']['Pr'][0][0][0]
    Ra = file['tasks']['Ra'][0][0][0]
    Ek = file['tasks']['Ek'][0][0][0]
    Lx = int(file['tasks']['Lx'][0][0][0])
    Lz = int(file['tasks']['Lz'][0][0][0])
    Nx = int(file['tasks']['Nx'][0][0][0])
    Nz = int(file['tasks']['Nz'][0][0][0])
    x = np.linspace(-Lx/2,Lx/2,Nx)

    z_basis = de.Chebyshev('z', Nz, interval=(0,Lz), dealias=3/2)
    z = np.array(z_basis.grid(1))

    print("Ra = {}".format(Ra))
    print("Pr = {}".format(Pr))
    print("Ek = {}".format(Ek))
    print("Resolution = ({},{})".format(Nx,Nz))

# Importing analysis tasks
with h5py.File(direc + "analysis/analysis_"+ run_name +".h5", mode = 'r') as file:
    L_cond = np.array(file['tasks']['L_cond'])[:,0,:]
    L_conv = np.array(file['tasks']['L_conv'])[:,0,:]
    L_tot = L_cond + L_conv

    KE = np.array(file['tasks']['KE'])[:,0,0]
    t = np.array(file['scales']['sim_time'])

    ubar = np.array(file['tasks']['ubar_x'])
    vbar = np.array(file['tasks']['vbar_x'])
    wbar = np.array(file['tasks']['wbar_x'])

# Importing system snapshots
with h5py.File(direc + "snapshots/snapshots_"+ run_name +".h5", mode = 'r') as file:
    #Empty

# Plot routine for energy flux
if plot_fluxes:
    plt.plot(L_cond[-1],z, 'r', linestyle='-', label="$L_{cond}$")
    plt.plot(L_conv[-1],z, 'g', linestyle='-', label="$L_{conv}$")
    plt.plot(L_tot[-1],z, 'k', linestyle='-',  label="$L_{total}$")
    plt.xlabel("L")
    plt.ylabel("z")
    plt.title("Ra = {}, Pr = {}, Ek = {}, phi = {}".format(Ra, Pr, Ek, phi))
    plt.legend()
    plt.savefig(save_direc + "intE_fluxes_"+ run_name)
    plt.clf()
    plt.close()

# Plot routine for system kinetic energy
if plot_KE:
    fig = plt.figure(figsize = (10,10))
    plt.plot(t,KE)
    plt.xlabel("t")
    plt.ylabel("Kinetic Energy")
    plt.title("Ra = {}, Pr = {}, Ek = {}, phi = {}".format(Ra, Pr, Ek, phi))
    plt.savefig(save_direc + "KE_"+ run_name)
    plt.clf
    plt.close

# Plot routine for mean flow plots
if mean_flow:
    fig, ax = plt.subplots(2,sharex = True)
    ax[0].contourf(t,z,ubar[:,0,:].T, cmap = 'jet', levels = 200)
    ax[1].contourf(t,z,vbar[:,0,:].T, cmap = 'jet', levels = 200)
    ax[0].set_title('Horizontally averaged u')
    ax[0].set_ylabel('z')
    ax[0].set_xlabel('t')
    ax[1].set_title('Horizontally averaged v')
    ax[1].set_ylabel('z')
    ax[1].set_xlabel('t')
    plt.savefig(save_direc+"mean_flow_"+run_name)
    plt.clf()
    plt.close()
