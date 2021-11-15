#Plotting tool for 3D simulations
from dedalus import public as de
import matplotlib.pyplot as plt
import matplotlib.animation as ani

import h5py
import numpy as np
import sys
import pathlib
import os

import run_param_file as rpf

direc = "sim_data/3D/"
save_direc = "figures/3D/"
run_name = rpf.run_name

plot_fluxes = False
plot_KE = False
make_animation = False
velocity_snapshot = False
temp_snapshot = False
reynolds_stress = True
stress = True
mean_flow = True

if os.path.exists(save_direc) == False:
    pathlib.Path(save_direc).mkdir(parents=True)

with h5py.File(direc + "run_parameters/run_parameters_" + run_name + ".h5", mode = 'r') as file:
    Pr = file['tasks']['Pr'][0][0][0]
    Ra = file['tasks']['Ra'][0][0][0]
    Ek = file['tasks']['Ek'][0][0][0]
    phi = file['tasks']['phi'][0][0][0]
    Lx = int(file['tasks']['Lx'][0][0][0])
    Ly = int(file['tasks']['Ly'][0][0][0])
    Lz = int(file['tasks']['Lz'][0][0][0])
    Nx = int(file['tasks']['Nx'][0][0][0])
    Ny = int(file['tasks']['Ny'][0][0][0])
    Nz = int(file['tasks']['Nz'][0][0][0])
    x = np.linspace(-Lx/2,Lx/2,Nx)
    y = np.linspace(-Ly/2,Ly/2,Ny)

    z_basis = de.Chebyshev('z', Nz, interval=(0,Lz), dealias=3/2)
    z = np.array(z_basis.grid(1))

    print(np.shape(z))
    print("Ra = {}".format(Ra))
    print("Pr = {}".format(Pr))
    print("Ek = {}".format(Ek))
    print("Resolution = ({},{},{})".format(Nx,Ny,Nz))

with h5py.File(direc + "analysis/analysis_" + run_name + ".h5", mode = 'r') as file:
    L_cond = np.array(file['tasks']['L_cond'])[:,0,0,:]
    L_conv = np.array(file['tasks']['L_conv'])[:,0,0,:]
    L_tot = L_cond + L_conv

    KE = np.array(file['tasks']['KE'])[:,0,0,0]
    t = np.array(file['scales']['sim_time'])

with h5py.File(direc + "snapshots/snapshots_" + run_name + ".h5", mode = 'r') as file:
    T = np.array(file['tasks']['T'])
    w = np.array(file['tasks']['w'])

with h5py.File(direc + "averages/averages_" + run_name + ".h5", mode = 'r') as file:
    dz_stress_uw = np.array(file['tasks']['dz_stress_uw'])
    dz_stress_vw = np.array(file['tasks']['dz_stress_vw'])#
    stress_uw = np.array(file['tasks']['stress_uw'])
    stress_vw = np.array(file['tasks']['stress_vw'])
    u_avgt = np.array(file['tasks']['u_avgt'])
    v_avgt = np.array(file['tasks']['v_avgt'])
    w_avgt = np.array(file['tasks']['w_avgt'])

    print(stress_uw.shape, stress_vw.shape)
    print(u_avgt.shape,v_avgt.shape,w_avgt.shape)


if plot_KE:
    fig = plt.figure(figsize = (10,10))
    plt.plot(t,KE)
    plt.xlabel("t")
    plt.ylabel("Kinetic Energy")
    plt.title("Ra = {}, Pr = {}, Ek = {}, phi = {}".format(Ra, Pr, Ek, phi))
    plt.savefig(save_direc + "KE_"+ run_name)
    plt.clf
    plt.close

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

if make_animation:
    def animate(frame):
        quad.set_array(frame)

    fig = plt.figure(figsize=(10,10), dpi=100)
    quad = plt.contourf(x, y, T[0], cmap='CMRmap',levels=200)
    plt.colorbar(label="Temperature")
    plt.title("Ra = {}, Pr = {}, Ek = {}".format(Ra, Pr, Ek))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    animation = ani.FuncAnimation(fig, animate, frames=T[:])
    animation.save(save_direc + 'convection.htm')  

if velocity_snapshot:
    fig = plt.figure(figsize=(10,10), dpi=100)
    plt.contourf(x,y,w[-1,:,:,-2],cmap='bwr',levels=200)
    plt.colorbar(label="Velocity")
    plt.title("Ra = {}, Pr = {}, Ek = {}, phi = {}".format(Ra, Pr, Ek, phi))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig(save_direc+"velocity_slice_"+ run_name)
    plt.clf()
    plt.close()

if temp_snapshot:
    fig = plt.figure(figsize=(10,10), dpi=100)
    plt.contourf(x,y,T[-1,:,:,-2],cmap='hot',levels=200)
    plt.colorbar(label="Temperature")
    plt.title("Ra = {}, Pr = {}, Ek = {}, phi = {}".format(Ra, Pr, Ek, phi))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig(save_direc+"temp_slice_"+ run_name)
    plt.clf()
    plt.close()

if reynolds_stress:

    plot_shape = np.array((2,1))
    plot_size_each = np.array((10,5))

    fig = plt.figure(figsize=np.flip(plot_shape) * plot_size_each)
    fig.suptitle("Ra = {}, Pr = {}, Ek = {}, phi = {}".format(Ra, Pr, Ek, phi))

    ax = fig.add_subplot(*plot_shape, 1)
    ax.set_title("dz(<u'w'>)")
    ctf = ax.contourf(x,z,dz_stress_uw[-1,:,0,:].T,cmap='Spectral',levels=200)
    fig.colorbar(ctf,ax=ax)
    ax.set_xlabel('x')
    ax.set_ylabel('z')

    ax = fig.add_subplot(*plot_shape, 2)
    ax.set_title("dz(<v'w'>)")
    ctf = ax.contourf(x,z,dz_stress_vw[-1,:,0,:].T,cmap='Spectral',levels=200)
    fig.colorbar(ctf,ax=ax)
    ax.set_xlabel('x')
    ax.set_ylabel('z')

    plt.savefig(save_direc + "reynolds_stress_" + run_name)
    plt.tight_layout()
    plt.clf()
    plt.close()

if stress:
    plot_shape = np.array((2,1))
    plot_size_each = np.array((10,5))

    fig = plt.figure(figsize=np.flip(plot_shape) * plot_size_each)
    fig.suptitle("Ra = {}, Pr = {}, Ek = {}, phi = {}".format(Ra, Pr, Ek, phi))

    ax = fig.add_subplot(*plot_shape, 1)
    ax.set_title("<u'w'>")
    ctf = ax.contourf(x,z,stress_uw[-1,:,0,:].T,cmap='Spectral',levels=200)
    fig.colorbar(ctf,ax=ax)
    ax.set_xlabel('x')
    ax.set_ylabel('z')

    ax = fig.add_subplot(*plot_shape, 2)
    ax.set_title("<v'w'>")
    ctf = ax.contourf(x,z,stress_vw[-1,:,0,:].T,cmap='Spectral',levels=200)
    fig.colorbar(ctf,ax=ax)
    ax.set_xlabel('x')
    ax.set_ylabel('z')

    plt.savefig(save_direc + "stress_" + run_name)
    plt.tight_layout()
    plt.clf()
    plt.close()

if mean_flow:
    plot_shape = np.array((3,1))
    plot_size_each = np.array((10,5))

    fig = plt.figure(figsize=np.flip(plot_shape) * plot_size_each)
    fig.suptitle("Ra = {}, Pr = {}, Ek = {}, phi = {}".format(Ra, Pr, Ek, phi))

    ax = fig.add_subplot(*plot_shape, 1)
    ax.set_title("<u>")
    ctf = ax.contourf(x,z,u_avgt[-1,:,0,:].T,cmap='Spectral',levels=200)
    fig.colorbar(ctf,ax=ax)
    ax.set_xlabel('x')
    ax.set_ylabel('z')

    ax = fig.add_subplot(*plot_shape, 2)
    ax.set_title("<v>")
    ctf = ax.contourf(x,z,v_avgt[-1,:,0,:].T,cmap='Spectral',levels=200)
    fig.colorbar(ctf,ax=ax)
    ax.set_xlabel('x')
    ax.set_ylabel('z')

    ax = fig.add_subplot(*plot_shape, 3)
    ax.set_title("<w>")
    ctf = ax.contourf(x,z,w_avgt[-1,:,0,:].T,cmap='Spectral',levels=200)
    fig.colorbar(ctf,ax=ax)
    ax.set_xlabel('x')
    ax.set_ylabel('z')

    plt.savefig(save_direc + "mean_flow_" + run_name)
    plt.tight_layout()
    plt.clf()
    plt.close()

