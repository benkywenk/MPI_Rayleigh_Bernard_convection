"""
Parameter file for use in the Dedalus convection script.
"""

import numpy as np

Lx, Ly, Lz = 4., 4., 1.                       # Domain size
Nx, Ny, Nz = 64, 64, 32                    # Number of
Pr = 1.                             # Prandtl number
Pm = 1.                             # Magnetic Prandtl number
Ra = 1e7                          # Rayleigh number
Ek = 2.5e-3                           #Ekman Number
phi = np.pi/6.                      #Angle corresponding to latitude

initial_timestep = 1e-6                 # Initial timestep
max_dt = 1e-4                         # max dt

snapshot_freq = 1.5e-3              # Frequency snapshot files are outputted
analysis_freq = 2e-4              # Frequency analysis files are outputted

end_sim_time = 0.5                 # Stop time in simulations units
end_wall_time = np.inf              # Stop time in wall time
end_iterations = np.inf             # Stop time in iterations

run_name = "1"
