"""
Parameter file for use in the Dedalus convection script.
"""

import numpy as np

Lx, Ly, Lz = 4., 4., 2.         # Domain size
Nx, Ny, Nz = 48, 48, 24           # Number of gridpoints
Pr = 1.                           # Prandtl number
Ra = 1e6                          # Rayleigh number
Ek = 1e-2                         # Ekman Number
phi = np.pi/4.                    # Angle corresponding to latitude

initial_timestep = 1e-6           # Initial timestep
max_dt = 1e-4                     # max dt

tau = 0.1                        # Time averaging interval

snapshot_freq = 1.5e-3            # Frequency snapshot files are outputted
analysis_freq = 2e-4              # Frequency analysis files are outputted

end_sim_time = 1.                # Stop time in simulations units
end_wall_time = np.inf            # Stop time in wall time
end_iterations = np.inf           # Stop time in iterations

run_name = "1"
