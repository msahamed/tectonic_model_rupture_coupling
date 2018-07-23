# Geodynamic model description :
import fdfault
import fdfault1 as fdfault_out
import sys
import numpy as np
import pandas as pd
from get_surface import read_surface_data
import mesh as mesh
import matplotlib.pyplot as plt
import get_fault_data as ifault

# latex parameter
font = {
    'family': 'serif', 
    'serif': ['Computer Modern Roman'],
    'weight' : 'regular',
    'size'   : 14
    }

plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.figure(figsize=(9, 3))
#-------- Grid set up and import surface data ----------
refine = mesh.refine
nt = mesh.nt*refine # Number of timesteps
nx1 = mesh.nx1*refine+1 # number of grid point in left block
nx2 = mesh.nx2*refine+1 # number of grid point in right block
nx = nx1+nx2 # total for x coordinates
nby0 = mesh.nby0*refine+1 # for 000 and 100 blocks
nby1 = mesh.nby1*refine+1 # for 000 and 100 blocks
ny = nby0 + nby1

[x1, y1, x2, y2, block010_length, block110_length] = read_surface_data(nx1, nx2)
added_layer_thickness = 60.

# # # fault surface
xcord = x1[-1] # fault bottom x coordinate
ycord = y1[-1] # fault bottom y coordinate

x_fault = np.linspace(block010_length, xcord, nby1)
y_fault = np.linspace(added_layer_thickness, ycord, nby1)

ipstrain_file = 'interpolated_plastic_strain.csv'
pstrain = (np.loadtxt(ipstrain_file, delimiter=','))
pstrain = pstrain.reshape(nx, nby1)

result = fdfault_out.output('longterm_initial_mesh','vxbody')
result.load()

im = plt.pcolor(result.x[:, nby0:nby0+nby1], result.y[:, nby0:nby0+nby1]-70, pstrain,
vmin=0, vmax=4.0, cmap= 'inferno')

aa, bb, new_x, new_y, org_x, org_y, plstrain = ifault.plot_fault()

plt.scatter(new_x, new_y, s = 10, edgecolors='b', alpha=0.7)
plt.plot(x_fault, y_fault-70, color='w', linestyle='-', linewidth=2, alpha=0.5)

plt.xlabel('Distance across fault (km)')
plt.ylabel('Depth (km)')
plt.xticks([0, 10, 20, 30, 40, 50, 60, 70])
plt.yticks([0, -5, -10])
plt.colorbar(im, orientation="horizontal", pad = 0.3,
    ticks=[0,1, 2, 3, 4], label='$2^{nd}$ invariant of plastic strain')
plt.tight_layout()
plt.show()

