# Geodynamic model description :
import fdfault
import fdfault1 as fdfault2
import sys
import numpy as np
import pandas as pd
from get_surface import read_surface_data
import mesh as mesh
import matplotlib.pyplot as plt

# latex parameter
font = {
    'family': 'serif', 
    'serif': ['Computer Modern Roman'],
    'weight' : 'regular',
    'size'   : 14
    }

plt.rc('font', **font)
plt.rc('text', usetex=True)

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

import pickle
mydict = {'x1':x1, 'x2':x2, 'y1':y1-70, 'y2':y2-70}
output = open('surface.pkl', 'wb')


surf_height = max([max(y1), max(y2)])

arg = str(50)
problem_name = 'longterm_' + arg

p = fdfault.problem(problem_name)
# set rk and fd order
p.set_rkorder(4)
p.set_sbporder(4)

# set time step info
p.set_nt(nt)
p.set_cfl(0.3)
p.set_ninfo(100*refine)

# set number of blocks and coordinate information
p.set_nblocks((2,2,1))
p.set_nx_block(([nx1, nx2], [nby0, nby1], [1]))

# set block dimensions
added_layer_thickness = 60.
block000_height = added_layer_thickness
block100_height = added_layer_thickness
block010_height = mesh.height/mesh.km - added_layer_thickness + np.absolute((mesh.height/mesh.km) - surf_height)
block110_height = mesh.height/mesh.km - added_layer_thickness + np.absolute((mesh.height/mesh.km) - surf_height)

p.set_block_lx((0,0,0),(block010_length, block000_height))
p.set_block_lx((1,0,0),(block110_length, block100_height))
p.set_block_lx((0,1,0),(block010_length, block010_height))
p.set_block_lx((1,1,0),(block110_length, block110_height))

# set block boundary conditions
p.set_bounds((0,0,0),['absorbing', 'none', 'absorbing', 'none'])
p.set_bounds((1,0,0),['none', 'absorbing', 'absorbing', 'none'])
p.set_bounds((0,1,0),['absorbing', 'none', 'none', 'free'])
p.set_bounds((1,1,0),['none', 'absorbing', 'none', 'free'])


# ========== set block surface ===============
# top left surface 010
surf1 = fdfault.curve(nx1, 'y', x1, y1)
p.set_block_surf((0,1,0), 3, surf1)

plt.plot(x1, y1-70, 'b-')

# # top right surface 110
surf2 = fdfault.curve(nx2, 'y', x2, y2)
p.set_block_surf((1,1,0), 3, surf2)

plt.plot(x2, y2-70, 'b-')

# fault surface
xcord = x1[-1] # fault bottom x coordinate
ycord = y1[-1] # fault bottom y coordinate

x_fault = np.linspace(block010_length, xcord, nby1)
y_fault = np.linspace(added_layer_thickness, ycord, nby1)
# y_fault1 = 0.7686680318511926 * x_fault + 36.506554679323905

plt.plot(x_fault, y_fault-70, color='red', linestyle='--', linewidth=3)
# mydict['xfault'] = x_fault
# mydict['yfault'] = y_fault

surf3 = fdfault.curve(nby1, 'x', x_fault, y_fault)
p.set_block_surf((0,1,0), 1, surf3)
p.set_block_surf((1,1,0), 0, surf3)

# # block 010 left boundary surface
x = np.zeros(nby1)
y = np.linspace(added_layer_thickness, y1[0], nby1)
surf4 = fdfault.curve(nby1, 'x', x, y)
p.set_block_surf((0,1,0), 0, surf4)

plt.plot(x, y-70, 'k')

# # block 110 right boundary surface
x = np.ones(nby1)*(block010_length + block110_length)
y = np.linspace(added_layer_thickness, y2[-1], nby1)
surf5 = fdfault.curve(nby1, 'x', x, y)
p.set_block_surf((1,1,0), 1, surf5)

plt.plot(x, y-70, 'k')

# # block 010 bottom boundary
x = np.linspace(0, np.amin(x_fault), nx1) 
y = np.ones(nx1) * y_fault[0]
surf6 = fdfault.curve(nx1,'y', x, y)
p.set_block_surf((0,1,0), 2, surf6)
p.set_block_surf((0,0,0), 3, surf6)

plt.plot(x, y-70, 'g')


# # block 110 bottom boundary
x = np.linspace(np.amin(x_fault), block010_length + block110_length, nx2) 
y = np.ones(nx2) * y_fault[0]
surf7 = fdfault.curve(nx2,'y', x, y)
p.set_block_surf((1,1,0), 2, surf7)
p.set_block_surf((1,0,0), 3, surf7)

plt.plot(x, y-70, 'g')

# block 000 left boundary surface
x = np.zeros(nby0)
y = np.linspace(0, added_layer_thickness, nby0)
surf8 = fdfault.curve(nby0, 'x', x, y)
p.set_block_surf((0,0,0), 0, surf8)

plt.plot(x, y-70, 'k')

# block 000 and 100 shared boundary surface
x = np.ones(nby0) * x_fault[0]
y = np.linspace(0, added_layer_thickness, nby0)
surf9 = fdfault.curve(nby0, 'x', x, y)
p.set_block_surf((0,0,0), 1, surf9)
p.set_block_surf((1,0,0), 0, surf9)

# plt.plot(x, y, 'k')

# # block 100 right boundary surface
x = np.ones(nby0)*(block010_length + block110_length)
y = y = np.linspace(0, added_layer_thickness, nby0)
surf10 = fdfault.curve(nby0, 'x', x, y)
p.set_block_surf((1,0,0), 1, surf10)

plt.plot(x, y-70, 'k')

# # block 000 bottom boundary
x = np.linspace(0, x_fault[0], nx1) 
y = np.zeros(nx1)
surf11 = fdfault.curve(nx1,'y', x, y)
p.set_block_surf((0,0,0), 2, surf11)

plt.plot(x, y-70, 'k')

# # block 000 bottom boundary
x = np.linspace(x_fault[0], block010_length + block110_length, nx2)
y = np.zeros(nx2)
surf12 = fdfault.curve(nx2,'y', x, y)
p.set_block_surf((1,0,0), 2, surf12)


plt.plot(x, y-70, 'k')
plt.xlabel('Distance across fault (km)')
plt.ylabel('Depth (km)')

plt.text(20, 0.90, 'Free surface', color= 'blue')

plt.text(20, -15, 'Lower boundary of TM', color= 'green')

plt.text(20, -68, 'Absorbing boundary', color= 'black')
plt.text(1.5, -30, 'Absorbing', rotation=90, color= 'black')
plt.text(66, -30, 'Absorbing', rotation=90, color= 'black')

plt.text(55, -5, 'Footwall')
plt.text(5, -5, 'Hanging wall')
plt.annotate('Fault', xy=(35, -6.5), xytext=(45, -8),
            arrowprops=dict(facecolor='black', shrink=0.005), color='red')
plt.tight_layout()
plt.show()