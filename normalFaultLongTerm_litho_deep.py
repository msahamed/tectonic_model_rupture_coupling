# Geodynamic model description :
import fdfault
import fdfault1 as fdfault2
import sys
import numpy as np
import pandas as pd
from get_surface import read_surface_data
from get_stress import lithostatic_stress
import mesh as mesh
import matplotlib.pyplot as plt


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

surf_height = max([max(y1), max(y2)])

problem_name = 'litho_deep_initial'
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

# # top right surface 110
surf2 = fdfault.curve(nx2, 'y', x2, y2)
p.set_block_surf((1,1,0), 3, surf2)


# fault surface
xcord = x1[-1] # fault bottom x coordinate
ycord = y1[-1] # fault bottom y coordinate

x_fault = np.linspace(block010_length, xcord, nby1)
y_fault = np.linspace(added_layer_thickness, ycord, nby1)

surf3 = fdfault.curve(nby1, 'x', x_fault, y_fault)
p.set_block_surf((0,1,0), 1, surf3)
p.set_block_surf((1,1,0), 0, surf3)


# # block 010 left boundary surface
x = np.zeros(nby1)
y = np.linspace(added_layer_thickness, y1[0], nby1)
surf4 = fdfault.curve(nby1, 'x', x, y)
p.set_block_surf((0,1,0), 0, surf4)

# # block 110 right boundary surface
x = np.ones(nby1)*(block010_length + block110_length)
y = np.linspace(added_layer_thickness, y2[-1], nby1)
surf5 = fdfault.curve(nby1, 'x', x, y)
p.set_block_surf((1,1,0), 1, surf5)

# # block 010 bottom boundary
x = np.linspace(0, np.amin(x_fault), nx1) 
y = np.ones(nx1) * y_fault[0]
surf6 = fdfault.curve(nx1,'y', x, y)
p.set_block_surf((0,1,0), 2, surf6)
p.set_block_surf((0,0,0), 3, surf6)


# # block 010 bottom boundary
x = np.linspace(np.amin(x_fault), block010_length + block110_length, nx2) 
y = np.ones(nx2) * y_fault[0]
surf7 = fdfault.curve(nx2,'y', x, y)
p.set_block_surf((1,1,0), 2, surf7)
p.set_block_surf((1,0,0), 3, surf7)

# block 000 left boundary surface
x = np.zeros(nby0)
y = np.linspace(0, added_layer_thickness, nby0)
surf8 = fdfault.curve(nby0, 'x', x, y)
p.set_block_surf((0,0,0), 0, surf8)

# block 000 and 100 shared boundary surface
x = np.ones(nby0) * x_fault[0]
y = np.linspace(0, added_layer_thickness, nby0)
surf9 = fdfault.curve(nby0, 'x', x, y)
p.set_block_surf((0,0,0), 1, surf9)
p.set_block_surf((1,0,0), 0, surf9)


# # block 100 right boundary surface
x = np.ones(nby0)*(block010_length + block110_length)
y = np.linspace(0, added_layer_thickness, nby0)
surf10 = fdfault.curve(nby0, 'x', x, y)
p.set_block_surf((1,0,0), 1, surf10)

# # block 000 bottom boundary
x = np.linspace(0, x_fault[0], nx1) 
y = np.zeros(nx1)
surf11 = fdfault.curve(nx1,'y', x, y)
p.set_block_surf((0,0,0), 2, surf11)

# # block 000 bottom boundary
x = np.linspace(x_fault[0], block010_length + block110_length, nx2)
y = np.zeros(nx2)
surf12 = fdfault.curve(nx2,'y', x, y)
p.set_block_surf((1,0,0), 2, surf12)

density = 2.7
first_lames_param = 30
shear_modulus = 30

p.set_material(fdfault.material('elastic', density, first_lames_param, shear_modulus))

# set interface types
p.set_iftype(1,'slipweak')

# Add stress
sn, s2, s3 = lithostatic_stress(y_fault)

p.set_loadfile(1, fdfault.loadfile(1, nby1, sn, s2, s3))

# set slip weakening parameters for fault interface-1
p.add_pert(fdfault.swparam('constant', dc = 0.4, mus = 0.681, mud = 0.25),1)
# p.add_pert(fdfault.swparam('boxcar', 0., x0 = 65.0, dx = 0.25, mus = -0.20, c0 = 0.),1)
# p.add_pert(fdfault.swparam('boxcar', 0., x0 = 61.0, dx = 1.0, mus = 10000., c0 = 500.),1)

# add cohesionion to top 5 km of the fault
# cohesion = np.zeros((nby1,1))
# zer = np.zeros((nby1,1))

# def get_cohesion(depth):
# 	c1 = 3.0
# 	c2 = 10.0
# 	y2 = 70.79
# 	y1 = 65.0
# 	# Interpolated cohesion
# 	icohesion = c1 + (depth - y1) * (c2-c1)/(y2-y1)
# 	return icohesion

# for i in range(nby1):
# 	if y_fault[i] >= 65.0:
# 		cohesion[i,0] = get_cohesion(y_fault[i])

# p.set_paramfile(1,fdfault.swparamfile(nby1,1,zer,zer,zer,cohesion,zer,zer))


# Full domain properties
p.add_output(fdfault.output('vxbody','vx',0, nt-1, 50*refine, 0, nx-1, 2*refine, 0, ny-1, 2*refine, 0, 0, 1))
p.add_output(fdfault.output('vybody','vy',0, nt-1, 50*refine, 0, nx-1, 2*refine, 0, ny-1, 2*refine, 0, 0, 1))
p.add_output(fdfault.output('sxbody','sxx',0, nt-1, 50*refine, 0, nx-1, 2*refine, 0, ny-1, 1*refine, 0, 0, 1))
p.add_output(fdfault.output('sxybody','sxy',0, nt-1, 50*refine, 0, nx-1, 2*refine, 0, ny-1, 2*refine, 0, 0, 1))
p.add_output(fdfault.output('sybody','syy',0, nt-1, 50*refine, 0, nx-1, 2*refine, 0, ny-1, 2*refine, 0, 0, 1))

# Only fault properties
p.add_output(fdfault.output('ufault','U', 0,  nt-1, 50*refine, nx1, nx1, 1*refine, nby0, nby0+nby1-1, 1*refine, 0, 0, 1))
p.add_output(fdfault.output('vfault','V', 0,  nt-1, 50*refine, nx1, nx1, 1*refine, nby0, nby0+nby1-1, 1*refine, 0, 0, 1))
p.add_output(fdfault.output('sfault','S', 0,  nt-1, 50*refine, nx1, nx1, 1*refine, nby0, nby0+nby1-1, 1*refine, 0, 0, 1))
p.add_output(fdfault.output('snfault','Sn', 0,  nt-1, 50*refine, nx1, nx1, 1*refine, nby0, nby0+nby1-1, 1*refine, 0, 0, 1))

p.write_input()