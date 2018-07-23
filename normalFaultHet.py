# Geodynamic model description :
import fdfault
import numpy as np
from get_surface import read_surface_data
# from get_stress import set_initial_stress
# import stressrotate_v2
import matplotlib.pylab as plt

# create problem

#-------- Grid set up and import surface data ----------

# to increase model resolution, increase refine (refine = 2 doubles resolution)

refine = 1

nt = 7000*refine # Number of timesteps
nx1 = 600*refine+1 # number of grid point in left block
nx2 = 600*refine+1 # number of grid point in right block
nx = nx1+nx2 # total for x coordinates
nby0 = 200*refine+1 # for 000 and 100 blocks

fault_surface = 51.
fault_bottom = 42.
surf_height = 10.

[x1, y1, x2, y2, block000_length, block100_length] = read_surface_data(nx1, nx2, fault_bottom,
                                                                       fault_surface, surf_height)
p = fdfault.problem('normalfault')

# set rk and fd order
p.set_rkorder(4)
p.set_sbporder(4)

# set time step info
p.set_nt(nt)
p.set_cfl(0.3)
p.set_ninfo(100*refine)

# set number of blocks and coordinate information
p.set_nblocks((2,1,1))
p.set_nx_block(([nx1, nx2], [nby0], [1]))

# set block dimensions
p.set_block_lx((0,0,0),(block000_length, surf_height))
p.set_block_lx((1,0,0),(block100_length, surf_height))


# set block boundary conditions
p.set_bounds((0,0,0),['absorbing', 'none', 'absorbing', 'free'])
p.set_bounds((1,0,0),['none', 'absorbing', 'absorbing', 'free'])


# ========== set block surface ===============
# top left surface
surf1 = fdfault.curve(nx1, 'y', x1, y1)
p.set_block_surf((0,0,0), 3, surf1)


# top right surface
surf2 = fdfault.curve(nx2, 'y', x2, y2)
p.set_block_surf((1,0,0), 3, surf2)

# fault surface
xcord = x1[-1]
ycord = y1[-1]

x_fault = np.linspace(block000_length, xcord, nby0)
y_fault = np.linspace(0., ycord, nby0)
fault_bottom = x_fault[0]

surf3 = fdfault.curve(nby0, 'x', x_fault, y_fault)
p.set_block_surf((0,0,0), 1, surf3)
p.set_block_surf((1,0,0), 0, surf3)

# block 000 left boundary surface
x = np.zeros(nby0)
y = np.linspace(0, y1[0], nby0)
surf4 = fdfault.curve(nby0, 'x', x, y)
p.set_block_surf((0,0,0), 0, surf4)

# block 100 right boundary surface
x = np.ones(nby0)*(block000_length + block100_length)
y = np.linspace(0, y2[-1], nby0)
surf5 = fdfault.curve(nby0, 'x', x, y)
p.set_block_surf((1,0,0), 1, surf5)

# ------------ set initial stress -------------
# sn, st = stressrotate_v2.calculate_fault_normal_shear_stress(sxx[nx1,:], sxy[nx1,:], syy[nx1,:], x_fault, y_fault, 'x')
# sxx, sxy, syy = set_initial_stress(20,90,20, 10)
#
# stress = np.zeros((3,(nx1+nx2),nby0))
# stress[0,:,:] = sxx
# stress[1,:,:] = sxy
# stress[2,:,:] = syy
#
p.set_het_stress(stress)

# set interface types
p.set_iftype(0,'slipweak')

# set slip weakening parameters
p.add_pert(fdfault.swparam('constant', dc = 0.4, mus = 0.677, mud = 0.44),0)
p.add_pert(fdfault.swparam('boxcar', x0 = 0.75, dx = 0.75, mus = 10000.),0)
#p.add_pert(fdfault.swparam('boxcar',0., 5., 0.25, 0., 0., 0.,-0.15, 0.),0)

# add load perturbations
p.add_load(fdfault.load('boxcar',0., 5.,0.25, 0., 0., 0.,-16.75, 0.),0)

# add cohesion to free surface
cohes = np.zeros((nby0,1))
zer = np.zeros((nby0,1))

dcohes = 3.
cmax = 20.

for i in range(nby0):
	if y_fault[i] > ycord-dcohes:
		cohes[i,0] = cmax*(y_fault[i]-(ycord-dcohes))/dcohes
	if y_fault[i] < 0.4:
		cohes[i,0] = cmax

for i in range(nby0):
	if y_fault[i] > 7.6:
		cohes[i,0] += cmax

p.set_paramfile(0,fdfault.swparamfile(nby0,1,zer,zer,zer,cohes,zer,zer))
p.set_stress

plt.plot(cohes[:,0], y_fault)
plt.show()

fdfault.loadfile()


# add output unit
p.add_output(fdfault.output('vxbody','vx',0, nt-1, 50*refine, 0, nx-1, 2*refine, 0, nby0-1, 2*refine, 0, 0, 1))
p.add_output(fdfault.output('vybody','vy',0, nt-1, 50*refine, 0, nx-1, 2*refine, 0, nby0-1, 2*refine, 0, 0, 1))
p.add_output(fdfault.output('vfault','V', 0,  nt-1, 50*refine, nx1, nx1, refine, 0, nby0-1, 2*refine, 0, 0, 1))
p.add_output(fdfault.output('sxbody','sxx',0, nt-1, 50*refine, 0, nx-1, 2*refine, 0, nby0-1, 1*refine, 0, 0, 1))
p.add_output(fdfault.output('sxybody','sxy',0, nt-1, 50*refine, 0, nx-1, 2*refine, 0, nby0-1, 2*refine, 0, 0, 1))
p.add_output(fdfault.output('sybody','syy',0, nt-1, 50*refine, 0, nx-1, 2*refine, 0, nby0-1, 2*refine, 0, 0, 1))
p.add_output(fdfault.output('sfault','S', 0,  nt-1, 50*refine, nx1, nx1, refine, 0, nby0-1, 2*refine, 0, 0, 1))
p.add_output(fdfault.output('snfault','Sn', 0,  nt-1, 50*refine, nx1, nx1, refine, 0, nby0-1, 2*refine, 0, 0, 1))

# -- add output for seismogram
# off fault stations
offfault = [('000', '000', '000'), ('100', '000', '000'), ('200', '000', '000'), ('300', '000', '000'), ('400', '000', '000'),('500', '000', '000'),('600', '000', '000'), ('700', '000', '000'), ('800', '000', '000')]

fields = ['h-vel','n-vel']
fname = ['vx','vy']

for station in offfault:
    xcoord = float(station[0])/10.
    ycoord = float(station[1])/10.
    zcoord = float(station[2])/10.
    xpt, ypt, zpt = p.find_nearest_point((xcoord, ycoord, zcoord))
    for fld, fn in zip(fields, fname):
        p.add_output(fdfault.output('body'+station[0]+'km_off'+'-'+fld, fn, 0, nt, 1, xpt, xpt, 1,
                                    ypt, ypt, 1, zpt, zpt, 1))



p.write_input()
