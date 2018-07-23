
#Rupture params
refine = 1
# nt = 15000 * refine # Number of timesteps
nt = 100 * refine # Number of timesteps
nx1 = 400*refine+1 # number of grid point in left block
nx2 = 400*refine+1 # number of grid point in right block
nby0 = 600*refine+1 # for 000 and 100 blocks
nby1 = 200*refine+1 # for 010 and 110 blocks

# fault location
# y_fault1 = 0.7686680318511926 * x_fault + 36.506554679323905
fault_surface = 54601. # x coordinate of fault starting point at the surface
fault_bottom = 40580.  # x coordinate of fault starting point at the bottom
# surf_height = 11.43787601

# geodynamic model slice params
start = 10000
end = 80000
km = 1000
pa2mpa = 1e6
height = 70000
