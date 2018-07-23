
# Sabber Ahamed
# Graduate Research Assistant
# Center for Earthquake Research and Information (CERI)
# The University of Memphis
# msahamed@memphis.edu

# modified by Eric Daub to create regular grid for the free surface

# Purpose of this code :
# This code can import all the coordinates generated from
# dynearthsol3D-> geodynamic and the convert it to fdfault
# format-> rupture simulation

# importing all the required modules
import numpy as np
from scipy import interpolate
import matplotlib.pylab as plt
import mesh as mesh
import sys

# this function reads the data from csv file and prepare for simulation
def read_surface_data(nx1, nx2):
	"""
	create surfaces for rupture model by interpolating geodynamic data

	parameters:
	nx1 (int) = number of grid points in left block
	nx2 (int) = number of grid points in right block
	fault_bottom (float) = horizontal coordinate of down dip end of fault
	surface_intersection (float) = horizontal coordinate of fault surface intersection
	surf_height (float) = constant offset to be added to free surface height

	returns:
	x1, y1 (arrays of length nx1) = free surface coordinates of left block
	x2, y2 (arrays of length nx2) = free surface coordinates of right block
	seg1_length, seg2_length (floats) = horizontal length of left and right blocks
	"""

	# This section reads geodynamic model surface
	datax, dataz = np.loadtxt('long_surface_50.csv', delimiter = ',', skiprows = 1, unpack = True)
	where_geo_model_begins_x = mesh.start
	where_geo_model_ends_x = mesh.end
	seg1 = (mesh.fault_bottom - where_geo_model_begins_x)/mesh.km
	seg2 = (where_geo_model_ends_x - mesh.fault_bottom)/mesh.km

	xx = [0.0]
	zz = [mesh.height/mesh.km]
	length = len(datax)

	# print(seg2_length)

	for i in range (length):
		if datax[i] >= (where_geo_model_begins_x) and datax[i] <= (where_geo_model_ends_x) :
			x = (datax[i] - where_geo_model_begins_x)/mesh.km
			z = (dataz[i] + mesh.height)/mesh.km
			xx.append(x)
			zz.append(z)

	xx.append(seg1+seg2)
	zz.append(mesh.height/mesh.km)
	# plt.scatter(xx, zz)

	f = interpolate.interp1d(xx, zz)
	x1 = np.linspace(min(xx), (mesh.fault_surface - where_geo_model_begins_x)/mesh.km, nx1)
	x2 = np.linspace((mesh.fault_surface - where_geo_model_begins_x)/mesh.km, seg1+seg2, nx2)

	y1 = f(x1)
	y2 = f(x2)
	
	# plt.plot(x1, y1, x2, y2)
	# plt.show()
	# sys.exit()

	return x1, y1, x2, y2, seg1, seg2
