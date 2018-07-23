
# importing all the required modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
import errno, sys
import fdfault1 as fdfault
import mesh as mesh


def get_geodynamic_stress():

	stress_file = 'long_stress_'+sys.argv[1]+'.csv'
	data = pd.read_csv(stress_file)
	xx  = []
	yy  = []
	sxx = []
	sxy = []
	syy = []

	where_geo_model_begins_x = mesh.start
	where_geo_model_ends_x = mesh.end

	for i in range (len(data.x)):
		if data.x[i] >= (where_geo_model_begins_x) and data.x[i] <= (where_geo_model_ends_x) :
			x = (data.x[i] - where_geo_model_begins_x)/mesh.km
			z = (data.z[i] + mesh.height)/mesh.km
			xx.append(x)
			yy.append(z)
			sxx.append(data.stress_XX[i]/mesh.pa2mpa)
			syy.append(data.stress_ZZ[i]/mesh.pa2mpa)
			sxy.append(data.stress_XZ[i]/mesh.pa2mpa)

	return np.asarray(xx).T , np.asarray(yy).T, np.asarray(sxx).T, np.asarray(sxy).T, np.asarray(syy).T


def set_initial_stress():

	'''------------- Paramters for this functions -------------
	nx = number of grid points in x-direction
	ny = number of grid points in y-direction
	start = initial position of the domian in Geodynamic model
	end = end position of the domain in geodynamic model
	height  = height of the geodynamic model
	offset = offset is used to start everything from (0,0)
	-----------------------------------------------------------
	'''
	[x , y, sxx, sxy, syy] = get_geodynamic_stress()

	# 'fdfault' generated grid (curvilinear)
	result = fdfault.output('longterm_initial_mesh','vxbody')
	result.load()

	# added_layer_thickness = 60.
	
	# print(nby1, nby0, nby0+nby1)
	
	# plt.pcolor(result.x[:, nby0:nby0+nby1], result.y[:, nby0:nby0+nby1], result.vx[:, nby0:nby0+nby1])
	# plt.scatter(x,y)
	# plt.show()

	# print(np.shape(result.x))
	# print(np.shape(result.y[:,0:nby1]))
	# print(result.y[0,0:nby1])
	# sys.exit()

	print("====== done step 1 =======")

	npts = 30
	sxxi = []
	sxyi = []
	syyi = []
	
	new_xx = result.x[:, nby0:nby0+nby1]
	new_yy = result.y[:, nby0:nby0+nby1]
	new_result = result.vx[:, nby0:nby0+nby1]

	for (xpt, ypt) in zip(new_xx.flatten(), new_yy.flatten()):
		dist = np.sqrt((x-xpt)**2+(y-ypt)**2)
		order = dist.argsort()
		f1 = interpolate.SmoothBivariateSpline(x[order[:npts]], y[order[:npts]], sxx[order[:npts]], kx=3, ky=3)
		f2 = interpolate.SmoothBivariateSpline(x[order[:npts]], y[order[:npts]], sxy[order[:npts]], kx=3, ky=3)
		f3 = interpolate.SmoothBivariateSpline(x[order[:npts]], y[order[:npts]], syy[order[:npts]], kx=3, ky=3)

		sxxres = f1(xpt, ypt)
		sxxi.append(sxxres[0][0])

		sxyres = f2(xpt, ypt)
		sxyi.append(sxyres[0][0])

		syyres = f3(xpt, ypt)
		syyi.append(syyres[0][0])

	print("====== done step 2 =======")
	sxxi = np.reshape(np.array(sxxi), np.shape(new_result))
	sxyi = np.reshape(np.array(sxyi), np.shape(new_result))
	syyi = np.reshape(np.array(syyi), np.shape(new_result))

	return sxxi, sxyi, syyi


# def calculate_lithostatic_stress():
# 	result = fdfault.output('longterm_initial_mesh','vxbody')
# 	result.load()
# 	rho = 2700
# 	gravity = 9.80

# 	ycoord = result.y[:, nby0:nby0+nby1]
# 	xcoord = result.x[:, nby0:nby0+nby1]
# 	# max_depth = np.amax(depth_info)
# 	# actual_depth = depth_info - max_depth
	
# 	xsurf, ysurf = [result.x[:, -1], result.y[:, -1]]
# 	shape = np.shape(result.vx[:, nby0:nby0+nby1])
# 	szz = np.zeros(shape)

# 	for i in range(shape[0]):
# 		for j in range(shape[1]):
# 			current_y_coord = ycoord[i][j]
# 			current_x_coord = xcoord[i][j]
# 			surface_y_coord = yinterp = np.interp(current_x_coord, xsurf, ysurf)
# 			depth = np.absolute(surface_y_coord - current_y_coord)
# 			szz[i][j] = -(rho * gravity * depth*1000)/mesh.pa2mpa

# 	return szz

def lithostatic_stress(yfault):

	rho = 2700
	gravity = 9.80
	fric_coff = 0.35

	sn = np.zeros([1, np.shape(yfault)[0]])
	s2 = np.zeros([1, np.shape(yfault)[0]])
	s3 = np.zeros([1, np.shape(yfault)[0]])

	for i in range(np.shape(yfault)[0]):
		depth = (yfault[i] - np.amax(yfault))*1000
		sn[0][i] = (rho * gravity * depth)/mesh.pa2mpa
		s2[0][i] = sn[0][i] * fric_coff
		s3[0][i] = s2[0][i]
	return sn, s2, s3


# refine = mesh.refine
# nx1 = mesh.nx1*refine+1 # number of grid point in left block
# nx2 = mesh.nx2*refine+1 # number of grid point in right block
# nx = nx1+nx2 #total for x coordinates
# nby0 = mesh.nby0*refine+1 # for 000 and 100 blocks
# nby1 = mesh.nby1*refine+1 # for 000 and 100 blocks
# ny = nby1 # for 000 and 100 blocks

# # sxx, sxy, syy = set_initial_stress()

# litho_szz = calculate_lithostatic_stress()
# print(np.amin(litho_szz), np.amax(litho_szz))

# result = fdfault.output('longterm_initial_mesh','vxbody')
# result.load()
# ycoord = result.y[:, nby0:nby0+nby1]
# xcoord = result.x[:, nby0:nby0+nby1]


# plt.pcolormesh(xcoord, ycoord- np.amax(ycoord), litho_szz, vmin = np.amin(litho_szz), vmax = np.amax(litho_szz))
# plt.colorbar(orientation="vertical")
# plt.xlabel('distance accross fault (km)')
# plt.ylabel('Depth (km)')
# plt.title('lithostatic stress (MPa)')
# plt.tight_layout()
# plt.show()
# sys.exit()

# stress = np.zeros((3, nx, ny))
# stress[0,:,:] = sxx
# stress[1,:,:] = sxy
# stress[2,:,:] = syy
# new_stress = stress.reshape(1, 3*nx*ny)

# interpolated_stress = 'interpolated_stress_'+sys.argv[1]+'.csv'
# np.savetxt(interpolated_stress, new_stress, delimiter=",")
