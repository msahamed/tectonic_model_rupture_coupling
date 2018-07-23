
# importing all the required modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
import errno, sys
import fdfault1 as fdfault
import mesh as mesh


def get_geodynamic_pstrain():

	stress_file = 'long_plastic_strain.csv'
	data = pd.read_csv(stress_file, sep='\t')
	xx  = []
	yy  = []
	pstrain = []

	where_geo_model_begins_x = mesh.start
	where_geo_model_ends_x = mesh.end

	for i in range (len(data.x)):
		if data.x[i] >= (where_geo_model_begins_x) and data.x[i] <= (where_geo_model_ends_x) :
			x = (data.x[i] - where_geo_model_begins_x)/mesh.km
			z = (data.z[i] + mesh.height)/mesh.km
			xx.append(x)
			yy.append(z)
			pstrain.append(data.plastic_strain[i])


	return np.asarray(xx).T , np.asarray(yy).T, np.asarray(pstrain).T

def get_interpolated_pstrain():

	[x , y, pstrain] = get_geodynamic_pstrain()

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
	pstraini = []
	
	new_xx = result.x[:, nby0:nby0+nby1]
	new_yy = result.y[:, nby0:nby0+nby1]
	new_result = result.vx[:, nby0:nby0+nby1]

	for (xpt, ypt) in zip(new_xx.flatten(), new_yy.flatten()):
		dist = np.sqrt((x-xpt)**2+(y-ypt)**2)
		order = dist.argsort()
		f1 = interpolate.SmoothBivariateSpline(x[order[:npts]], y[order[:npts]], pstrain[order[:npts]], kx=3, ky=3)

		pstrainres = f1(xpt, ypt)
		pstraini.append(pstrainres[0][0])


	print("====== done step 2 =======")
	pstraini = np.reshape(np.array(pstraini), np.shape(new_result))
	print(np.shape(pstraini))
	sys.exit()

	return pstraini


refine = mesh.refine
nx1 = mesh.nx1*refine+1 # number of grid point in left block
nx2 = mesh.nx2*refine+1 # number of grid point in right block
nx = nx1+nx2 #total for x coordinates
nby0 = mesh.nby0*refine+1 # for 000 and 100 blocks
nby1 = mesh.nby1*refine+1 # for 000 and 100 blocks
ny = nby1 # for 000 and 100 blocks

plstrain = get_interpolated_pstrain()
interpolated_pstrain = 'interpolated_plastic_strain.csv'
np.savetxt(interpolated_pstrain, plstrain, delimiter=",")
