#!/usr/bin/env python
import sys,math
import numpy as np
import scipy.optimize
import ogp_util

np.set_printoptions(precision=3)

if (len(sys.argv)!=3):
	print sys.argv[0]+' (front|back) file.PRT'
	sys.exit()
side = sys.argv[1]
if (side=='front'):
	near = 'Stereo'
	far = 'Axial'
elif (side=='back'):
	near = 'Axial'
	far = 'Stereo'
else:
	print 'Invalid side'
	sys.exit()

stepdict = ogp_util.getsteps(sys.argv[2])

fixture_basis = ogp_util.l456_basis()
meas_basis = ogp_util.get_measbasis(stepdict, fixture_basis)

def check_alignment(moduleend):
	sensorname = near+' '+moduleend
	points = np.empty([3,11])
	points = ogp_util.get_points(stepdict, sensorname, 4, 7)

	#print points
	for i in range(0,11):
		points[:,i] = ogp_util.transform_pt(meas_basis,fixture_basis,points[:,i])

	#print points

	[fiducials_nominal,p0] = ogp_util.l456_fiducials(moduleend,side)

	sol = scipy.optimize.leastsq(ogp_util.sensor_fitfunc,p0,args=(points,fiducials_nominal))[0]
	print 'Sensor fit parameters: '+str(sol)
	print 'Sensor fit residuals (um): '+str(1000.0*ogp_util.sensor_fitfunc(sol, points, fiducials_nominal))

	sensor_origin = np.array(sol[0:3])
	sensor_rotation = ogp_util.make_rotation(sol[3],sol[4],sol[5])
	print 'Sensor origin: {0}'.format(sensor_origin)
	norm_vec = sensor_rotation.T.dot([0.0, 0.0, 1.0])
	print 'Normal vector: {0}'.format(norm_vec)
	meas_vec = sensor_rotation.T.dot([1.0, 0.0, 0.0])
	print 'Strip vector: {0}'.format(meas_vec)
	angle = math.atan2(meas_vec[1],meas_vec[0])
	print 'Strip angle: {0}'.format(angle)

	if (moduleend=='block'):
		fixture_x1 = 5.875
		fixture_x2 = 2.500
		if (side=='front'):
			nominal1 = 2.176
			nominal2 = 2.008
		else:
			nominal1 = 2.176
			nominal2 = 2.176
	else:
		fixture_x1 = 9.250
		fixture_x2 = 5.875
		if (side=='front'):
			nominal1 = 2.346
			nominal2 = 2.176
		else:
			nominal1 = 2.176
			nominal2 = 2.176

	topedge1 = sensor_origin[1] + 20.17/abs(math.cos(angle)) + (fixture_x1*25.4 - sensor_origin[0])*math.tan(angle)
	topedge2 = sensor_origin[1] + 20.17/abs(math.cos(angle)) + (fixture_x2*25.4 - sensor_origin[0])*math.tan(angle)

	print near+' '+moduleend+' edge 1: nominal {0:.4f}, actual {1:.4f} (deviation {2:.2f} um)'.format(nominal1*25.4, topedge1, 1000*(topedge1-nominal1*25.4))
	print near+' '+moduleend+' edge 2: nominal {0:.4f}, actual {1:.4f} (deviation {2:.2f} um)'.format(nominal2*25.4, topedge2, 1000*(topedge2-nominal2*25.4))

	corner1 = sol[6]*(50.0)**2 + sol[7]*(20.17)**2 + sol[8]*50.0*20.17
	corner2 = sol[6]*(50.0)**2 + sol[7]*(20.17)**2 - sol[8]*50.0*20.17
	print 'Deviation from plane at corners (um): {0:.2f}, {1:.2f}'.format(corner1*1000,corner2*1000)

	hess = np.array([[2*sol[6],sol[8]], [sol[8],2*sol[7]]])
	w, v = np.linalg.eigh(hess)
	print 'Eigenvalues of curvature matrix:'
	print w
	print 'Eigenvectors of curvature matrix:'
	print v


	if (moduleend=='block'):
		if (side=='front'):
			nominal1 = 2.176
			nominal2 = 2.176
			visible = 'upper'
			edgedistance = 0.0
		else:
			nominal1 = 2.176
			nominal2 = 2.008
			visible = 'lower'
			edgedistance = 40.34
	else:
		if (side=='front'):
			nominal1 = 2.176
			nominal2 = 2.176
			visible = 'lower'
			edgedistance = 40.34
		else:
			nominal1 = 2.346
			nominal2 = 2.176
			visible = 'upper'
			edgedistance = 0.0

	edgeangle = math.radians(stepdict[far+' '+moduleend+' '+visible+' edge']['XY Angle'][1])
	edgeel = math.radians(stepdict[far+' '+moduleend+' '+visible+' edge']['Elevation'][1])
	edgevec = np.array([math.cos(edgeangle)*math.cos(edgeel),
		math.sin(edgeangle)*math.cos(edgeel),
		math.sin(edgeel)])
	edgepos = np.array([stepdict[far+' '+moduleend+' '+visible+' edge']['X Location'][1],
		stepdict[far+' '+moduleend+' '+visible+' edge']['Y Location'][1]+edgedistance/abs(math.cos(edgeangle)),
		stepdict[far+' '+moduleend+' '+visible+' edge']['Z Location'][1]])
	edgevec = ogp_util.transform_vec(meas_basis,fixture_basis,edgevec)
	edgepos = ogp_util.transform_pt(meas_basis,fixture_basis,edgepos)

	topedge1 = edgepos[1] + (fixture_x1*25.4 - edgepos[0])*(edgevec[1]/edgevec[0])
	topedge2 = edgepos[1] + (fixture_x2*25.4 - edgepos[0])*(edgevec[1]/edgevec[0])

	print far+' '+moduleend+' edge 1: nominal {0:.4f}, actual {1:.4f} (deviation {2:.2f} um)'.format(nominal1*25.4, topedge1, 1000*(topedge1-nominal1*25.4))
	print far+' '+moduleend+' edge 2: nominal {0:.4f}, actual {1:.4f} (deviation {2:.2f} um)'.format(nominal2*25.4, topedge2, 1000*(topedge2-nominal2*25.4))

check_alignment('block')
check_alignment('lever')
