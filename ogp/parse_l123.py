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

fixture_basis = ogp_util.l123_basis()
meas_basis = ogp_util.get_measbasis(stepdict, fixture_basis)

points = ogp_util.get_points(stepdict, near, 7, 4)

#print points
points = ogp_util.transform_pts(meas_basis,fixture_basis,points)
#print points

[fiducials_nominal,p0] = ogp_util.l123_fiducials(side)

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

topedge1 = sensor_origin[1] + 20.17/math.cos(angle) + (6.500*25.4 - sensor_origin[0])*math.tan(angle)
topedge2 = sensor_origin[1] + 20.17/math.cos(angle) + (3.125*25.4 - sensor_origin[0])*math.tan(angle)

if (side=='front'):
	nominal1 = 2.505
	nominal2 = 2.166
else:
	nominal1 = 2.337
	nominal2 = 2.337

print near+' edge 1: nominal {0:.4f}, actual {1:.4f} (deviation {2:.2f} um)'.format(nominal1*25.4, topedge1, 1000*(topedge1-nominal1*25.4))
print near+' edge 2: nominal {0:.4f}, actual {1:.4f} (deviation {2:.2f} um)'.format(nominal2*25.4, topedge2, 1000*(topedge2-nominal2*25.4))

print 'Fit to z = a*x^2 + b*y^2 + c*x*y: a={0}, b={1}, c={2}'.format(sol[6],sol[7],sol[8])

corner1 = sol[6]*(50.0)**2 + sol[7]*(20.17)**2 + sol[8]*50.0*20.17
corner2 = sol[6]*(50.0)**2 + sol[7]*(20.17)**2 - sol[8]*50.0*20.17
print 'Deviation from plane at corners (um): {0:.2f}, {1:.2f}'.format(corner1*1000,corner2*1000)

hess = np.array([[2*sol[6],sol[8]], [sol[8],2*sol[7]]])
w, v = np.linalg.eigh(hess)
print 'Eigenvalues of curvature matrix: {0}'.format(w)
print 'Eigenvectors of curvature matrix: \n{0}'.format(v)

edgeangle = math.radians(stepdict[far+' upper edge']['XY Angle'][1])
edgeel = math.radians(stepdict[far+' upper edge']['Elevation'][1])
edgevec = np.array([math.cos(edgeangle)*math.cos(edgeel),
	math.sin(edgeangle)*math.cos(edgeel),
	math.sin(edgeel)])
edgepos = np.array([stepdict[far+' upper edge']['X Location'][1],
	stepdict[far+' upper edge']['Y Location'][1],
	stepdict[far+' upper edge']['Z Location'][1]])
edgevec = ogp_util.transform_vec(meas_basis,fixture_basis,edgevec)
edgepos = ogp_util.transform_pt(meas_basis,fixture_basis,edgepos)

topedge1 = edgepos[1] + (6.500*25.4 - edgepos[0])*(edgevec[1]/edgevec[0])
topedge2 = edgepos[1] + (3.125*25.4 - edgepos[0])*(edgevec[1]/edgevec[0])

if (side=='front'):
	nominal1 = 2.337
	nominal2 = 2.337
else:
	nominal1 = 2.505
	nominal2 = 2.166

print far+' edge 1: nominal {0:.4f}, actual {1:.4f} (deviation {2:.2f} um)'.format(nominal1*25.4, topedge1, 1000*(topedge1-nominal1*25.4))
print far+' edge 2: nominal {0:.4f}, actual {1:.4f} (deviation {2:.2f} um)'.format(nominal2*25.4, topedge2, 1000*(topedge2-nominal2*25.4))
