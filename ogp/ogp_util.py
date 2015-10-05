#!/usr/bin/env python
import re
import numpy as np
import math
import scipy
import csv

def countruns(filename):
	runcount=0
	f = open(filename)
	line = f.readline()
	while len(line)!=0:
		if line[49:54]=='Run #':
			runcount+=1
		line = f.readline()
	f.close()
	print '{0} contains {1} runs'.format(filename,runcount)
	return runcount


def getsteps(filename):
	runcount = countruns(filename)
	return readrun(filename,runcount)

def readrun(filename,runnum):
	runcount = 0

	f = open(filename)
	line = f.readline()
	while True:
		line = f.readline()
		if line[49:54]=='Run #':
			runcount+=1
			if runcount==runnum:
				break

	line = f.readline()
	line = f.readline()
	run = int(line[49:54])
	print 'Run {0}, {1}'.format(run,line[55:79])

	stepdict = dict()

	while True:
		line = f.readline()
		while (len(line)!=0 and (line[0:5] != 'Step ' and line[49:54]!='Run #')):
			line = f.readline()
		if (len(line)==0 or line[49:54]=='Run #'):
			break
		tokens = line.split(None,3)
		step=int(tokens[1])
		if (len(tokens)<4):
			stepname='unnamed'
		else:
			stepname = tokens[3].strip()
		#print 'Step {0}: {1}'.format(step,stepname)
		measdict = dict()
		while True:
			line = f.readline()
			if re.match(' \w\w ',line[13:17]) == None:
				break
			name = line[0:14].rstrip()
			unit = line[13:17].strip()
			nominal = line[17:27].strip()
			if len(nominal)==0:
				nominal='NaN'
			nominal = float(nominal)
			actual = line[28:38]
			if len(actual)==0:
				actual='NaN'
			actual = float(actual)
			#print '{0} ({1}): nominal {2}, actual {3}'.format(name,unit,nominal,actual)
			measdict[name] = [nominal,actual]
		if (stepname in stepdict and stepname != 'unnamed'):
			print 'error!'
		stepdict[stepname]=measdict
	f.close()
	return stepdict

	#print stepdict['Cold block face']['Elevation']

def getsteps_box_survey(filename):
	f = open(filename)
	stepdict = dict()
	reader = csv.DictReader(f)
	for row in reader:
		name = row['Name']
		x = float(row['X (mm)'])
		y = float(row['Y (mm)'])
		z = float(row['Z (mm)'])
		#print [name,x,y,z]
		stepdict[name] = np.array([x,y,z])
	return stepdict


def orthogonalize(vec1,vec2): #return vector orthogonal to vec1, in the plane of vec1 and vec2
	return vec2 - (vec1.dot(vec2)/vec1.dot(vec1))*vec1

def make_basis(oriball,axiball,diaball):
	basis = np.empty([3,3])

	basis[0] = axiball-oriball
	basis[1] = diaball-oriball
	basis[1] = orthogonalize(basis[0],basis[1])
	basis[2] = np.cross(basis[0],basis[1])
	#basis[1] = np.cross(basis[2],basis[0])

	basis[0] = basis[0]/np.linalg.norm(basis[0])
	basis[1] = basis[1]/np.linalg.norm(basis[1])
	basis[2] = basis[2]/np.linalg.norm(basis[2])
	return [oriball,basis,oriball,axiball,diaball]

def make_pin_basis(holepin,slotpin,normal):
	basis = np.empty([3,3])
	basis[0] = slotpin-holepin
	basis[1] = normal
	basis[1] = orthogonalize(basis[0],basis[1])
	basis[2] = np.cross(basis[0],basis[1])

	basis[0] = basis[0]/np.linalg.norm(basis[0])
	basis[1] = basis[1]/np.linalg.norm(basis[1])
	basis[2] = basis[2]/np.linalg.norm(basis[2])
	return [holepin,basis,holepin,slotpin,normal]

def make_uch_basis(holeupst,holedownst,slotupst,slotdownst):
	upst_mid = (holeupst+slotupst)/2.0
	downst_mid = (holedownst+slotdownst)/2.0
	upst_vec = (slotupst-holeupst)
	downst_vec = (slotdownst-holedownst)
	upst_vec = upst_vec/np.linalg.norm(upst_vec)
	downst_vec = downst_vec/np.linalg.norm(downst_vec)
	basis = np.empty([3,3])

	basis[0] = upst_vec+downst_vec
	basis[2] = downst_mid-upst_mid
	basis[2] = orthogonalize(basis[0],basis[2])
	basis[1] = np.cross(basis[2],basis[0])

	basis[0] = basis[0]/np.linalg.norm(basis[0])
	basis[1] = basis[1]/np.linalg.norm(basis[1])
	basis[2] = basis[2]/np.linalg.norm(basis[2])
	return [upst_mid,basis,holeupst,holedownst,slotupst,slotdownst]

def make_rotation(ux, uy, uz):
	rotvec = np.array([ux,uy,uz])
	angle = np.linalg.norm(rotvec)
	if (angle==0):
		return np.identity(3)
	rotvec = rotvec/angle
	matrix = math.cos(angle)*np.identity(3) +\
			math.sin(angle)*np.cross(np.identity(3),rotvec) +\
			(1.0-math.cos(angle))*np.outer(rotvec, rotvec)
	#print np.cross(np.identity(3),rotvec)
	#print np.outer(rotvec, rotvec)
	return matrix

def get_uchbasis(stepdict,uch,pos_name):
	uch_basis = make_uch_basis(stepdict[uch+'1'+pos_name],
				stepdict[uch+'2'+pos_name],
				stepdict[uch+'3'+pos_name],
				stepdict[uch+'4'+pos_name])
	return uch_basis


def get_measbasis(stepdict,fixture_basis):
	oriball = np.empty([3])
	axiball = np.empty([3])
	diaball = np.empty([3])
	oriball[0]=stepdict['Origin ball reference']['X Location'][1]
	oriball[1]=stepdict['Origin ball reference']['Y Location'][1]
	oriball[2]=stepdict['Origin ball reference']['Z Location'][1]
	axiball[0]=stepdict['Axis ball reference']['X Location'][1]
	axiball[1]=stepdict['Axis ball reference']['Y Location'][1]
	axiball[2]=stepdict['Axis ball reference']['Z Location'][1]
	diaball[0]=stepdict['Diagonal ball reference']['X Location'][1]
	diaball[1]=stepdict['Diagonal ball reference']['Y Location'][1]
	diaball[2]=stepdict['Diagonal ball reference']['Z Location'][1]

	meas_basis = make_basis(oriball,axiball,diaball)
	#print fixture_basis
	#print meas_basis



	#print 'Nominal and measured origin, axis and diagonal ball positions in fixture basis:'

	#print 'Origin ball: nominal {0}, actual {1} (deviation {2} um)'.format(fixture_basis[2], transform_pt(meas_basis,fixture_basis,oriball), (transform_pt(meas_basis,fixture_basis,oriball)-fixture_basis[2])*1000)
	#print 'Axis ball: nominal {0}, actual {1} (deviation {2} um)'.format(fixture_basis[3], transform_pt(meas_basis,fixture_basis,axiball), (transform_pt(meas_basis,fixture_basis,axiball)-fixture_basis[3])*1000)
	#print 'Diagonal ball: nominal {0}, actual {1} (deviation {2} um)'.format(fixture_basis[4], transform_pt(meas_basis,fixture_basis,diaball), (transform_pt(meas_basis,fixture_basis,diaball)-fixture_basis[4])*1000)
	return meas_basis

def get_points(stepdict, name, n_fiducial, n_focus):
	points = np.empty([n_fiducial+n_focus,3])
	for i in range(0, n_fiducial):
		points[i,0] = stepdict[name+' fiducial '+str(i+1)+', XY']['X Location'][1]
		points[i,1] = stepdict[name+' fiducial '+str(i+1)+', XY']['Y Location'][1]
		points[i,2] = stepdict[name+' fiducial '+str(i+1)+', Z']['Z Location'][1]
	for i in range(0, n_focus):
		points[n_fiducial+i,0] = stepdict[name+' focus '+str(i+1)]['X Location'][1]
		points[n_fiducial+i,1] = stepdict[name+' focus '+str(i+1)]['Y Location'][1]
		points[n_fiducial+i,2] = stepdict[name+' focus '+str(i+1)]['Z Location'][1]
	return points

def l123_fiducials(side):
	if (side=='front'):
		fiducials_nominal = np.array([[0.0, 22.32, 49.8, 49.8, 49.8, 22.32, 0.0],
			[-19.97, -19.97, -19.97, 0.0, 19.97, 19.97, 19.97],
			[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]).T
		p0 = [124.0, 0.0, 39.0, 1.57, 0.0, 0.0, 0.0, 0.0, 0.0]
	else:
		fiducials_nominal = np.array([[0.0, 22.32, 49.8, 49.8, 49.8, 22.32, 0.0],
			[19.97, 19.97, 19.97, 0.0, -19.97, -19.97, -19.97],
			[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]).T
		p0 = [124.0, 0.0, 39.0, -1.57, 0.0, 0.0, 0.0, 0.0, 0.0]
	return [fiducials_nominal, p0]

def l456_fiducials(moduleend, side):
	if (moduleend=='block'):
		if (side=='front'):
			fiducials_nominal = np.array([[0.0, 22.32, 22.32, 0.0],
				[-19.97, -19.97, 19.97, 19.97],
				[0.0, 0.0, 0.0, 0.0]]).T
			p0 = [124.0, 0.0, 39.0, math.pi/2, 0.0, 0.0, 0.0, 0.0, 0.0]
		else:
			fiducials_nominal = np.array([[22.32, 0.0, 0.0, 22.32],
				[19.97, 19.97, -19.97, -19.97],
				[0.0, 0.0, 0.0, 0.0]]).T
			p0 = [124.0, 0.0, 39.0, -math.pi/2, 0.0, 0.0, 0.0, 0.0, 0.0]
	else:
		if (side=='front'):
			fiducials_nominal = np.array([[22.32, 0.0, 0.0, 22.32],
				[19.97, 19.97, -19.97, -19.97],
				[0.0, 0.0, 0.0, 0.0]]).T
			p0 = [124.0, 0.0, 39.0, 0.0, math.pi/math.sqrt(2), math.pi/math.sqrt(2), 0.0, 0.0, 0.0]
		else:
			fiducials_nominal = np.array([[0.0, 22.32, 22.32, 0.0],
				[-19.97, -19.97, 19.97, 19.97],
				[0.0, 0.0, 0.0, 0.0]]).T
			p0 = [124.0, 0.0, 39.0, 0.0, -math.pi/math.sqrt(2), math.pi/math.sqrt(2), 0.0, 0.0, 0.0]
	return [fiducials_nominal, p0]

def get_l123_sensor(filename,side):
	if (side=='front'):
		near = 'Stereo'
		far = 'Axial'
	elif (side=='back'):
		near = 'Axial'
		far = 'Stereo'
	stepdict = getsteps(filename)

	fixture_basis = l123_basis()
	meas_basis = get_measbasis(stepdict, fixture_basis)

	points = get_points(stepdict, near, 7, 4)
	points = transform_pts(meas_basis,fixture_basis,points)

	[fiducials_nominal,p0] = l123_fiducials(side)

	sol = scipy.optimize.leastsq(sensor_fitfunc,p0,args=(points,fiducials_nominal))[0]
	#print sol
	#print sensor_fitfunc(sol,points,fiducials_nominal)
	return sol

def get_l456_sensor(stepdict,side,moduleend):
	if (side=='front'):
		near = 'Stereo'
		far = 'Axial'
	elif (side=='back'):
		near = 'Axial'
		far = 'Stereo'

	fixture_basis = l456_basis()
	meas_basis = get_measbasis(stepdict, fixture_basis)

	sensorname = near+' '+moduleend
	points = get_points(stepdict, sensorname, 4, 7)
	points = transform_pts(meas_basis,fixture_basis,points)

	[fiducials_nominal,p0] = l456_fiducials(moduleend,side)

	sol = scipy.optimize.leastsq(sensor_fitfunc,p0,args=(points,fiducials_nominal))[0]
	#print sensor_fitfunc(sol,points,fiducials_nominal)
	return sol

def get_l456_sensors(filename,side):
	stepdict = getsteps(filename)
	return [get_l456_sensor(stepdict,side,'block'), get_l456_sensor(stepdict,side,'lever')]

def get_uch_mount(stepdict,layer):
	holepin = np.empty([3])
	slotpin = np.empty([3])
	normal = np.empty([3])
	holepin[0]=stepdict[layer+' hole pin intersection']['X Location'][1]
	holepin[1]=stepdict[layer+' hole pin intersection']['Y Location'][1]
	holepin[2]=stepdict[layer+' hole pin intersection']['Z Location'][1]
	slotpin[0]=stepdict[layer+' slot pin intersection']['X Location'][1]
	slotpin[1]=stepdict[layer+' slot pin intersection']['Y Location'][1]
	slotpin[2]=stepdict[layer+' slot pin intersection']['Z Location'][1]
	az = math.radians(stepdict[layer+' plane, touch']['XY Angle'][1])
	el = math.radians(stepdict[layer+' plane, touch']['Elevation'][1])
	if (el<0):
		az = az+math.pi
		el = -1.0*el
	normal = np.array([math.cos(az)*math.cos(el),
		math.sin(az)*math.cos(el),
		math.sin(el)])
	#print l1normal
	basis = make_pin_basis(holepin, slotpin, normal)
	return basis

def get_sensor_points(stepdict,sensor,n_points):
	points = np.empty([n_points,3])
	for i in range(0,n_points):
		points[i,0] = stepdict[sensor+' '+str(i+1)+', XY']['X Location'][1]
		points[i,1] = stepdict[sensor+' '+str(i+1)+', XY']['Y Location'][1]
		points[i,2] = stepdict[sensor+' '+str(i+1)+', Z']['Z Location'][1]
	return points



def sensor_fitfunc(p, points, nominals):
	num_fiducials = nominals.shape[0]
	points = np.copy(points)
	#print points
	#origin = np.array([p[0],p[1],p[2]])
	points[:,0] = points[:,0] - p[0]
	points[:,1] = points[:,1] - p[1]
	points[:,2] = points[:,2] - p[2]
	#print points
	rotation = make_rotation(p[3],p[4],p[5])
	points = points.dot(rotation)
	#print points
	xres = points[0:num_fiducials,0]-nominals[:,0]
	yres = points[0:num_fiducials,1]-nominals[:,1]
	zres = points[:,2]-\
			p[6]*points[:,0]*points[:,0]+\
			p[7]*points[:,1]*points[:,1]+\
			p[8]*points[:,0]*points[:,1]

	#print points
	#print xres
	#print yres
	#print zres
	return np.concatenate((xres,yres,zres))

def l123_basis():
	#oriball = np.array([252.8540, 22.3107, 53.4938])
	#axiball = np.array([-62.0937, 22.4465, 53.4165])
	#diaball = np.array([252.8277, -2.8950, 28.0921])
	oriball = np.array([252.8540, 53.4938, -22.3107])
	axiball = np.array([-62.0937, 53.4165, -22.4465])
	diaball = np.array([252.8277, 28.0921, 2.8950])
	basis = make_basis(oriball,axiball,diaball)
	return basis

def l456_basis():
	#oriball = np.array([360.6875, 22.4470, 49.2446])
	#axiball = np.array([-62.1860, 22.2939, 49.3544])
	#diaball = np.array([360.7389, -2.8264, 23.8597])
	oriball = np.array([360.6875, 49.2446, -22.4470])
	axiball = np.array([-62.1860, 49.3544, -22.2939])
	diaball = np.array([360.7389, 23.8597, 2.8264])
	basis = make_basis(oriball,axiball,diaball)
	return basis

def t13_global():
	holeupst = np.array([-6.5061, 0.3326, 14.8955])*25.4
	holedownst = np.array([-6.2632, 0.3329, 22.7687])*25.4
	slotupst = np.array([2.8118, 0.3304, 14.6071])*25.4
	slotdownst = np.array([3.0541, 0.3361, 22.4787])*25.4
	return make_uch_basis(holeupst,holedownst,slotupst,slotdownst)

def b13_global():
	holeupst = np.array([-6.4869, -0.3316, 15.5152])*25.4
	holedownst = np.array([-6.2483, -0.3326, 23.3874])*25.4
	slotupst = np.array([2.8211, -0.3315, 15.2338])*25.4
	slotdownst = np.array([3.0593, -0.3298, 23.1039])*25.4
	return make_uch_basis(holeupst,holedownst,slotupst,slotdownst)

def t46_global():
	holeupst = np.array([-7.0295, 0.3328, 30.6691])*25.4
	holedownst = np.array([-6.5571, 0.3320, 46.4145])*25.4
	slotupst = np.array([6.5326, 0.3335, 30.2593])*25.4
	slotdownst = np.array([7.0089, 0.3333, 46.0044])*25.4
	return make_uch_basis(holeupst,holedownst,slotupst,slotdownst)

def b46_global():
	holeupst = np.array([-7.0144, -0.3307, 31.2933])*25.4
	holedownst = np.array([-6.5379, -0.3301, 47.0366])*25.4
	slotupst = np.array([6.5437, -0.3316, 30.8827])*25.4
	slotdownst = np.array([7.0210, -0.3313, 46.6311])*25.4
	return make_uch_basis(holeupst,holedownst,slotupst,slotdownst)

def make_parametrized_basis(x,y,z,ux,uy,uz):
	origin = np.array([x,y,z])
	rotation = make_rotation(ux,uy,uz)
	return [origin,rotation]

def make_height_basis(y,tilt,roll):
	origin = np.array([0.0,y,0.0])
	#rotation = make_rotation(0.0,0.0,roll).dot(make_rotation(tilt,0.0,0.0))
	#rotation = make_rotation(tilt,0.0,0.0).dot(make_rotation(0.0,0.0,roll))
	rotation = make_rotation(tilt,0.0,roll)
	return [origin,rotation]

def transform_pt(meas_basis,fixture_basis,meas_point):
	if meas_point.size!=3:
		return
	point = np.copy(meas_point)
	#print point.shape
	point = point-meas_basis[0]
	point = meas_basis[1].dot(point)
	point = point.dot(fixture_basis[1])
	point = point+fixture_basis[0]
	return point

def transform_pts(meas_basis,fixture_basis,meas_points):
	points = np.copy(meas_points)
	for i in range(0,points.shape[0]):
		points[i,:] = transform_pt(meas_basis,fixture_basis,points[i,:])
	return points

def transform_vec(meas_basis,fixture_basis,meas_vec):
	if meas_vec.size!=3:
		return
	vec = np.copy(meas_vec)
	vec = meas_basis[1].dot(vec)
	vec = vec.dot(fixture_basis[1])
	#vec = fixture_basis[1].dot(vec)
	return vec

def transform_vecs(meas_basis,fixture_basis,meas_vecs):
	vecs = np.copy(meas_vecs)
	for i in range(0,vecs.shape[0]):
		vecs[i,:] = transform_vec(meas_basis,fixture_basis,vecs[i,:])
	return vecs

def transform_basis(meas_basis,fixture_basis,basis_meas):
	origin = transform_pt(meas_basis,fixture_basis,basis_meas[0])
	rotation = transform_vecs(meas_basis,fixture_basis,basis_meas[1])
	#rotation = transform_vec(meas_basis,fixture_basis,basis_meas[1].T).T
	#rotation = np.copy(basis_meas[1].T)
	#for i in range(0,rotation.shape[1]):
	#	rotation[:,i] = transform_vec(meas_basis,fixture_basis,rotation[:,i])
	return [origin,rotation]

def reverse_basis(input_basis):
	#origin = transform_pt(meas_basis,fixture_basis,np.array([0.0, 0.0, 0.0]))
	basis_vecs = input_basis[1].T
	#print basis_vecs
	#print -1.0*input_basis[0]
	basis_origin = input_basis[1].dot(-1.0*input_basis[0])
	#print basis_origin
	return [basis_origin,basis_vecs]

def transform_plane(meas_basis,fixture_basis,meas_plane):
	point = meas_plane[0]*meas_plane[1]
	transformed_point = transform_pt(meas_basis,fixture_basis,point)
	vec = transform_vec(meas_basis,fixture_basis,meas_plane[0])
	distance = transformed_point.dot(vec)
	return [vec, distance]

def make_sensor_basis(sensor):
	rotation = make_rotation(sensor[3],sensor[4],sensor[5])
	origin = rotation.T.dot(-sensor[0:3])
	curvature = np.copy(sensor[6:9])
	return [origin,rotation,curvature]

def find_sensor_intersection(sensor,module_basis,edge_y,x):
	sensor_x = 0;
	sensor_y = edge_y;
	sensor_z = sensor[2][0]*sensor_x*sensor_x + sensor[2][1]*sensor_y*sensor_y + sensor[2][2]*sensor_x*sensor_y

	global_pt = transform_pt(sensor,module_basis,np.array([sensor_x,sensor_y,sensor_z]))
	return global_pt

def find_sensor_residuals(module_basis,sensor,points):
	npoints = points.shape[0]
	#print npoints
	#print points
	transformed_points = transform_pts(module_basis,sensor,points)
	residuals = np.empty([npoints,2])
	for i in range(0,npoints):
		sensor_x = transformed_points[i,0]
		sensor_y = transformed_points[i,1]
		residuals[i,0] = sensor_y - math.copysign(20.17,sensor_y)
		residuals[i,1] = transformed_points[i,2]- ( sensor[2][0]*sensor_x*sensor_x + sensor[2][1]*sensor_y*sensor_y + sensor[2][2]*sensor_x*sensor_y)

	#print transformed_points
	residuals_surveyframe = np.zeros([npoints,3])
	residuals_surveyframe[:,1] = residuals[:,0]
	residuals_surveyframe[:,2] = residuals[:,1]
	#print residuals
	#temp = np.array([[0,0,1]])
	#print temp
	#print transform_pts(sensor,module_basis,temp)
	#print residuals_surveyframe
	residuals_surveyframe = transform_vecs(sensor,module_basis,residuals_surveyframe)
	#print residuals_surveyframe
	print residuals_surveyframe[:,1:3].T
	return residuals

#def sensor_plane_intersection(sensor_basis,sensor_basis,plane):
#	return 

