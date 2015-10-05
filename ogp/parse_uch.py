#!/usr/bin/env python
import sys,math
import numpy as np
import scipy.optimize
import ogp_util

np.set_printoptions(precision=3)

if (len(sys.argv)!=5):
	print sys.argv[0]+' t13.PRT b13.PRT t46.PRT b46.PRT'
	sys.exit()

t13dict = ogp_util.getsteps(sys.argv[1])
b13dict = ogp_util.getsteps(sys.argv[2])
t46dict = ogp_util.getsteps(sys.argv[3])
b46dict = ogp_util.getsteps(sys.argv[4])

def ball_basis(dict, upst, downst):

	holeupst = np.empty([3])
	holedownst = np.empty([3])
	slotupst = np.empty([3])
	slotdownst = np.empty([3])
	holeupst[0]=dict[upst+' hole ball reference']['X Location'][1]
	holeupst[1]=dict[upst+' hole ball reference']['Y Location'][1]
	holeupst[2]=dict[upst+' hole ball reference']['Z Location'][1]
	holedownst[0]=dict[downst+' hole ball reference']['X Location'][1]
	holedownst[1]=dict[downst+' hole ball reference']['Y Location'][1]
	holedownst[2]=dict[downst+' hole ball reference']['Z Location'][1]
	slotupst[0]=dict[upst+' slot ball reference']['X Location'][1]
	slotupst[1]=dict[upst+' slot ball reference']['Y Location'][1]
	slotupst[2]=dict[upst+' slot ball reference']['Z Location'][1]
	slotdownst[0]=dict[downst+' slot ball reference']['X Location'][1]
	slotdownst[1]=dict[downst+' slot ball reference']['Y Location'][1]
	slotdownst[2]=dict[downst+' slot ball reference']['Z Location'][1]

	basis = ogp_util.make_uch_basis(holeupst,holedownst,slotupst,slotdownst)

	return basis

print 'top 1-3:'
t13_basis = ball_basis(t13dict,'L1','L3')
#print t13_basis
t13_global = ogp_util.t13_global()
#print t13_global
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(t13_basis,t13_global,t13_basis[2]),t13_global[2],ogp_util.transform_pt(t13_basis,t13_global,t13_basis[2])-t13_global[2])
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(t13_basis,t13_global,t13_basis[3]),t13_global[3],ogp_util.transform_pt(t13_basis,t13_global,t13_basis[3])-t13_global[3])
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(t13_basis,t13_global,t13_basis[4]),t13_global[4],ogp_util.transform_pt(t13_basis,t13_global,t13_basis[4])-t13_global[4])
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(t13_basis,t13_global,t13_basis[5]),t13_global[5],ogp_util.transform_pt(t13_basis,t13_global,t13_basis[5])-t13_global[5])

print 'bottom 1-3:'
b13_basis = ball_basis(b13dict,'L1','L3')
#print b13_basis
b13_global = ogp_util.b13_global()
#print b13_global
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(b13_basis,b13_global,b13_basis[2]),b13_global[2],ogp_util.transform_pt(b13_basis,b13_global,b13_basis[2])-b13_global[2])
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(b13_basis,b13_global,b13_basis[3]),b13_global[3],ogp_util.transform_pt(b13_basis,b13_global,b13_basis[3])-b13_global[3])
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(b13_basis,b13_global,b13_basis[4]),b13_global[4],ogp_util.transform_pt(b13_basis,b13_global,b13_basis[4])-b13_global[4])
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(b13_basis,b13_global,b13_basis[5]),b13_global[5],ogp_util.transform_pt(b13_basis,b13_global,b13_basis[5])-b13_global[5])

print 'top 4-6:'
t46_basis = ball_basis(t46dict,'L4','L6')
#print b13_basis
t46_global = ogp_util.t46_global()
#print b13_global
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(t46_basis,t46_global,t46_basis[2]),t46_global[2],ogp_util.transform_pt(t46_basis,t46_global,t46_basis[2])-t46_global[2])
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(t46_basis,t46_global,t46_basis[3]),t46_global[3],ogp_util.transform_pt(t46_basis,t46_global,t46_basis[3])-t46_global[3])
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(t46_basis,t46_global,t46_basis[4]),t46_global[4],ogp_util.transform_pt(t46_basis,t46_global,t46_basis[4])-t46_global[4])
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(t46_basis,t46_global,t46_basis[5]),t46_global[5],ogp_util.transform_pt(t46_basis,t46_global,t46_basis[5])-t46_global[5])

print 'bottom 4-6:'
b46_basis = ball_basis(b46dict,'L4','L6')
#print b13_basis
b46_global = ogp_util.b46_global()
#print b13_global
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(b46_basis,b46_global,b46_basis[2]),b46_global[2],ogp_util.transform_pt(b46_basis,b46_global,b46_basis[2])-b46_global[2])
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(b46_basis,b46_global,b46_basis[3]),b46_global[3],ogp_util.transform_pt(b46_basis,b46_global,b46_basis[3])-b46_global[3])
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(b46_basis,b46_global,b46_basis[4]),b46_global[4],ogp_util.transform_pt(b46_basis,b46_global,b46_basis[4])-b46_global[4])
print 'ball {0}, nominal {1}, residual {2}'.format(ogp_util.transform_pt(b46_basis,b46_global,b46_basis[5]),b46_global[5],ogp_util.transform_pt(b46_basis,b46_global,b46_basis[5])-b46_global[5])


def axial_residuals(dict, layer, nom_height, ball_basis, global_basis):

	az = math.radians(dict[layer]['XY Angle'][1])
	el = math.radians(dict[layer]['Elevation'][1])
	vec = np.array([math.cos(az)*math.cos(el),
		math.sin(az)*math.cos(el),
		math.sin(el)])
	pos = np.array([dict[layer]['X Location'][1],
		dict[layer]['Y Location'][1],
		dict[layer]['Z Location'][1]])
	vec = ogp_util.transform_vec(ball_basis,global_basis,vec)
	pos = ogp_util.transform_pt(ball_basis,global_basis,pos)
	angle = math.atan2(vec[1],vec[0])
	#height = pos[1] - pos[0]*(vec[1]/vec[0])
	height = pos[1]
	if (abs(angle)>math.pi/2):
		angle = angle-math.copysign(math.pi,angle)
	#print pos
	#print vec
	#print angle
	height_residual = height-nom_height
	print 'Y residual {0}, angle residual {1}'.format(height_residual,angle)
	return np.array([height_residual,angle])


#def f_min(p):
#	global_basis = ogp_util.make_parametrized_basis(0.0,p[0],0.0,p[1],0.0,p[2])
#	heights = np.array([])
#	avg_angle = 0
#	[height, angle] = axial_residuals(t13dict, 'L1 axial', 0.5, t13_basis, global_basis)
#	heights = np.append(heights,height)
#	avg_angle = avg_angle+angle
#	[height, angle] = axial_residuals(t13dict, 'L2 axial', 2.0, t13_basis, global_basis)
#	heights = np.append(heights,height)
#	avg_angle = avg_angle+angle
#	[height, angle] = axial_residuals(t13dict, 'L3 axial', 3.5, t13_basis, global_basis)
#	heights = np.append(heights,height)
#	avg_angle = avg_angle+angle
#	return np.append(heights,avg_angle*1e20)
#p0 = [8.5, -math.pi/2, 0.0]
#sol = scipy.optimize.leastsq(f_min,p0,args=())[0]
#print sol

axial_residuals(t13dict, 'L1 axial', 0.5, t13_basis, t13_global)
axial_residuals(t13dict, 'L2 axial', 2.0, t13_basis, t13_global)
axial_residuals(t13dict, 'L3 axial', 3.5, t13_basis, t13_global)
#global_basis = ogp_util.make_parametrized_basis(.0,8.5,0.0,0.0,0.0,0.0)
#print ogp_util.transform_pt(t13_basis,global_basis,t13_basis[2])
#print ogp_util.transform_pt(t13_basis,global_basis,t13_basis[3])
#print ogp_util.transform_pt(t13_basis,global_basis,t13_basis[4])
#print ogp_util.transform_pt(t13_basis,global_basis,t13_basis[5])


axial_residuals(b13dict, 'L1 axial', -0.5, b13_basis, b13_global)
axial_residuals(b13dict, 'L2 axial', -2.0, b13_basis, b13_global)
axial_residuals(b13dict, 'L3 axial', -3.5, b13_basis, b13_global)

axial_residuals(t46dict, 'L4 axial block', 6.5, t46_basis, t46_global)
axial_residuals(t46dict, 'L4 axial lever', 6.5, t46_basis, t46_global)
axial_residuals(t46dict, 'L5 axial block', 9.5, t46_basis, t46_global)
axial_residuals(t46dict, 'L5 axial lever', 9.5, t46_basis, t46_global)
axial_residuals(t46dict, 'L6 axial block', 12.5, t46_basis, t46_global)
axial_residuals(t46dict, 'L6 axial lever', 12.5, t46_basis, t46_global)

axial_residuals(b46dict, 'L4 axial block', -6.5, b46_basis, b46_global)
axial_residuals(b46dict, 'L4 axial lever', -6.5, b46_basis, b46_global)
axial_residuals(b46dict, 'L5 axial block', -9.5, b46_basis, b46_global)
axial_residuals(b46dict, 'L5 axial lever', -9.5, b46_basis, b46_global)
axial_residuals(b46dict, 'L6 axial block', -12.5, b46_basis, b46_global)
axial_residuals(b46dict, 'L6 axial lever', -12.5, b46_basis, b46_global)

