#!/usr/bin/env python
import sys,math
import numpy as np
import scipy.optimize
import ogp_util

np.set_printoptions(precision=4)


if (len(sys.argv)!=2):
	print sys.argv[0]+' box_survey_data.csv'
	sys.exit()

stepdict = ogp_util.getsteps_box_survey(sys.argv[1])

print 'upstream center (box origin): '+str((stepdict['TBA']+(stepdict['TBB']+stepdict['TBC'])/2)/2)
print 'downstream center (box +Z):'+str((stepdict['TB1']+(stepdict['TB2']+stepdict['TB3'])/2)/2)
print 'upstream vector (box +X):'+str((stepdict['TBC']-stepdict['TBB'])/2)
print 'downstream vector (box twist):'+str((stepdict['TB3']-stepdict['TB2'])/2)

b13_nom = ogp_util.get_uchbasis(stepdict,'M_BOT13','')
t13_nom = ogp_util.get_uchbasis(stepdict,'M_TOP13','')
b46_nom = ogp_util.get_uchbasis(stepdict,'M_BOT46','')
t46_nom = ogp_util.get_uchbasis(stepdict,'M_TOP46','_CHK2')
#print b13_nom
#print t13_nom
print "B46 ball basis in box fiducial frame:"
print b46_nom
print "T46 ball basis in box fiducial frame:"
print t46_nom


b13_positions = []
b13_positions.append(ogp_util.get_uchbasis(stepdict,'M_BOT13','_-1'))
b13_positions.append(ogp_util.get_uchbasis(stepdict,'M_BOT13','_-2'))
b13_positions.append(ogp_util.get_uchbasis(stepdict,'M_BOT13','_-3'))
b13_positions.append(ogp_util.get_uchbasis(stepdict,'M_BOT13','_-5'))
b13_positions.append(ogp_util.get_uchbasis(stepdict,'M_BOT13','_END'))
b13_positions.append(ogp_util.get_uchbasis(stepdict,'M_BOT13','_0'))

t13_positions = []
t13_positions.append(ogp_util.get_uchbasis(stepdict,'M_TOP13','_+1'))
t13_positions.append(ogp_util.get_uchbasis(stepdict,'M_TOP13','_+2'))
t13_positions.append(ogp_util.get_uchbasis(stepdict,'M_TOP13','_+3'))
t13_positions.append(ogp_util.get_uchbasis(stepdict,'M_TOP13','_+5'))
t13_positions.append(ogp_util.get_uchbasis(stepdict,'M_TOP13','_END'))
t13_positions.append(ogp_util.get_uchbasis(stepdict,'M_TOP13','_0'))
t13_positions.append(ogp_util.get_uchbasis(stepdict,'M_TOP13','_END2'))
t13_positions.append(ogp_util.get_uchbasis(stepdict,'M_TOP13','_02'))
#print b13_positions

null_basis = ogp_util.make_basis(np.array([0.0, 0.0, 0.0]), #ball basis in ball frame
		np.array([1.0, 0.0, 0.0]),
		np.array([0.0, 1.0, 0.0]))

p0 = [-50.0, 300.0,0.0,300.0]

def get_axis(p):
	#return np.array([[0.0,p[0],p[1]],
	#	[50.0,p[2],p[3]]])
	return np.array([[0.0,p[0],p[1]],
		[50.0,p[0],p[1]]])
def get_basis(p):
	axis = get_axis(p)
	diaball = axis[0]+np.array([0.0,1.0,0.0])
	return ogp_util.make_basis(axis[0],axis[1],diaball)

def transform_axis(null_basis,uch_basis,p):
	axis = get_axis(p)
	return ogp_util.transform_pts(null_basis,uch_basis,axis)

def fitfunc(p,positions,nominal):
	axis_nom = transform_axis(null_basis,nominal,p)
	residuals = np.empty(0)
	for position in positions:
		#print position
		axis = transform_axis(null_basis,position,p)
		residuals = np.append(residuals,(axis-axis_nom)[:,1:3])
		#print axis-axis_nom
		#print residuals
	return residuals
#print fitfunc(p0,b13_positions,b13_nom)

sol = scipy.optimize.leastsq(fitfunc,p0,args=(b13_positions,b13_nom))[0]
#print sol
print "B13 ball basis in pivot frame:"
print ogp_util.reverse_basis(get_basis(sol))
#print transform_axis(null_basis,b13_nom,sol)
print "B13 pivot basis in box frame:"
print ogp_util.transform_basis(null_basis,b13_nom,get_basis(sol))
#print fitfunc(sol,b13_positions,b13_nom).reshape(-1,4)

sol = scipy.optimize.leastsq(fitfunc,p0,args=(t13_positions,t13_nom))[0]
#print sol
print "T13 ball basis in pivot frame:"
print ogp_util.reverse_basis(get_basis(sol))
#print transform_axis(null_basis,t13_nom,sol)
print "T13 pivot basis in box frame:"
print ogp_util.transform_basis(null_basis,t13_nom,get_basis(sol))
#print fitfunc(sol,t13_positions,t13_nom).reshape(-1,4)


# tweaking/repeatability tests of b46
#print ogp_util.get_uchbasis(stepdict,'M_TOP46','')
#print ogp_util.get_uchbasis(stepdict,'M_TOP46','_CHK')
#print ogp_util.get_uchbasis(stepdict,'M_TOP46','_MOVE')
#print ogp_util.get_uchbasis(stepdict,'M_TOP46','_ADJ')
#print ogp_util.get_uchbasis(stepdict,'M_TOP46','_MOVE2')
#print ogp_util.get_uchbasis(stepdict,'M_TOP46','_CHK2')
