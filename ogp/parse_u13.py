#!/usr/bin/env python
import sys,math
import numpy as np
import scipy.optimize
import ogp_util

np.set_printoptions(precision=4)

if (len(sys.argv)!=11):
	print sys.argv[0]+' (top|bottom) noballs.PRT balls.PRT modules.PRT L1front.PRT L1back.PRT L2front.PRT L2back.PRT L3front.PRT L3back.PRT '
	sys.exit()

null_basis = ogp_util.make_basis(np.array([0.0, 0.0, 0.0]), #ball basis in ball frame
		np.array([1.0, 0.0, 0.0]),
		np.array([0.0, 1.0, 0.0]))

side = sys.argv[1]
if (side=='top'):
	dz = (11.973-0.5-0.375-0.3125)*25.4
	p0 = [-47.0,8.5,96.0,0.0,0.0,0.0]
	l2_basis = ogp_util.make_pin_basis(np.array([0.0, 0.0, 0.0]), #pin basis in U-channel survey frame
			np.array([1.0, 0.0, 0.0]),
			np.array([0.0, 0.0, 1.0]))
	modulesurvey_basis = ogp_util.make_basis(np.array([0.0, 0.0, 0.0]), #ball basis in module survey frame
			np.array([1.0, 0.0, 0.0]),
			np.array([0.0, 0.0, -1.0]))
	is_top = 1
elif (side=='bottom'):
	dz = (11.973-0.5-0.375+0.3125)*25.4
	p0 = [-45.0,-8.5,104.0,0.0,0.0,0.0]
	l2_basis = ogp_util.make_pin_basis(np.array([0.0, 0.0, 0.0]),
			np.array([-1.0, 0.0, 0.0]),
			np.array([0.0, 0.0, 1.0]))
	modulesurvey_basis = ogp_util.make_basis(np.array([0.0, 0.0, 0.0]),
			np.array([-1.0, 0.0, 0.0]),
			np.array([0.0, 0.0, 1.0]))
	is_top = -1
else:
	print 'Invalid side'
	sys.exit()

l1s = ogp_util.make_sensor_basis(ogp_util.get_l123_sensor(sys.argv[5],'front')) #pin basis in sensor frame
#print l1s
l1a = ogp_util.make_sensor_basis(ogp_util.get_l123_sensor(sys.argv[6],'back'))
#print l1a
l2s = ogp_util.make_sensor_basis(ogp_util.get_l123_sensor(sys.argv[7],'front'))
l2a = ogp_util.make_sensor_basis(ogp_util.get_l123_sensor(sys.argv[8],'back'))
l3s = ogp_util.make_sensor_basis(ogp_util.get_l123_sensor(sys.argv[9],'front'))
l3a = ogp_util.make_sensor_basis(ogp_util.get_l123_sensor(sys.argv[10],'back'))

noballsdict = ogp_util.getsteps(sys.argv[2])
ballsdict = ogp_util.getsteps(sys.argv[3])
modulesdict = ogp_util.getsteps(sys.argv[4])




holeupst = np.empty([3])
holedownst = np.empty([3])
slotupst = np.empty([3])
slotdownst = np.empty([3])
holeupst[0]=ballsdict['L1 hole ball']['X Location'][1]
holeupst[1]=ballsdict['L1 hole ball']['Y Location'][1]
holeupst[2]=ballsdict['L1 hole ball']['Z Location'][1]
holedownst[0]=ballsdict['L3 hole ball']['X Location'][1]
holedownst[1]=ballsdict['L3 hole ball']['Y Location'][1]
holedownst[2]=ballsdict['L3 hole ball']['Z Location'][1]
slotupst[0]=ballsdict['L1 slot ball']['X Location'][1]
slotupst[1]=ballsdict['L1 slot ball']['Y Location'][1]
slotupst[2]=ballsdict['L1 slot ball']['Z Location'][1]
slotdownst[0]=ballsdict['L3 slot ball']['X Location'][1]
slotdownst[1]=ballsdict['L3 slot ball']['Y Location'][1]
slotdownst[2]=ballsdict['L3 slot ball']['Z Location'][1]

ball_basis = ogp_util.make_uch_basis(holeupst,holedownst,slotupst,slotdownst) #ball basis in U-channel survey frame: ball +Z is survey +Y, ball +Y is survey +/-Z, survey origin is L2 hole pin
#print ball_basis
#print modulesurvey_basis

l1_basis = ogp_util.get_uch_mount(noballsdict,'L1') #pin basis in U-channel survey frame: l1, l2, l3 all the same, pin +Y is survey +Z, pin +Z is survey +/-Y, survey origin is L2 hole pin
#print l1_basis
l3_basis = ogp_util.get_uch_mount(noballsdict,'L3')

l1_basis_trans = ogp_util.transform_basis(ball_basis,modulesurvey_basis,l1_basis) #pin basis in module survey frame: pin +Y is survey +Z, pin +Z is survey +/-Y
#print l1_basis_trans

l2_basis_trans = ogp_util.transform_basis(ball_basis,modulesurvey_basis,l2_basis)
#print l2_basis_trans

l3_basis_trans = ogp_util.transform_basis(ball_basis,modulesurvey_basis,l3_basis)
#print l3_basis_trans


print "residuals between survey of mounted modules, and where we expect sensors to be from individual surveys of modules and U-channel:"
l1s_pts = ogp_util.get_sensor_points(modulesdict,'L1 stereo',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l1_basis_trans,l1s,l1s_pts)

l1a_pts = ogp_util.get_sensor_points(modulesdict,'L1 axial',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l1_basis_trans,l1a,l1a_pts)

l2s_pts = ogp_util.get_sensor_points(modulesdict,'L2 stereo',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l2_basis_trans,l2s,l2s_pts)

l2a_pts = ogp_util.get_sensor_points(modulesdict,'L2 axial',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l2_basis_trans,l2a,l2a_pts)

l3s_pts = ogp_util.get_sensor_points(modulesdict,'L3 stereo',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l3_basis_trans,l3s,l3s_pts)

l3a_pts = ogp_util.get_sensor_points(modulesdict,'L3 axial',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l3_basis_trans,l3a,l3a_pts)


l1_basis_ball = ogp_util.transform_basis(ball_basis,null_basis,l1_basis) #pin basis in ball frame: pin +X is ball +X, pin +Y is ball +/-Y
print "L1 pin basis in U-channel fiducial frame:"
print l1_basis_ball

l2_basis_ball = ogp_util.transform_basis(ball_basis,null_basis,l2_basis) #pin basis in ball frame: pin +X is ball +X, pin +Y is ball +/-Y
print "L2 pin basis in U-channel fiducial frame:"
print l2_basis_ball

l3_basis_ball = ogp_util.transform_basis(ball_basis,null_basis,l3_basis) #pin basis in ball frame: pin +X is ball +X, pin +Y is ball +/-Y
print "L3 pin basis in U-channel fiducial frame:"
print l3_basis_ball

print "L1S sensor basis in pin frame:"
print ogp_util.reverse_basis(l1s)
print "L1A sensor basis in pin frame:"
print ogp_util.reverse_basis(l1a)
print "L2S sensor basis in pin frame:"
print ogp_util.reverse_basis(l2s)
print "L2A sensor basis in pin frame:"
print ogp_util.reverse_basis(l2a)
print "L3S sensor basis in pin frame:"
print ogp_util.reverse_basis(l3s)
print "L3A sensor basis in pin frame:"
print ogp_util.reverse_basis(l3a)
sys.exit()

print l2s #stuff in frame 1
print l2_basis #stuff in frame 2
print ogp_util.transform_basis(l2_basis,l2s,l2_basis) #transform stuff in frame 2 to frame 1
print ogp_util.transform_basis(l2s,l2_basis,l2s) #transform stuff in frame 1 to frame 2

print l2s #stuff in frame 1
print l2a #stuff in frame 2
print ogp_util.transform_basis(l2a,l2s,l2a) #transform stuff in frame 2 to frame 1
print ogp_util.transform_basis(l2s,l2a,l2s) #transform stuff in frame 1 to frame 2


#mat = np.array([[-1.,  0.,  0.],
	#[ 0.,  0.,  -1.],
	#[ 0.,  1., -0.]])
#vec = np.array([[0],[1],[0]])
#print mat
#print vec
#print mat.dot(vec)
#print l1s_uch
#print ogp_util.transform_pt(l1s,l1_basis_trans,np.array([10,0,0]))
#l1a_uch = ogp_util.transform_pt(l1a,l1_basis_trans,np.array([0,0,0]))
#print l1a_uch
#print ogp_util.transform_pt(l1a,l1_basis_trans,np.array([10,0,0]))

#print l2s
#print l2a
#l2s_uch = ogp_util.transform_pt(l2s,l2_basis_trans,np.array([0,0,0]))
#print l2s_uch
#l2a_uch = ogp_util.transform_pt(l2a,l2_basis_trans,np.array([0,0,0]))
#print l2a_uch

#l3s_uch = ogp_util.transform_pt(l3s,l3_basis_trans,np.array([0,0,0]))
#print l3s_uch
#l3a_uch = ogp_util.transform_pt(l3a,l3_basis_trans,np.array([0,0,0]))
#print l3a_uch

null_basis = ogp_util.make_parametrized_basis(0.0, 0.0, 0.0, 0.0, 0.0,0.0)
#axial_origin = np.array([124.083, 9.525, 39.1898])
#stereo_origin = np.array([124.067, 13.776, 39.1898])
#axial_edge = np.array([124.083, 9.525, 39.1898+20.17])
axial_edge = np.array([123.66, 0.375*25.4, 2.337*25.4])
axial_vec = np.array([1.0, 0.0, 0.0])


#dz = -0.3125*25.4
#z = 100
def f_min(p):
	global_basis =ogp_util.make_parametrized_basis(p[0],p[1],p[2],p[3],p[4],p[5])
	#print p
	def residual(mod_basis,uch_basis,y,z,dz):
		edge = ogp_util.transform_pt(ball_basis,uch_basis,ogp_util.transform_pt(null_basis,mod_basis,axial_edge))
		#print edge
		nominal = [math.sin(0.0305)*(z+dz-700)-8.789, y, math.cos(0.0305)*(z+dz)]
		vec = ogp_util.transform_vec(ball_basis,uch_basis,ogp_util.transform_vec(null_basis,mod_basis,axial_vec))
		#print vec
		angle = math.atan2(vec[1],vec[0])
		if (abs(angle)>math.pi/2):
			angle = angle-math.copysign(math.pi,angle)
		#print edge-nominal
		#print angle
		return np.append(edge-nominal,angle*1000)
	return np.concatenate((residual(l1_basis,global_basis,is_top*0.5,100,dz),
		residual(l2_basis,global_basis,is_top*2.0,200,dz),
		residual(l3_basis,global_basis,is_top*3.5,300,dz)))

sol = scipy.optimize.leastsq(f_min,p0,args=())[0]
#print sol
#print f_min(sol)

global_basis =ogp_util.make_parametrized_basis(sol[0],sol[1],sol[2],sol[3],sol[4],sol[5])
#if (is_top==1):
#	print 'top 1-3:'
#else:
#	print 'bottom 1-3:'
#print ogp_util.transform_pt(ball_basis,global_basis,ball_basis[2])/25.4
#print ogp_util.transform_pt(ball_basis,global_basis,ball_basis[3])/25.4
#print ogp_util.transform_pt(ball_basis,global_basis,ball_basis[4])/25.4
#print ogp_util.transform_pt(ball_basis,global_basis,ball_basis[5])/25.4
