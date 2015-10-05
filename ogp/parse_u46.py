#!/usr/bin/env python
import sys,math
import numpy as np
import scipy.optimize
import ogp_util

np.set_printoptions(precision=4)

if (len(sys.argv)!=11):
	print sys.argv[0]+' (top|bottom) noballs.PRT balls.PRT modules.PRT L1front.PRT L1back.PRT L2front.PRT L2back.PRT L3front.PRT L3back.PRT'
	sys.exit()

null_basis = ogp_util.make_basis(np.array([0.0, 0.0, 0.0]), #ball basis in ball frame
		np.array([1.0, 0.0, 0.0]),
		np.array([0.0, 1.0, 0.0]))

side = sys.argv[1]
if (side=='top'):
	dz = (11.973-0.5-0.375-0.3125)*25.4
	p0 = [-47.0,8.5,96.0,0.0,0.0,0.0]
	l5_basis = ogp_util.make_pin_basis(np.array([0.0, 0.0, 0.0]), #pin basis in U-channel survey frame
			np.array([1.0, 0.0, 0.0]),
			np.array([0.0, 0.0, 1.0]))
	modulesurvey_basis = ogp_util.make_basis(np.array([0.0, 0.0, 0.0]), #ball basis in module survey frame
			np.array([1.0, 0.0, 0.0]),
			np.array([0.0, 0.0, -1.0]))
	is_top = 1
elif (side=='bottom'):
	dz = (11.973-0.5-0.375+0.3125)*25.4
	p0 = [-45.0,-8.5,104.0,0.0,0.0,0.0]
	l5_basis = ogp_util.make_pin_basis(np.array([0.0, 0.0, 0.0]),
			np.array([-1.0, 0.0, 0.0]),
			np.array([0.0, 0.0, 1.0]))
	modulesurvey_basis = ogp_util.make_basis(np.array([0.0, 0.0, 0.0]),
			np.array([-1.0, 0.0, 0.0]),
			np.array([0.0, 0.0, 1.0]))
	is_top = -1
else:
	print 'Invalid side'
	sys.exit()


[l4se, l4sp] = ogp_util.get_l456_sensors(sys.argv[5],'front')
l4se = ogp_util.make_sensor_basis(l4se)
l4sp = ogp_util.make_sensor_basis(l4sp)
[l4ae, l4ap] = ogp_util.get_l456_sensors(sys.argv[6],'back')
l4ae = ogp_util.make_sensor_basis(l4ae)
l4ap = ogp_util.make_sensor_basis(l4ap)
[l5se, l5sp] = ogp_util.get_l456_sensors(sys.argv[7],'front')
l5se = ogp_util.make_sensor_basis(l5se)
l5sp = ogp_util.make_sensor_basis(l5sp)
[l5ae, l5ap] = ogp_util.get_l456_sensors(sys.argv[8],'back')
l5ae = ogp_util.make_sensor_basis(l5ae)
l5ap = ogp_util.make_sensor_basis(l5ap)
[l6se, l6sp] = ogp_util.get_l456_sensors(sys.argv[9],'front')
l6se = ogp_util.make_sensor_basis(l6se)
l6sp = ogp_util.make_sensor_basis(l6sp)
[l6ae, l6ap] = ogp_util.get_l456_sensors(sys.argv[10],'back')
l6ae = ogp_util.make_sensor_basis(l6ae)
l6ap = ogp_util.make_sensor_basis(l6ap)

#print l4se
#print l4sp
#print l4ae
#print l4ap
#print l5se
#print l5sp
#print l5ae
#print l5ap
#print l6se
#print l6sp
#print l6ae
#print l6ap


noballsdict = ogp_util.getsteps(sys.argv[2])
ballsdict = ogp_util.getsteps(sys.argv[3])
modulesdict = ogp_util.getsteps(sys.argv[4])

holeupst = np.empty([3])
holedownst = np.empty([3])
slotupst = np.empty([3])
slotdownst = np.empty([3])
holeupst[0]=ballsdict['L4 hole ball']['X Location'][1]
holeupst[1]=ballsdict['L4 hole ball']['Y Location'][1]
holeupst[2]=ballsdict['L4 hole ball']['Z Location'][1]
holedownst[0]=ballsdict['L6 hole ball']['X Location'][1]
holedownst[1]=ballsdict['L6 hole ball']['Y Location'][1]
holedownst[2]=ballsdict['L6 hole ball']['Z Location'][1]
slotupst[0]=ballsdict['L4 slot ball']['X Location'][1]
slotupst[1]=ballsdict['L4 slot ball']['Y Location'][1]
slotupst[2]=ballsdict['L4 slot ball']['Z Location'][1]
slotdownst[0]=ballsdict['L6 slot ball']['X Location'][1]
slotdownst[1]=ballsdict['L6 slot ball']['Y Location'][1]
slotdownst[2]=ballsdict['L6 slot ball']['Z Location'][1]

ball_basis = ogp_util.make_uch_basis(holeupst,holedownst,slotupst,slotdownst) #ball basis in U-channel survey frame: ball +Z is survey +Y, ball +Y is survey +/-Z, survey origin is L5 hole pin
#print ball_basis
#print modulesurvey_basis

l4_basis = ogp_util.get_uch_mount(noballsdict,'L4') #pin basis in U-channel survey frame: l4, l5, l6 all the same, pin +Y is survey +Z, pin +Z is survey +/-Y, survey origin is L2 hole pin
#print l4_basis
l6_basis = ogp_util.get_uch_mount(noballsdict,'L6')


l4_basis_trans = ogp_util.transform_basis(ball_basis,modulesurvey_basis,l4_basis) #pin basis in module survey frame: pin +Y is survey +Z, pin +Z is survey +/-Y
#print l1_basis_trans

l5_basis_trans = ogp_util.transform_basis(ball_basis,modulesurvey_basis,l5_basis)
#print l2_basis_trans

l6_basis_trans = ogp_util.transform_basis(ball_basis,modulesurvey_basis,l6_basis)
#print l3_basis_trans

print "residuals between survey of mounted modules, and where we expect sensors to be from individual surveys of modules and U-channel:"
l4se_pts = ogp_util.get_sensor_points(modulesdict,'L4 stereo block',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l4_basis_trans,l4se,l4se_pts)

l4sp_pts = ogp_util.get_sensor_points(modulesdict,'L4 stereo lever',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l4_basis_trans,l4sp,l4sp_pts)

l4ae_pts = ogp_util.get_sensor_points(modulesdict,'L4 axial block',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l4_basis_trans,l4ae,l4ae_pts)

l4ap_pts = ogp_util.get_sensor_points(modulesdict,'L4 axial lever',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l4_basis_trans,l4ap,l4ap_pts)

l5se_pts = ogp_util.get_sensor_points(modulesdict,'L5 stereo block',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l5_basis_trans,l5se,l5se_pts)

l5sp_pts = ogp_util.get_sensor_points(modulesdict,'L5 stereo lever',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l5_basis_trans,l5sp,l5sp_pts)

l5ae_pts = ogp_util.get_sensor_points(modulesdict,'L5 axial block',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l5_basis_trans,l5ae,l5ae_pts)

l5ap_pts = ogp_util.get_sensor_points(modulesdict,'L5 axial lever',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l5_basis_trans,l5ap,l5ap_pts)

l6se_pts = ogp_util.get_sensor_points(modulesdict,'L6 stereo block',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l6_basis_trans,l6se,l6se_pts)

l6sp_pts = ogp_util.get_sensor_points(modulesdict,'L6 stereo lever',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l6_basis_trans,l6sp,l6sp_pts)

l6ae_pts = ogp_util.get_sensor_points(modulesdict,'L6 axial block',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l6_basis_trans,l6ae,l6ae_pts)

l6ap_pts = ogp_util.get_sensor_points(modulesdict,'L6 axial lever',3) #sensor edge points in module survey frame
ogp_util.find_sensor_residuals(l6_basis_trans,l6ap,l6ap_pts)


l4_basis_ball = ogp_util.transform_basis(ball_basis,null_basis,l4_basis) #pin basis in ball frame: pin +X is ball +X, pin +Y is ball +/-Y
print "L4 pin basis in U-channel fiducial frame:"
print l4_basis_ball

l5_basis_ball = ogp_util.transform_basis(ball_basis,null_basis,l5_basis) #pin basis in ball frame: pin +X is ball +X, pin +Y is ball +/-Y
print "L5 pin basis in U-channel fiducial frame:"
print l5_basis_ball

l6_basis_ball = ogp_util.transform_basis(ball_basis,null_basis,l6_basis) #pin basis in ball frame: pin +X is ball +X, pin +Y is ball +/-Y
print "L6 pin basis in U-channel fiducial frame:"
print l6_basis_ball


print "L4SE sensor basis in pin frame:"
print ogp_util.reverse_basis(l4se)
print "L4SP sensor basis in pin frame:"
print ogp_util.reverse_basis(l4sp)
print "L4AE sensor basis in pin frame:"
print ogp_util.reverse_basis(l4ae)
print "L4AP sensor basis in pin frame:"
print ogp_util.reverse_basis(l4ap)
print "L5SE sensor basis in pin frame:"
print ogp_util.reverse_basis(l5se)
print "L5SP sensor basis in pin frame:"
print ogp_util.reverse_basis(l5sp)
print "L5AE sensor basis in pin frame:"
print ogp_util.reverse_basis(l5ae)
print "L5AP sensor basis in pin frame:"
print ogp_util.reverse_basis(l5ap)
print "L6SE sensor basis in pin frame:"
print ogp_util.reverse_basis(l6se)
print "L6SP sensor basis in pin frame:"
print ogp_util.reverse_basis(l6sp)
print "L6AE sensor basis in pin frame:"
print ogp_util.reverse_basis(l6ae)
print "L6AP sensor basis in pin frame:"
print ogp_util.reverse_basis(l6ap)

sys.exit()

#print l2_basis

l1holepin = np.empty([3])
l1slotpin = np.empty([3])
l1normal = np.empty([3])
l1holepin[0]=noballsdict['L4 hole pin intersection']['X Location'][1]
l1holepin[1]=noballsdict['L4 hole pin intersection']['Y Location'][1]
l1holepin[2]=noballsdict['L4 hole pin intersection']['Z Location'][1]
l1slotpin[0]=noballsdict['L4 slot pin intersection']['X Location'][1]
l1slotpin[1]=noballsdict['L4 slot pin intersection']['Y Location'][1]
l1slotpin[2]=noballsdict['L4 slot pin intersection']['Z Location'][1]
l1az = math.radians(noballsdict['L4 plane, touch']['XY Angle'][1])
l1el = math.radians(noballsdict['L4 plane, touch']['Elevation'][1])
if (l1el<0):
	l1az = l1az+math.pi
	l1el = -1.0*l1el
l1normal = np.array([math.cos(l1az)*math.cos(l1el),
	math.sin(l1az)*math.cos(l1el),
	math.sin(l1el)])
#print l1normal
l1_basis = ogp_util.make_pin_basis(l1holepin, l1slotpin, l1normal)
#print l1_basis

l3holepin = np.empty([3])
l3slotpin = np.empty([3])
l3normal = np.empty([3])
l3holepin[0]=noballsdict['L6 hole pin intersection']['X Location'][1]
l3holepin[1]=noballsdict['L6 hole pin intersection']['Y Location'][1]
l3holepin[2]=noballsdict['L6 hole pin intersection']['Z Location'][1]
l3slotpin[0]=noballsdict['L6 slot pin intersection']['X Location'][1]
l3slotpin[1]=noballsdict['L6 slot pin intersection']['Y Location'][1]
l3slotpin[2]=noballsdict['L6 slot pin intersection']['Z Location'][1]
l3az = math.radians(noballsdict['L6 plane, touch']['XY Angle'][1])
l3el = math.radians(noballsdict['L6 plane, touch']['Elevation'][1])
if (l3el<0):
	l3az = l3az+math.pi
	l3el = -1.0*l3el
l3normal = np.array([math.cos(l3az)*math.cos(l3el),
	math.sin(l3az)*math.cos(l3el),
	math.sin(l3el)])
#print l3normal
l3_basis = ogp_util.make_pin_basis(l3holepin, l3slotpin, l3normal)
#print l3_basis


null_basis = ogp_util.make_parametrized_basis(0.0, 0.0, 0.0, 0.0, 0.0,0.0)
#axial_origin = np.array([138.125, 9.525, 35.1004])
#stereo_origin = np.array([124.067, 13.776, 39.1898])
axial_edge = np.array([149.225, 0.375*25.4, 35.1004+20.17])
#axial_edge = np.array([138.125, 9.525, 35.1004+20.17])
axial_vec = np.array([1.0, 0.0, 0.0])

#dz = -0.3125*25.4
#z = 100
def f_min(p):
	global_basis =ogp_util.make_parametrized_basis(p[0],p[1],p[2],p[3],p[4],p[5])
	#print p
	def residual(mod_basis,uch_basis,y,z,dz):
		edge = ogp_util.transform_pt(ball_basis,uch_basis,ogp_util.transform_pt(null_basis,mod_basis,axial_edge))
		#print edge
		#nominal = [math.sin(0.0305)*(z+dz-700), y, math.cos(0.0305)*(z+dz-700)+700]
		nominal = [math.sin(0.0305)*(z+dz-700)-8.789, y, math.cos(0.0305)*(z+dz)]
		vec = ogp_util.transform_vec(ball_basis,uch_basis,ogp_util.transform_vec(null_basis,mod_basis,axial_vec))
		#print vec
		angle = math.atan2(vec[1],vec[0])
		if (abs(angle)>math.pi/2):
			angle = angle-math.copysign(math.pi,angle)
		#print edge-nominal
		#print angle
		return np.append(edge-nominal,angle*1000)
	return np.concatenate((residual(l1_basis,global_basis,is_top*6.5,500,dz),
		residual(l2_basis,global_basis,is_top*9.5,700,dz),
		residual(l3_basis,global_basis,is_top*12.5,900,dz)))

sol = scipy.optimize.leastsq(f_min,p0,args=())[0]
print sol
print f_min(sol)

global_basis =ogp_util.make_parametrized_basis(sol[0],sol[1],sol[2],sol[3],sol[4],sol[5])
if (is_top==1):
	print 'top 4-6:'
else:
	print 'bottom 4-6:'
print ogp_util.transform_pt(ball_basis,global_basis,ball_basis[2])/25.4
print ogp_util.transform_pt(ball_basis,global_basis,ball_basis[3])/25.4
print ogp_util.transform_pt(ball_basis,global_basis,ball_basis[4])/25.4
print ogp_util.transform_pt(ball_basis,global_basis,ball_basis[5])/25.4
