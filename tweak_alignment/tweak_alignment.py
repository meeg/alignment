#!/usr/bin/env python
import numpy as np
import math
import scipy.optimize


def fitfunc(p, initial_du, dudx, dudy):
    z = np.array([100,100,
        200,200,
        300,300,
        500,500,500,500,
        700,700,700,700,
        900,900,900,900])
    dOmega = p[0]
    dLambda = p[1]
    dPhi = p[2]
    dD0 = p[3]
    dZ0 = p[4]
    du = dZ0*dudy + dLambda*z*dudy + dD0*dudx + dPhi*z*dudx - dOmega*0.5*z*z*dudx
    #print du
    return du+initial_du

def fitfunc_fixomega(p, initial_du, dudx, dudy,dOmega):
    z = np.array([100,100,
        200,200,
        300,300,
        500,500,500,500,
        700,700,700,700,
        900,900,900,900])
    dLambda = p[0]
    dPhi = p[1]
    dD0 = p[2]
    dZ0 = p[3]
    du = dZ0*dudy + dLambda*z*dudy + dD0*dudx + dPhi*z*dudx - dOmega*0.5*z*z*dudx
    #print du
    return du+initial_du

def fitfunc_fixomega_d0_z0(p, initial_du, dudx, dudy, dOmega, dD0, dZ0):
    z = np.array([100,100,
        200,200,
        300,300,
        500,500,500,500,
        700,700,700,700,
        900,900,900,900])
    dLambda = p[0]
    dPhi = p[1]
    du = dZ0*dudy + dLambda*z*dudy + dD0*dudx + dPhi*z*dudx - dOmega*0.5*z*z*dudx
    #print du
    return du+initial_du

def fitfunc_constraintarget(p, initial_du1, dudx1, dudy1, initial_du2, dudx2, dudy2, dD01, dD02, dZ01, dZ02, dOmega1, dOmega2):
    z = np.array([100,100,
        200,200,
        300,300,
        500,500,500,500,
        700,700,700,700,
        900,900,900,900])
    dLambda1 = p[0]
    dPhi1 = p[1]
    dLambda2 = p[2]
    dPhi2 = p[3]
    dD0 = p[4]
    dZ0 = p[5]
    du1 = (dZ0+(dZ01-dZ02)/2)*dudy1 + dLambda1*z*dudy1 + (dD0+(dD01-dD02)/2)*dudx1 + dPhi1*z*dudx1 - dOmega1*0.5*z*z*dudx1
    du2 = (dZ0+(dZ02-dZ01)/2)*dudy2 + dLambda2*z*dudy2 + (dD0+(dD02-dD01)/2)*dudx2 + dPhi2*z*dudx2 - dOmega2*0.5*z*z*dudx2
    #print du
    return np.concatenate((du1+initial_du1,du2+initial_du2))

def print_du(du):
	print du
	print np.sqrt(np.mean(np.square(du)))

def print_tweak(tweak,half):
	print tweak
	print '<constant name="dOmega_'+half+'" value="'+str(tweak[0])+'" />'
	print '<constant name="dLambda_'+half+'" value="'+str(tweak[1])+'" />'
	print '<constant name="dPhi_'+half+'" value="'+str(tweak[2])+'" />'
	print '<constant name="dD0_'+half+'" value="'+str(tweak[3])+'" />'
	print '<constant name="dZ0_'+half+'" value="'+str(tweak[4])+'" />'


dudy_top = np.array([1,-1/math.cos(0.1),
    1,-1/math.cos(0.1),
    1,-1/math.cos(0.1),
    1,-1/math.cos(0.05),-1,1/math.cos(0.05),
    1,-1/math.cos(0.05),-1,1/math.cos(0.05),
    1,-1/math.cos(0.05),-1,1/math.cos(0.05)])
dudx_top = np.array([0,-math.sin(0.1),
    0,-math.sin(0.1),
    0,-math.sin(0.1),
    0,-math.sin(0.05),0,math.sin(0.05),
    0,-math.sin(0.05),0,math.sin(0.05),
    0,-math.sin(0.05),0,math.sin(0.05)])
dudy_bot = np.array([1/math.cos(0.1),-1,
    1/math.cos(0.1),-1,
    1/math.cos(0.1),-1,
    1/math.cos(0.05),-1,-1/math.cos(0.05),1,
    1/math.cos(0.05),-1,-1/math.cos(0.05),1,
    1/math.cos(0.05),-1,-1/math.cos(0.05),1])
dudx_bot = np.array([-math.sin(0.1),0,
    -math.sin(0.1),0,
    -math.sin(0.1),0,
    -math.sin(0.05),0,math.sin(0.05),0,
    -math.sin(0.05),0,math.sin(0.05),0,
    -math.sin(0.05),0,math.sin(0.05),0])

step5_top = np.array([0.142557,
    -0.108300,
    0.129410,
    -0.071933,
    0.111141,
    -0.048276,
    -0.000909,
    0.001063,
    0.009966,
    -0.010293,
    -0.016174,
    0.013722,
    -0.001454,
    -0.000033,
    0.008125,
    -0.012739,
    -0.051231,
    0.060625])

step5_bot = np.array([-0.530096,
    0.035614,
    -0.381392,
    0.020565,
    -0.277673,
    0.023975,
    -0.006136,
    0.007321,
    0.005669,
    0.000505,
    -0.005321,
    0.009127,
    0.014415,
    -0.016824,
    0.015446,
    -0.022457,
    -0.036677,
    0.036675])

params_top = np.array([0.4e-5,
    0.5e-3,
    1.5e-3,
    0.8,
    -0.15])
params_bot = np.array([2.5e-5,
    1.5e-3,
    17e-3,
    -5.7,
    0.15])

print "top"
print "global tweak parameters from external constraints (omega,lambda,phi,d0,z0):"
print_tweak(params_top,'top')
print "millepede constants after external constraints:"
print_du(fitfunc(params_top,step5_top,dudx_top,dudy_top))

sol = scipy.optimize.leastsq(fitfunc,params_top,args=(step5_top,dudx_top,dudy_top))[0]
print "global tweak parameters from survey constraints (omega,lambda,phi,d0,z0):"
print sol
print_tweak(sol,'top')
print "millepede constants after survey constraints:"
print_du(fitfunc(sol,step5_top,dudx_top,dudy_top))

sol = scipy.optimize.leastsq(fitfunc_fixomega,params_top[1:5],args=(step5_top,dudx_top,dudy_top,params_top[0]))[0]
print "global tweak parameters from mixed (external omega + survey) constraints (omega,lambda,phi,d0,z0):"
tweak = np.concatenate((np.array([params_top[0]]),sol))
print_tweak(tweak,'top')
print "millepede constants after mixed constraints:"
print_du(fitfunc_fixomega(sol,step5_top,dudx_top,dudy_top,params_top[0]))

sol = scipy.optimize.leastsq(fitfunc_fixomega_d0_z0,params_top[1:3],args=(step5_top,dudx_top,dudy_top,params_top[0],params_top[3],params_top[4]))[0]
print "global tweak parameters from mixed (external omega, d0, z0 + survey) constraints (omega,lambda,phi,d0,z0):"
tweak = np.array([params_top[0],sol[0],sol[1],params_top[3],params_top[4]])
print_tweak(tweak,'top')
print "millepede constants after mixed constraints:"
print_du(fitfunc(tweak,step5_top,dudx_top,dudy_top))

print

print "bottom"
print "global tweak parameters from external constraints (omega,lambda,phi,d0,z0):"
print_tweak(params_bot,'bot')
print "millepede constants after external constraints:"
print_du(fitfunc(params_bot,step5_bot,dudx_bot,dudy_bot))

sol = scipy.optimize.leastsq(fitfunc,params_bot,args=(step5_bot,dudx_bot,dudy_bot))[0]
print "global tweak parameters from survey constraints (omega,lambda,phi,d0,z0):"
print_tweak(sol,'bot')
print "millepede constants after survey constraints:"
print_du(fitfunc(sol,step5_bot,dudx_bot,dudy_bot))

sol = scipy.optimize.leastsq(fitfunc_fixomega,params_bot[1:5],args=(step5_bot,dudx_bot,dudy_bot,params_bot[0]))[0]
print "global tweak parameters from mixed (external omega + survey) constraints (omega,lambda,phi,d0,z0):"
tweak = np.concatenate((np.array([params_bot[0]]),sol))
print_tweak(tweak,'bot')
print "millepede constants after mixed constraints:"
print_du(fitfunc_fixomega(sol,step5_bot,dudx_bot,dudy_bot,params_bot[0]))

sol = scipy.optimize.leastsq(fitfunc_fixomega_d0_z0,params_bot[1:3],args=(step5_bot,dudx_bot,dudy_bot,params_bot[0],params_bot[3],params_bot[4]))[0]
print "global tweak parameters from mixed (external omega, d0, z0 + survey) constraints (omega,lambda,phi,d0,z0):"
tweak = np.array([params_bot[0],sol[0],sol[1],params_bot[3],params_bot[4]])
print_tweak(tweak,'bot')
print "millepede constants after mixed constraints:"
print_du(fitfunc(tweak,step5_bot,dudx_bot,dudy_bot))


print
sol = scipy.optimize.leastsq(fitfunc_constraintarget,np.array([params_top[1],params_top[2],params_bot[1],params_bot[2],(params_top[3]+params_bot[3])/2,(params_top[4]+params_bot[4])/2]),args=(step5_top,dudx_top,dudy_top,step5_bot,dudx_bot,dudy_bot,params_top[3],params_bot[3],params_top[4],params_bot[4],params_top[0],params_bot[0]))[0]
#print sol
print "global tweak parameters and millepede constants from mixed2 (external omega + d0, z0 constraint + survey) constraints (omega,lambda,phi,d0,z0):"
print "top:"
top = np.array([params_top[0],sol[0],sol[1],sol[4]+(params_top[3]-params_bot[3])/2,sol[5]+(params_top[4]-params_bot[4])/2])
print_tweak(top,'top')
print_du(fitfunc(top,step5_top,dudx_top,dudy_top))
print "bottom:"
bot = np.array([params_bot[0],sol[2],sol[3],sol[4]+(params_bot[3]-params_top[3])/2,sol[5]+(params_bot[4]-params_top[4])/2])
print_tweak(bot,'bot')
#print fitfunc_constraintarget(sol,step5_top,dudx_top,dudy_top,step5_bot,dudx_bot,dudy_bot,params_top[3],params_bot[3],params_top[4],params_bot[4],params_top[0],params_bot[0])
print_du(fitfunc(bot,step5_bot,dudx_bot,dudy_bot))
