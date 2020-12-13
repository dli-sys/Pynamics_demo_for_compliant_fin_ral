# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 21:30:03 2020

@author: dongting
"""

import pynamics
from pynamics.frame import Frame
from pynamics.variable_types import Differentiable,Constant
from pynamics.system import System
from pynamics.body import Body
from pynamics.dyadic import Dyadic
from pynamics.output import Output,PointsOutput
from pynamics.particle import Particle
import pynamics.integration
import logging
import sympy
#import sympy
import numpy
import matplotlib.pyplot as plt
from math import pi
from scipy import optimize
from sympy import sin
import pynamics.tanh as tanh
import time
from fit_qs import exp_fit
import fit_qs



def Cal_system(initial_states,drag_direction,tinitial,tstep,tfinal,fit_vel,f1,f2):
    
    g_k,g_b_damping,g_b_damping1=  [0.30867935, 1.42946955, 1.08464536]
    system = System()
    pynamics.set_system(__name__,system)
    
    global_q = True
        
    lO = Constant(7/1000,'lO',system)
    lA = Constant(33/1000,'lA',system)
    lB = Constant(33/1000,'lB',system)
    lC = Constant(33/1000,'lC',system)
    
    mO = Constant(10/1000,'mA',system)
    mA = Constant(2.89/1000,'mA',system)
    mB = Constant(2.89/1000,'mB',system)
    mC = Constant(2.89/1000,'mC',system)
    k = Constant(g_k,'k',system)
    k1 = Constant(0.4,'k1',system)
    
    friction_perp = Constant(f1,'f_perp',system)
    friction_par = Constant(f2,'f_par',system)
    b_damping = Constant(g_b_damping,'b_damping',system)
    b_damping1 = Constant(g_b_damping1,'b_damping1',system)
       
    preload0 = Constant(0*pi/180,'preload0',system)
    preload1 = Constant(0*pi/180,'preload1',system)
    preload2 = Constant(0*pi/180,'preload2',system)
    preload3 = Constant(0*pi/180,'preload3',system)
    
    Ixx_O = Constant(1,'Ixx_O',system)
    Iyy_O = Constant(1,'Iyy_O',system)
    Izz_O = Constant(1,'Izz_O',system)
    Ixx_A = Constant(1,'Ixx_A',system)
    Iyy_A = Constant(1,'Iyy_A',system)
    Izz_A = Constant(1,'Izz_A',system)
    Ixx_B = Constant(1,'Ixx_B',system)
    Iyy_B = Constant(1,'Iyy_B',system)
    Izz_B = Constant(1,'Izz_B',system)
    Ixx_C = Constant(1,'Ixx_C',system)
    Iyy_C = Constant(1,'Iyy_C',system)
    Izz_C = Constant(1,'Izz_C',system)
    
    y,y_d,y_dd = Differentiable('y',system)
    qO,qO_d,qO_dd = Differentiable('qO',system)
    qA,qA_d,qA_dd = Differentiable('qA',system)
    qB,qB_d,qB_dd = Differentiable('qB',system)
    qC,qC_d,qC_dd = Differentiable('qC',system)
     
    fit_states = initial_states
    
    initialvalues = {}
    initialvalues[y]    = fit_states[0]
    initialvalues[y_d]  = fit_states[5]
    initialvalues[qO]   = 0
    initialvalues[qO_d] = 0
    initialvalues[qA]   = fit_states[2]
    initialvalues[qA_d] = fit_states[7]
    initialvalues[qB]   = fit_states[3]
    initialvalues[qB_d] = fit_states[8]
    initialvalues[qC]   = fit_states[4]
    initialvalues[qC_d] = fit_states[9]
   
    statevariables = system.get_state_variables()
    ini = [initialvalues[item] for item in statevariables]
    
    N = Frame('N')
    O = Frame('O')
    A = Frame('A')
    B = Frame('B')
    C = Frame('C')
      
    system.set_newtonian(N)
    if not global_q:
        O.rotate_fixed_axis_directed(N,[0,0,1],qO,system)
        A.rotate_fixed_axis_directed(O,[0,0,1],qA,system)
        B.rotate_fixed_axis_directed(A,[0,0,1],qB,system)
        C.rotate_fixed_axis_directed(B,[0,0,1],qC,system)
    else:
        O.rotate_fixed_axis_directed(N,[0,0,1],qO,system)
        A.rotate_fixed_axis_directed(N,[0,0,1],qA,system)
        B.rotate_fixed_axis_directed(N,[0,0,1],qB,system)
        C.rotate_fixed_axis_directed(N,[0,0,1],qC,system)
    
    pNO=    0*N.x + y*N.y
    pOA=    lO*N.x + y*N.y
    pAB=    pOA+lA*A.x
    pBC =   pAB + lB*B.x
    pCtip = pBC + lC*C.x
    
    pOcm= pNO +lO/2*N.x
    pAcm= pOA+lA/2*A.x
    pBcm= pAB+lB/2*B.x
    pCcm= pBC+lC/2*C.x
    
    wNO = N.getw_(O)
    wOA = N.getw_(A)
    wAB = A.getw_(B)
    wBC = B.getw_(C)
    
    IO = Dyadic.build(O,Ixx_O,Iyy_O,Izz_O)
    IA = Dyadic.build(A,Ixx_A,Iyy_A,Izz_A)
    IB = Dyadic.build(B,Ixx_B,Iyy_B,Izz_B)
    IC = Dyadic.build(C,Ixx_C,Iyy_C,Izz_C)
    
    
    BodyO = Body('BodyO',O,pOcm,mO,IO,system)
    BodyA = Body('BodyA',A,pAcm,mA,IA,system)
    BodyB = Body('BodyB',B,pBcm,mB,IB,system)
    BodyC = Body('BodyC',C,pCcm,mC,IC,system)
    
    vAcm = pAcm.time_derivative()
    vBcm = pBcm.time_derivative()
    vCcm = pCcm.time_derivative()
    
    system.add_spring_force1(k1+10000*(qA+abs(qA)),(qA-qO-preload1)*N.z,wOA) 
    system.add_spring_force1(k+10000*(qB+abs(qB)),(qB-qA-preload2)*N.z,wAB)
    system.add_spring_force1(k+10000*(qC+abs(qC)),(qC-qB-preload3)*N.z,wBC)

    
    #new Method use nJoint
    nvAcm = 1/vAcm.length()*vAcm
    nvBcm = 1/vBcm.length()*vBcm
    nvCcm = 1/vCcm.length()*vCcm
    
    vSoil = drag_direction*1*N.y
    nSoil = 1/vSoil.length()*vSoil
    
    if fit_vel ==0:
        vSoil = 1*1*N.y
        nSoil = 1/vSoil.length()*vSoil
        
        faperp  = friction_perp*nSoil.dot(A.y)*A.y
        fapar   = friction_par*nSoil.dot(A.x)*A.x
        system.addforce(-(faperp+fapar),vAcm) 
        
        fbperp = friction_perp*nSoil.dot(B.y)*B.y
        fbpar  = friction_par*nSoil.dot(B.x)*B.x
        system.addforce(-(fbperp+fbpar),vBcm) 
        
        fcperp = friction_perp*nSoil.dot(C.y)*C.y
        fcpar  = friction_par*nSoil.dot(C.x)*C.x
        system.addforce(-(fcperp+fcpar),vCcm)  
    else:
        faperp  = friction_perp*nvAcm.dot(A.y)*A.y
        fapar   = friction_par*nvAcm.dot(A.x)*A.x
        system.addforce(-(faperp+fapar),vAcm) 
        
        fbperp = friction_perp*nvBcm.dot(B.y)*B.y
        fbpar  = friction_par*nvBcm.dot(B.x)*B.x
        system.addforce(-(fbperp+fbpar),vBcm) 
        
        fcperp = friction_perp*nvCcm.dot(C.y)*C.y
        fcpar  = friction_par*nvCcm.dot(C.x)*C.x
        system.addforce(-(fcperp+fcpar),vCcm) 
    
    system.addforce(-b_damping1*wOA,wOA)
    system.addforce(-b_damping*wAB,wAB)
    system.addforce(-b_damping*wBC,wBC)
    
    eq = []
    eq_d=[(system.derivative(item)) for item in eq]
    
    eq_d.append(y_d-fit_vel)
    eq_dd=[(system.derivative(item)) for item in eq_d]
        
    f,ma = system.getdynamics()
    func1 = system.state_space_post_invert(f,ma,eq_dd)    

    
    points = [pNO,pOA,pAB,pBC,pCtip]
         
    constants = system.constant_values
    states=pynamics.integration.integrate_odeint(func1,ini,t, args=({'constants':constants},))

    
    points_output = PointsOutput(points,system,constant_values = constants)
    y = points_output.calc(states)
    final = numpy.asarray(states[-1,:])
    time1 = time.time()
    points_output.animate(fps = 30,movie_name = str(time1)+'video_1.mp4',lw=2,marker='o',color=(1,0,0,1),linestyle='-')
    return final,states,y,system

plt.close('all')

drag_direction = 1
[t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1
     ,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3] = fit_qs.fit_30_amount(500)
fit_t = t
fit_qA = exp_fit(fit_t,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3)
fit_qB = exp_fit(fit_t,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3)
fit_qC = exp_fit(fit_t,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3)
fit_qAd1 = numpy.diff(fit_qA)/numpy.diff(fit_t)
fit_qAd = numpy.append(fit_qAd1[0],fit_qAd1)
fit_qBd1 = numpy.diff(fit_qB)/numpy.diff(fit_t)
fit_qBd = numpy.append(fit_qBd1[0],fit_qBd1)
fit_qCd1 = numpy.diff(fit_qC)/numpy.diff(fit_t)
fit_qCd = numpy.append(fit_qCd1[0],fit_qCd1)

fit_states1 = numpy.stack((fit_qA,fit_qB,fit_qC,fit_qAd,fit_qBd,fit_qCd),axis=1)
fit_states1[:,0:3] = fit_states1[:,0:3]-fit_states1[0,0:3]
fit_states2 = -drag_direction*numpy.deg2rad(fit_states1)

distance = 200/1000
tinitial = 0
tfinal = 10
tstep = 1/30
time_step = 1/30
t = numpy.r_[tinitial:tfinal:tstep]
vel1 = 20/1000 
fit_vel = vel1*drag_direction

fit_states = numpy.array([0,0, 0, 0, 0, fit_vel ,0 , 0,  0 ,  0])
fit_states[2:5] = fit_states2[0,0:3]
fit_states[7::] = fit_states2[0,3:]
final1,states1,y1,system1  = Cal_system(fit_states,drag_direction,tinitial,tstep,tfinal,fit_vel,1.2,-0.2)

states = states1
y = y1
plt.figure()
plt.plot(*(y[::int(len(y)/20)].T)*1000)    
plt.axis('equal')
plt.title("Plate Configuration vs Distance")
plt.xlabel("Configuration")
plt.ylabel("Distance (mm)")

plt.figure()
q_states = numpy.c_[(states[:,2],states[:,3],states[:,4])]    
plt.plot(t,numpy.rad2deg(q_states) )
plt.title("Joint Angels over Time")
plt.xlabel("Time (s)")
plt.ylabel("Joint Angles (deg)")
plt.legend(["Joint 1", "Joint 2", "Joint 3"])


drag_direction = -1
vel1 = 20/1000
fit_vel = vel1*drag_direction
final1[5] = fit_vel
final1[7::] = 0
final3,states3,y3,system3  = Cal_system(final1,drag_direction,tinitial,tstep,tfinal,fit_vel,1.2,-0.2)
states = states3
y = y3
plt.figure()
plt.axis('equal')
plt.plot(*(y[::int(len(y)/20)].T)*1000)    
plt.axis('equal')
plt.title("Plate Configuration vs Distance")
plt.xlabel("Configuration")
plt.ylabel("Distance (mm)")

plt.figure()
q_states = numpy.c_[(states[:,2],states[:,3],states[:,4])]    
plt.plot(t,numpy.rad2deg(q_states) )
plt.title("Joint Angels over Time")
plt.xlabel("Time (s)")
plt.ylabel("Joint Angles (deg)")
plt.legend(["Joint 1", "Joint 2", "Joint 3"])

