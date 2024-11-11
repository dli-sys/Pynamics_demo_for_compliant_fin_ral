# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 21:30:03 2020

@author: dongting
"""

import pynamics
from pynamics.output import PointsOutput
import pynamics.integration
import logging
import numpy
import matplotlib.pyplot as plt
import cma

def init_system(v,drag_direction,time_step):
    from pynamics.frame import Frame
    from pynamics.variable_types import Differentiable,Constant
    from pynamics.system import System
    from pynamics.body import Body
    from pynamics.dyadic import Dyadic
    import pynamics.integration
    import numpy
    from math import pi

    from fit_qs import exp_fit
    import fit_qs
    
    # time_step = tstep
    x = numpy.zeros((7,1))
    friction_perp=   x[0]
    friction_par=    x[1]
    given_b        = x[2]
    given_k     =    x[3]
    given_k1     =   x[4]
    given_b1     =   x[4]
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
    k = Constant(0.209,'k',system)
    k1 = Constant(0.209,'k1',system)
    
    friction_perp = Constant(1.2,'f_perp',system)
    friction_par = Constant(-0.2,'f_par',system)
    b_damping = Constant(given_b,'b_damping',system)
    
    # time_step = 1/00
   
    if v == 0:
        [t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1
     ,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3] = fit_qs.fit_0_amount(time_step)
    elif v == 10:
            [t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1
     ,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3] = fit_qs.fit_10_amount(time_step)
    elif v == 20:
            [t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1
     ,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3] = fit_qs.fit_20_amount(time_step)
    elif v == 30:
            [t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1
     ,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3] = fit_qs.fit_30_amount(time_step)
    elif v == 40:
        [t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1
     ,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3] = fit_qs.fit_40_amount(time_step)
    elif v == 50:
            [t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1
     ,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3] = fit_qs.fit_50_amount(time_step)
    
    
    distance = 200/1000
    
    nums = int(tfinal/tstep)
    array_num = numpy.arange(0,nums)
    array_num1 = numpy.repeat(array_num,nums,axis=0)
    array_num1.shape = (nums,nums)
    error_k = array_num1/8000+ numpy.ones((nums,nums))
    
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
    fit_states = -drag_direction*numpy.deg2rad(fit_states1)
    
    # plt.plot(t,fit_states)
    
    
    
    if drag_direction== -1:
        zero_shape = fit_states.shape
        fit_states = numpy.zeros(zero_shape)
    
    fit_vel = drag_direction*distance/(tfinal)
    
    if qAa1 ==0:
        fit_vel = 0
    fit_v = numpy.ones(t.shape)*fit_vel
    
    if qAa1 ==0:
        fit_d = numpy.ones(t.shape)*fit_vel
    else:
        fit_d = drag_direction*numpy.r_[tinitial:distance:tstep*abs(fit_vel)]
    
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
    
    initialvalues = {}
    initialvalues[y]=0 +1e-14
    initialvalues[y_d] = fit_vel +1e-14
    initialvalues[qO]   = 0 +1e-14
    initialvalues[qO_d] = 0 +1e-14
    initialvalues[qA]   =fit_states[0,0] +1e-14
    initialvalues[qA_d] =fit_states[0,3] +1e-14
    initialvalues[qB]   =fit_states[0,1] +1e-14
    initialvalues[qB_d] =fit_states[0,4] +1e-14
    initialvalues[qC]   =fit_states[0,2] +1e-14
    initialvalues[qC_d] =fit_states[0,5] +1e-14
    
    statevariables = system.get_state_variables()
    ini = [initialvalues[item] for item in statevariables]
    
    N = Frame('N')
    O = Frame('O')
    A = Frame('A')
    B = Frame('B')
    C = Frame('C')
    
    drag_direction =drag_direction
    velocity = 200/tfinal/1000
    vSoil = drag_direction*velocity*N.y
    nSoil = 1/vSoil.length()*vSoil
    
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
    # BodyC = Particle(pCcm,mC,'ParticleC',system)
    
    vOcm = pOcm.time_derivative()
    vAcm = pAcm.time_derivative()
    vBcm = pBcm.time_derivative()
    vCcm = pCcm.time_derivative()
    
    system.add_spring_force1(k1+10000*(qA+abs(qA)),(qA-qO-preload1)*N.z,wOA) 
    system.add_spring_force1(k +10000*(qB+abs(qB)),(qB-qA-preload2)*N.z,wAB)
    system.add_spring_force1(k +10000*(qC+abs(qC)),(qC-qB-preload3)*N.z,wBC)
    
    #new Method use nJoint
    nvAcm = 1/vAcm.length()*vAcm
    nvBcm = 1/vBcm.length()*vBcm
    nvCcm = 1/vCcm.length()*vCcm
    
    faperp  = friction_perp*nvAcm.dot(A.y)*A.y
    fapar   = friction_par*nvAcm.dot(A.x)*A.x
    system.addforce(-(faperp+fapar),vAcm) 
    
    fbperp = friction_perp*nvBcm.dot(B.y)*B.y
    fbpar  = friction_par*nvBcm.dot(B.x)*B.x
    system.addforce(-(fbperp+fbpar),vBcm) 
    
    fcperp = friction_perp*nvCcm.dot(C.y)*C.y
    fcpar  = friction_par*nvCcm.dot(C.x)*C.x
    system.addforce(-(fcperp+fcpar),vCcm) 
    
    system.addforce(-b_damping*wOA,wOA)
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

    return system,f,ma,func1,points,t,ini,constants,b_damping,k,k1,tstep,fit_states

def cal_system(system00,f00,ma00,func100,points00,t00,ini00,constants00,b_damping00,k00,x):
    import pynamics
    g_b_damping= x
    system = system00
    func1 = func100
    ini = ini00
    t = t00
    constants = system.constant_values
    constants[b_damping00] = g_b_damping
    constants[k00] =  0.209
    states=pynamics.integration.integrate_odeint(func1,ini,t, args=({'constants':constants},))
    return states, constants


def post_process(states10,constants10,points10,system10,t10,fit_states10,vel):
    states =        states10
    constants =  constants10
    points =        points10
    system  =       system10
    t =                  t10
    fit_states =fit_states10
    
    points_output = PointsOutput(points,system,constant_values = constants)
    y1 = points_output.calc(states)
    plt.figure()
    plt.plot(*(y1[::int(len(y1)/20)].T)*1000)
    plt.axis('equal')
    plt.title(str(vel))
    
    plt.figure()
    q_states = numpy.c_[(states[:,2],states[:,3],states[:,4],states[:,7],states[:,8],states[:,9])]    
    plt.plot(t,numpy.rad2deg(q_states) )
    plt.plot(t,numpy.rad2deg(fit_states),'--')
    print('final states:' +str(numpy.rad2deg(q_states[-1,:])))
    plt.title(str(vel))
    # plt.figure()
    points_output.animate(fps = 1/tstep,movie_name = 'render''.mp4',lw=2,marker='o',color=(1,0,0,1),linestyle='-')
    return y1
    
# def my_error_3var(x,constants10,points10,system10,fit_states10,b_damping10,k10,k110,func110,ini10,t10):
#     g_k,g_k1,g_b_damping= x
#     # g_k,g_b_damping= x   
#     constants =  constants10
#     points =        points10
#     system  =       system10
#     fit_states =fit_states10
#     b_damping =  b_damping10
#     k =                  k10
#     k1 =                 k110
#     func1 =          func110
#     ini =              ini10
#     t =                  t10
#     constants = system.constant_values
#     constants[b_damping] = g_b_damping
#     constants[k] =  g_k
#     constants[k1] =  g_k1
#     states=pynamics.integration.integrate_odeint(func1,ini,t, args=({'constants':constants},))     
#     q_states = numpy.c_[(states[:,2],states[:,3],states[:,4],states[:,7],states[:,8],states[:,9])] 
#     error_states = (q_states-fit_states)**2
#     mse1 = ((q_states - fit_states10)**2).mean(axis=0)
#     mse  =mse1.sum()
#     return mse


def my_error_single_justb(x,constants10,points10,system10,fit_states10,b_damping10,k10,func110,ini10,t10):
    g_b_damping= x   
    constants =  constants10
    points =        points10
    system  =       system10
    fit_states =fit_states10
    b_damping =  b_damping10
    k =                  k10
    func1 =          func110
    ini =              ini10
    t =                  t10
    constants = system.constant_values
    constants[b_damping] = g_b_damping
    constants[k] =  0.209
    states=pynamics.integration.integrate_odeint(func1,ini,t, args=({'constants':constants},))     
    q_states = numpy.c_[(states[:,2],states[:,3],states[:,4],states[:,7],states[:,8],states[:,9])] 
    mse1 = ((q_states - fit_states10)**2).mean(axis=0)
    mse  =mse1.sum()
    pos_error = mse1[0:3].sum()
    return mse

def my_error_sum(x):
    error10 = my_error_single_justb(x,constants10,points10,system10,fit_states10,b_damping10,k10,func110,ini10,t10)
    error20 = my_error_single_justb(x,constants20,points20,system20,fit_states20,b_damping20,k20,func120,ini20,t20)
    error30 = my_error_single_justb(x,constants30,points30,system30,fit_states30,b_damping30,k30,func130,ini30,t30)
    error40 = my_error_single_justb(x,constants40,points40,system40,fit_states40,b_damping40,k40,func140,ini40,t40)
    error50 = my_error_single_justb(x,constants50,points50,system50,fit_states50,b_damping50,k50,func150,ini50,t50)
    error_sum = (error10+error20+error30+error40+error50)/5
    
    return error_sum

tstep = 500
drag_direction = 1


system10,f10,ma10,func110,points10,t10,ini10,constants10,b_damping10,k10,k110,tstep10,fit_states10 = init_system(10,drag_direction,tstep)
system20,f20,ma20,func120,points20,t20,ini20,constants20,b_damping20,k20,k120,tstep20,fit_states20 = init_system(20,drag_direction,tstep)
system30,f30,ma30,func130,points30,t30,ini30,constants30,b_damping30,k30,k130,tstep30,fit_states30 = init_system(30,drag_direction,tstep)
system40,f40,ma40,func140,points40,t40,ini40,constants40,b_damping40,k40,k140,tstep40,fit_states40 = init_system(40,drag_direction,tstep)
system50,f50,ma50,func150,points50,t50,ini50,constants50,b_damping50,k50,k150,tstep50,fit_states50 = init_system(50,drag_direction,tstep)


logger1 = logging.getLogger('pynamics.system')
logger2 = logging.getLogger('pynamics.integration')
logger3 = logging.getLogger('pynamics.output')
logger1.disabled = True
logger2.disabled = True
logger3.disabled = True


# Single variable optimization
es = cma.CMAEvolutionStrategy([1.5,1.5], 1)
es.optimize(lambda x: my_error_sum(x[0]))
es.logger.plot(xsemilog=True)
ans = numpy.asarray(es.result.xbest)
print(ans)

x  = es.result.xbest[0]

states10,constants10 = cal_system(system10,f10,ma10,func110,points10,t10,ini10,constants10,b_damping10,k10,x)
states20,constants20 = cal_system(system20,f20,ma20,func120,points20,t20,ini20,constants20,b_damping20,k20,x)
states30,constants30 = cal_system(system30,f30,ma30,func130,points30,t30,ini30,constants30,b_damping30,k30,x)
states40,constants40 = cal_system(system40,f40,ma40,func140,points40,t40,ini40,constants40,b_damping40,k40,x)
states50,constants50 = cal_system(system50,f50,ma50,func150,points50,t50,ini50,constants50,b_damping50,k50,x)


plt.close('all')
y10 = post_process(states10,constants10,points10,system10,t10,fit_states10,10)
y20 = post_process(states20,constants20,points20,system20,t20,fit_states20,20)
y30 = post_process(states30,constants30,points30,system30,t30,fit_states30,30)
y40 = post_process(states40,constants40,points40,system40,t40,fit_states40,40)
y50 = post_process(states50,constants50,points50,system50,t50,fit_states50,50)
