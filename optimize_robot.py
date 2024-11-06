import sys
import pynamics
from pynamics.frame import Frame
from pynamics.variable_types import Differentiable,Constant
from pynamics.system import System
from pynamics.body import Body
from pynamics.dyadic import Dyadic
from pynamics.output import Output,PointsOutput
from pynamics.particle import Particle
import pynamics.integration
import pynamics.tanh as tanh

import logging
import sympy
#import sympy
import numpy
import matplotlib.pyplot as plt
from math import pi
import scipy
import scipy.optimize
from sympy import sin,cos
import cma
import time
from multiprocessing import Pool


def Cal_robot(direction,given_l,given_arm_l,omega1,t1,t2,ini_states,name1,video_on,x1,x2):
    system1 = System()
    time_a = time.time()
    pynamics.set_system(__name__,system1)
    given_k,given_b = x1
    f1,f2,f3 = x2
    global_q = True
    
    damping_r = 0
    tinitial = 0
    tfinal =  (t1-t2)/omega1
    tstep = 1/50
    t = numpy.r_[tinitial:tfinal:tstep]
    
    tol_1 = 1e-6
    tol_2 = 1e-6
    lO = Constant(27.5/1000,'lO',system1)
    # given_arm_l = 65
    lR = Constant(given_arm_l/1000,'lR',system1)
    # lR = Constant(40.5/1000,'lR',system11)
    lA = Constant(given_l/1000,'lA',system1)
    lB = Constant(given_l/1000,'lB',system1)
    lC = Constant(given_l/1000,'lC',system1)
    
    mO = Constant(1e0*154.5/1000 ,'mO',system1)
    # mR = Constant(9.282/1000   ,'mR',system1)
    mR = Constant((1.158+0.1445*given_arm_l)/1000,'mR',system1)    
    mA = Constant(given_l*2.75*0.14450000000000002/1000,'mA',system1)
    mB = Constant(given_l*2.75*0.14450000000000002/1000,'mB',system1)
    mC = Constant(given_l*2.75*0.14450000000000002/1000,'mC',system1)
    k  = Constant(given_k,'k',system1)
    
    friction_perp       = Constant(f1 ,'f_perp',system1)
    friction_par        = Constant(-1,'f_par',system1)
    friction_arm_perp   = Constant(2+given_arm_l*f3 ,'fr_perp',system1)
    friction_arm_par    = Constant(-0.3,'fr_par',system1)
    b_damping           = Constant(given_b ,'b_damping',system1)
    
    preload0 =  Constant(0*pi/180,'preload0',system1)
    preload1 =  Constant(0*pi/180,'preload1',system1)
    preload2 =  Constant(0*pi/180,'preload2',system1)
    preload3 =  Constant(0*pi/180,'preload3',system1)
    
    Ixx_O = Constant(1,'Ixx_O',system1)
    Iyy_O = Constant(1,'Iyy_O',system1)
    Izz_O = Constant(1,'Izz_O',system1)
    Ixx_R = Constant(1,'Ixx_R',system1)
    Iyy_R = Constant(1,'Iyy_R',system1)
    Izz_R = Constant(1,'Izz_R',system1)
    Ixx_A = Constant(1,'Ixx_A',system1)
    Iyy_A = Constant(1,'Iyy_A',system1)
    Izz_A = Constant(1,'Izz_A',system1)
    Ixx_B = Constant(1,'Ixx_B',system1)
    Iyy_B = Constant(1,'Iyy_B',system1)
    Izz_B = Constant(1,'Izz_B',system1)
    Ixx_C = Constant(1,'Ixx_C',system1)
    Iyy_C = Constant(1,'Iyy_C',system1)
    Izz_C = Constant(1,'Izz_C',system1)
    
    y,y_d,y_dd    = Differentiable('y',system1)
    qO,qO_d,qO_dd = Differentiable('qO',system1)
    qR,qR_d,qR_dd = Differentiable('qR',system1)
    qA,qA_d,qA_dd = Differentiable('qA',system1)
    qB,qB_d,qB_dd = Differentiable('qB',system1)
    qC,qC_d,qC_dd = Differentiable('qC',system1)
    
    initialvalues = {}
    initialvalues[y]    = ini_states[0] +tol_1
    initialvalues[qO]   = ini_states[1] +tol_1
    initialvalues[qR]   = ini_states[2] +tol_1
    initialvalues[qA]   = ini_states[3] +tol_1
    initialvalues[qB]   = ini_states[4] +tol_1
    initialvalues[qC]   = ini_states[5] +tol_1
    
    initialvalues[y_d]  = ini_states[6]  +tol_1
    initialvalues[qO_d] = ini_states[7]  +tol_1
    initialvalues[qR_d] = ini_states[8]  +tol_1
    initialvalues[qA_d] = ini_states[9]  +tol_1
    initialvalues[qB_d] = ini_states[10] +tol_1
    initialvalues[qC_d] = ini_states[11] +tol_1
    
    statevariables = system1.get_state_variables()
    ini = [initialvalues[item] for item in statevariables]
    N = Frame('N')
    O = Frame('O')
    R = Frame('R')
    A = Frame('A')
    B = Frame('B')
    C = Frame('C')
    
    system1.set_newtonian(N)
    if not global_q:
        O.rotate_fixed_axis_directed(N,[0,0,1],qO,system1)
        R.rotate_fixed_axis_directed(O,[0,0,1],qR,system1)
        A.rotate_fixed_axis_directed(O,[0,0,1],qA,system1)    
        B.rotate_fixed_axis_directed(A,[0,0,1],qB,system1)
        C.rotate_fixed_axis_directed(B,[0,0,1],qC,system1)
    else:
        O.rotate_fixed_axis_directed(N,[0,0,1],qO,system1)
        R.rotate_fixed_axis_directed(N,[0,0,1],qR,system1)    
        A.rotate_fixed_axis_directed(N,[0,0,1],qA,system1)
        B.rotate_fixed_axis_directed(N,[0,0,1],qB,system1)
        C.rotate_fixed_axis_directed(N,[0,0,1],qC,system1)
    
    pNO=    0*N.x + y*N.y
    pOR=    pNO   +lO*N.x
    # pOR=    pNO   +lO*N.x
    pRA=    pOR   +lR*R.x
    pAB=    pRA   +lA*A.x
    pBC =   pAB   +lB*B.x
    pCtip = pBC   +lC*C.x
    
    # pOcm=pNO +lO/2*N.x
    pOcm=pNO
    pRcm=pOR +lR/2*R.x
    pAcm=pRA +lA/2*A.x
    pBcm=pAB +lB/2*B.x
    pCcm=pBC +lC/2*C.x
    
    wNO = N.getw_(O)
    wOR = N.getw_(R)
    wRA = R.getw_(A)
    wAB = A.getw_(B)
    wBC = B.getw_(C)
    
    IO = Dyadic.build(O,Ixx_O,Iyy_O,Izz_O)
    IR = Dyadic.build(R,Ixx_R,Iyy_R,Izz_R)
    IA = Dyadic.build(A,Ixx_A,Iyy_A,Izz_A)
    IB = Dyadic.build(B,Ixx_B,Iyy_B,Izz_B)
    IC = Dyadic.build(C,Ixx_C,Iyy_C,Izz_C)
    
    
    BodyO = Body('BodyO',O,pOcm,mO,IO,system1)
    BodyR = Body('BodyR',R,pRcm,mR,IR,system1)
    BodyA = Body('BodyA',A,pAcm,mA,IA,system1)
    BodyB = Body('BodyB',B,pBcm,mB,IB,system1)
    BodyC = Body('BodyC',C,pCcm,mC,IC,system1)
    
    j_tol = 0*pi/180
    inv_k = 1e2
    alw = 1
    system1.add_spring_force1(k +inv_k*(qA-qR+alw*abs(qA-qR-j_tol)),(qA-qR-preload1)*N.z,wRA) 
    system1.add_spring_force1(k +inv_k*(qB-qA+alw*abs(qB-qA-j_tol)),(qB-qA-preload2)*N.z,wAB)
    system1.add_spring_force1(k +inv_k*(qC-qB+alw*abs(qC-qB-j_tol)),(qC-qB-preload3)*N.z,wBC)
    # system1.add_spring_force1(k,(qA-qR-preload1)*N.z,wRA) 
    # system1.add_spring_force1(k,(qB-qA-preload2)*N.z,wAB)
    # system1.add_spring_force1(k,(qC-qB-preload3)*N.z,wBC)   
    

    vOcm = y_d*N.y
    # vOcm = pOcm.time_derivative()
    vRcm = pRcm.time_derivative()
    vAcm = pAcm.time_derivative()
    vBcm = pBcm.time_derivative()
    vCcm = pCcm.time_derivative()
    
    nvRcm = 1/(vRcm.length()+tol_1)*vRcm
    nvAcm = 1/(vAcm.length()+tol_1)*vAcm
    nvBcm = 1/(vBcm.length()+tol_1)*vBcm
    nvCcm = 1/(vCcm.length()+tol_1)*vCcm
    
    vSoil = -direction*1*N.y
    nSoil = 1/vSoil.length()*vSoil
    # reff = abs( abs(y_d+0.01)-abs(y_d-0.01))*1/0.02*9.75
    # foperp = reff*nSoil
    foperp = f2*nSoil
    system1.addforce(foperp,vOcm)
    # system1.addforce(9.75*1*nSoil,vOcm)
    # system1.addforce(9.75*1*nSoil*y_d,vOcm)
    # system1.addforce(-100*N.x, y_d*N.x) 
    
    frperp = friction_arm_perp*nvRcm.dot(R.y)*R.y
    frpar  = friction_arm_par*nvRcm.dot(R.x)*R.x
    system1.addforce(-(frperp+frpar),vRcm) 
    
    faperp  = friction_perp*nvAcm.dot(A.y)*A.y
    fapar   = friction_par*nvAcm.dot(A.x)*A.x
    system1.addforce(-(faperp+fapar),vAcm) 
    
    fbperp = friction_perp*nvBcm.dot(B.y)*B.y
    fbpar  = friction_par*nvBcm.dot(B.x)*B.x
    system1.addforce(-(fbperp+fbpar),vBcm) 
    
    fcperp = friction_perp*nvCcm.dot(C.y)*C.y
    fcpar  = friction_par*nvCcm.dot(C.x)*C.x
    system1.addforce(-(fcperp+fcpar),vCcm) 
    
    system1.addforce(-b_damping*1*wRA,wRA)
    system1.addforce(-b_damping*1*wAB,wAB)
    system1.addforce(-b_damping*1*wBC,wBC)
    
    eq = []
    eq_d=[(system1.derivative(item)) for item in eq]
    eq_d.append(qR_d - omega1)
    # eq_dd= [(system1.derivative(item)) for item in eq_d]
    eq_dd = [system1.derivative(eq_d[0])]
    
    f,ma = system1.getdynamics()
    func1 = system1.state_space_post_invert(f,ma,eq_dd)
    points = [pNO,pOR,pRA,pAB,pBC,pCtip]
        
    constants = system1.constant_values
    states=pynamics.integration.integrate_odeint(func1,ini,t, args=({'constants':constants},))
    final = numpy.asarray(states[-1,:])
    
    logger1 = logging.getLogger('pynamics.system1')
    logger2 = logging.getLogger('pynamics.integration')
    logger3 = logging.getLogger('pynamics.output')
    logger1.disabled = True
    logger2.disabled = True
    logger3.disabled = True    
    points_output = PointsOutput(points,system1,constant_values = constants)
    
    
    y1 = points_output.calc(states)
    if video_on ==1:
        plt.figure()
        plt.plot(*(y1[::int(len(y1)/20)].T)*1000)
        plt.axis('equal')
        plt.axis('equal')
        plt.title("Plate Configuration vs Distance")
        plt.xlabel("Configuration")
        plt.ylabel("Distance (mm)")     
    
        plt.figure()
        plt.plot(t,numpy.rad2deg(states[:,2]))
        plt.plot(t,numpy.rad2deg(states[:,8]))
        plt.legend(["qR","qR_d"])
        plt.hlines(numpy.rad2deg(t1),tinitial,tfinal)
        plt.hlines(numpy.rad2deg(t2),tinitial,tfinal)
        plt.title("Robot Arm angle and velocitues (qR and qR_d) over Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Angles,Velocities (deg, deg/s)")
               
       
        plt.figure()
        q_states = numpy.c_[(states[:,2],states[:,3],states[:,4],states[:,5])]
        plt.plot(t,numpy.rad2deg(q_states) )
        plt.title("Joint Angule over Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Joint Angles (deg)")
        plt.legend(["Arm","Joint 1", "Joint 2", "Joint 3"])
        
        plt.figure()
        qd_states = numpy.c_[(states[:,8],states[:,9],states[:,10],states[:,11])]
        plt.plot(t,numpy.rad2deg(qd_states) )
        plt.legend(["qR_d","qA_d","qB_d","qC_d"])
        plt.title("Joint Angular Velocities over Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Joint Angular Velocities (deg/s)")
        plt.legend(["Arm","Joint 1", "Joint 2", "Joint 3"])
        
        plt.figure()
        plt.plot(t,states[:,0],'--') 
        plt.plot(t,states[:,6]) 
        y_d1 = states[:,6]
        # force1 = abs( abs(y_d1+0.1)-abs(y_d1-0.1))*40
        # plt.plot(t,force1) 
        plt.title("Robot Distance and Velocity over time")
        plt.xlabel("Time (s)")
        plt.ylabel("Distance (mm)")
        plt.legend(["Distance","Velocity of the robot"])       
        
        points_output.animate(fps = 1/tstep,movie_name = name1,lw=2,marker='o',color=(1,0,0,1),linestyle='-')
    else:
        pass
    return final,states,y1

def cal_eff(x,x1,x2,given_arm_l,video_flag):
    video_on = video_flag
    given_l = 20
    t1,t2 = x
    # if (0<=t1) & (t1<=90):
    #     error1 = 0
    # else:
    #     print("Theta_1 should <=90!!! and >=0!!!")
    #     error1 = 1e5
    # if (-90<=t2) & (t2<=0):
    #     error2 = 0
    # else:
    #     error2 = 1e5   
    #     print("Theta_2 should <=0!!! and >=-50!!!")
    # if error1+error2 == 0:
    logger1 = logging.getLogger('pynamics.system')
    logger2 = logging.getLogger('pynamics.integration')
    logger3 = logging.getLogger('pynamics.output')

    logger1.disabled = True
    logger2.disabled = True
    logger3.disabled = True   
    w = 66
    direction  = 1
    omega1 = w *pi/180
    omega = direction*omega1
    a1 = numpy.deg2rad(t1)
    a2 = numpy.deg2rad(t2)
    virtual_vel = 0.000001
    ini_states = numpy.array([0,0, a2, a2, a2, a2 ,virtual_vel*direction,omega, omega, 0,0,0])
    # system1 = System()
    final1,states1,y1 = Cal_robot(direction,given_l,given_arm_l,omega,a1,a2,ini_states,'robot_p1.mp4',video_on,x1,x2)
    
    direction  = -1
    omega = direction*omega1
    
    final = final1
    final[7]  = omega        
    final[8]  = omega
    # final[9::] = 0

    # system2 = System()
    final2,states2,y2 = Cal_robot(direction,given_l,given_arm_l,omega,a2,a1,final,'robot_p2.mp4',video_on,x1,x2)
    
    dis1 = states1[:,0]
    dis2 = states2[:,0]
    dis = numpy.append(dis1,dis2)
    real_dis = abs(dis[0]-dis[-1])
    forward_dis = abs(dis1[0]-dis1[-1])
    backward_dis = abs(dis2[0]-dis2[-1])
    ieta = real_dis/abs(dis2[0]-dis2[-1])
    # else:
    #     ieta = 1
    if video_on ==1:
        plt.figure()
        plt.plot(dis*1000)
        plt.title("Robot distance over time")
        plt.ylabel("Distance (mm)")
    else:
        pass
    total_eta = ieta
    return total_eta,forward_dis,backward_dis

def cal_eff1(x,x1,x2):
    total_eta,forward_dis,backward_dis = cal_eff(x,x1,x2,0)
    return total_eta


def swimmig_fit(x):
    # time.sleep(sleep_secs)
    x1 = x[0:2]
    x2 = x[2::]
    # print(x)
    
    total_eta1,forward_dis1,backward_dis1 = cal_eff([30,  0],x1,x2,65,0)    
    total_eta2,forward_dis2,backward_dis2 = cal_eff([50,  0],x1,x2,65,0)
    total_eta3,forward_dis3,backward_dis3 = cal_eff([90,  0],x1,x2,65,0)
    total_eta4,forward_dis4,backward_dis4 = cal_eff([ 0,-30],x1,x2,65,0)  
    total_eta5,forward_dis5,backward_dis5 = cal_eff([ 0,-50],x1,x2,65,0)
    total_eta6,forward_dis6,backward_dis6 = cal_eff([ 0,-90],x1,x2,65,0)
    total_eta7,forward_dis7,backward_dis7 = cal_eff([50,-50],x1,x2,65,0)   
    total_eta8,forward_dis8,backward_dis8 = cal_eff([30,-30],x1,x2,65,0)  
    forward_dis_total  = [forward_dis1,forward_dis2,forward_dis3,forward_dis4,forward_dis5,forward_dis6,forward_dis7,forward_dis8]
    backward_dis_total = [backward_dis1,backward_dis2,backward_dis3,backward_dis4,backward_dis5,backward_dis6,backward_dis7,backward_dis8]
    total_eta_total  = [total_eta1,total_eta3,total_eta4,total_eta4,total_eta5,total_eta6,total_eta7,total_eta8]
    net_dis= numpy.asarray([0.1440,0.1466,0.1250,0.1375,0.1272,0.1334,0.1457,0.1146])-numpy.asarray([ 0.1439,0.1457,0.1211,0.1316,0.1275,0.1336,0.1424,0.1132])
    eta_exp= net_dis/numpy.asarray([0.1440,0.1466,0.1250,0.1375,0.1272,0.1334,0.1457,0.1146])
    error_3 = total_eta_total - eta_exp
    total_error = numpy.sum(error_3**2)
    return total_error


def swimmig_fit1(para):
    a1,a2,given_arm_l = para
    x = [0.5878620470825995, 0.4399450573592256, 6.7792124058839685, 1.2306158325822283, 17.939121205840244]
    x1 = x[0:2]
    x2 = x[2::]
    # print(x)    
    total_eta1,forward_dis1,backward_dis1 = cal_eff([a1,  a2],x1,x2,given_arm_l,0)    
    total_eta  = 1-total_eta1
    return total_eta



options = {'bounds':[[1,-90,65],[90,0,85]]}
x0 = [50,-50,75]
es = cma.fmin(swimmig_fit1,x0,10,options)
