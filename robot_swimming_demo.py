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
import scipy
import scipy.optimize
from sympy import sin,cos
import pynamics.tanh as tanh
import cma
import time



def Cal_robot(direction,given_l,omega1,t1,t2,ini_states,name1,system,video_on,x1):
    time_a = time.time()
    pynamics.set_system(__name__,system)
    given_k,given_b = x1
    global_q = True
    
    damping_r = 0
    tinitial = 0
    tfinal =  (t1-t2)/omega1
    tstep = 1/30
    t = numpy.r_[tinitial:tfinal:tstep]
    
    tol_1 = 1e-6
    tol_2 = 1e-6
    lO = Constant(27.5/1000,'lO',system)
    lR = Constant(40.5/1000,'lR',system)
    lA = Constant(given_l/1000,'lA',system)
    lB = Constant(given_l/1000,'lB',system)
    lC = Constant(given_l/1000,'lC',system)
    
    mO = Constant(154.5/1000 ,'mO',system)
    mR = Constant(9.282/1000   ,'mR',system)
    mA = Constant(given_l*2.75*0.14450000000000002/1000,'mA',system)
    mB = Constant(given_l*2.75*0.14450000000000002/1000,'mB',system)
    mC = Constant(given_l*2.75*0.14450000000000002/1000,'mC',system)
    k  = Constant(given_k,'k',system)
    
    friction_perp   = Constant(13/3 ,'f_perp',system)
    friction_par    = Constant(-2/3,'f_par',system)
    friction_arm_perp   = Constant(5.6 ,'fr_perp',system)
    friction_arm_par    = Constant(-0.2,'fr_par',system)
    b_damping       = Constant(given_b ,'b_damping',system)
    
    preload0 =  Constant(0*pi/180,'preload0',system)
    preload1 =  Constant(0*pi/180,'preload1',system)
    preload2 =  Constant(0*pi/180,'preload2',system)
    preload3 =  Constant(0*pi/180,'preload3',system)
    
    Ixx_O = Constant(1,'Ixx_O',system)
    Iyy_O = Constant(1,'Iyy_O',system)
    Izz_O = Constant(1,'Izz_O',system)
    Ixx_R = Constant(1,'Ixx_R',system)
    Iyy_R = Constant(1,'Iyy_R',system)
    Izz_R = Constant(1,'Izz_R',system)
    Ixx_A = Constant(1,'Ixx_A',system)
    Iyy_A = Constant(1,'Iyy_A',system)
    Izz_A = Constant(1,'Izz_A',system)
    Ixx_B = Constant(1,'Ixx_B',system)
    Iyy_B = Constant(1,'Iyy_B',system)
    Izz_B = Constant(1,'Izz_B',system)
    Ixx_C = Constant(1,'Ixx_C',system)
    Iyy_C = Constant(1,'Iyy_C',system)
    Izz_C = Constant(1,'Izz_C',system)
    
    y,y_d,y_dd    = Differentiable('y',system)
    qO,qO_d,qO_dd = Differentiable('qO',system)
    qR,qR_d,qR_dd = Differentiable('qR',system)
    qA,qA_d,qA_dd = Differentiable('qA',system)
    qB,qB_d,qB_dd = Differentiable('qB',system)
    qC,qC_d,qC_dd = Differentiable('qC',system)
    
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
    
    statevariables = system.get_state_variables()
    ini = [initialvalues[item] for item in statevariables]
    N = Frame('N')
    O = Frame('O')
    R = Frame('R')
    A = Frame('A')
    B = Frame('B')
    C = Frame('C')
    
    system.set_newtonian(N)
    if not global_q:
        O.rotate_fixed_axis_directed(N,[0,0,1],qO,system)
        R.rotate_fixed_axis_directed(O,[0,0,1],qR,system)
        A.rotate_fixed_axis_directed(R,[0,0,1],qA,system)    
        B.rotate_fixed_axis_directed(A,[0,0,1],qB,system)
        C.rotate_fixed_axis_directed(B,[0,0,1],qC,system)
    else:
        O.rotate_fixed_axis_directed(N,[0,0,1],qO,system)
        R.rotate_fixed_axis_directed(N,[0,0,1],qR,system)    
        A.rotate_fixed_axis_directed(N,[0,0,1],qA,system)
        B.rotate_fixed_axis_directed(N,[0,0,1],qB,system)
        C.rotate_fixed_axis_directed(N,[0,0,1],qC,system)
    
    pNO=    0*N.x + y*N.y
    pOR=    pNO   +lO*N.x
    pRA=    pOR   +lR*R.x
    pAB=    pRA   +lA*A.x
    pBC =   pAB   +lB*B.x
    pCtip = pBC   +lC*C.x
    
    pOcm=pNO +lO/2*N.x
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
    
    
    BodyO = Body('BodyO',O,pOcm,mO,IO,system)
    BodyR = Body('BodyR',R,pRcm,mR,IR,system)
    BodyA = Body('BodyA',A,pAcm,mA,IA,system)
    BodyB = Body('BodyB',B,pBcm,mB,IB,system)
    BodyC = Body('BodyC',C,pCcm,mC,IC,system)
    
    j_tol = 3*pi/180
    inv_k = 10
    system.add_spring_force1(k +inv_k*(qA-qR+abs(qA-qR-j_tol)),(qA-qR-preload1)*N.z,wRA) 
    system.add_spring_force1(k +inv_k*(qB-qA+abs(qB-qA-j_tol)),(qB-qA-preload2)*N.z,wAB)
    system.add_spring_force1(k +inv_k*(qC-qB+abs(qC-qB-j_tol)),(qC-qB-preload3)*N.z,wBC)
    
    

    vOcm = y_d*N.y
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
    foperp = 8*nSoil
    system.addforce(-foperp,vOcm) 
    
    frperp = friction_arm_perp*nvRcm.dot(R.y)*R.y
    frpar  = friction_arm_par*nvRcm.dot(R.x)*R.x
    system.addforce(-(frperp+frpar),vRcm) 
    
    faperp  = friction_perp*nvAcm.dot(A.y)*A.y
    fapar   = friction_par*nvAcm.dot(A.x)*A.x
    system.addforce(-(faperp+fapar),vAcm) 
    
    fbperp = friction_perp*nvBcm.dot(B.y)*B.y
    fbpar  = friction_par*nvBcm.dot(B.x)*B.x
    system.addforce(-(fbperp+fbpar),vBcm) 
    
    fcperp = friction_perp*nvCcm.dot(C.y)*C.y
    fcpar  = friction_par*nvCcm.dot(C.x)*C.x
    system.addforce(-(fcperp+fcpar),vCcm) 
    
    system.addforce(-b_damping*1*wRA,wRA)
    system.addforce(-b_damping*1*wAB,wAB)
    system.addforce(-b_damping*1*wBC,wBC)
    
    eq = []
    eq_d=[(system.derivative(item)) for item in eq]
    eq_d.append(qR_d - omega1)
    eq_dd=[(system.derivative(item)) for item in eq_d]
    
    f,ma = system.getdynamics()
    func1 = system.state_space_post_invert(f,ma,eq_dd)
    points = [pNO,pOR,pRA,pAB,pBC,pCtip]
        
    constants = system.constant_values
    states=pynamics.integration.integrate_odeint(func1,ini,t, args=({'constants':constants},))
    final = numpy.asarray(states[-1,:])
    
    logger1 = logging.getLogger('pynamics.system')
    logger2 = logging.getLogger('pynamics.integration')
    logger3 = logging.getLogger('pynamics.output')
    logger1.disabled = True
    logger2.disabled = True
    logger3.disabled = True    
    points_output = PointsOutput(points,system,constant_values = constants)
    
    
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
        plt.title("Robot Distance and Velocity over time")
        plt.xlabel("Time (s)")
        plt.ylabel("Distance (mm)")
        plt.legend(["Distance","Velocity of the robot"])       
        
        points_output.animate(fps = 1/tstep,movie_name = name1,lw=2,marker='o',color=(1,0,0,1),linestyle='-')
    else:
        pass
    return final,states,y1

def cal_eff(x,x1,video_flag):
    video_on = video_flag
    given_l = 20
    t1,t2 = x
    if (0<=t1) & (t1<=90):
        error1 = 0
    else:
        print("Theta_1 should <=90!!! and >=0!!!")
        error1 = 1e5
    if (-50<=t2) & (t2<=0):
        error2 = 0
    else:
        error2 = 1e5   
        print("Theta_2 should <=0!!! and >=-50!!!")
    if error1+error2 == 0:
        logger1 = logging.getLogger('pynamics.system')
        logger2 = logging.getLogger('pynamics.integration')
        logger3 = logging.getLogger('pynamics.output')

        logger1.disabled = True
        logger2.disabled = True
        logger3.disabled = True   
        w = 80
        direction  = 1
        omega1 = w *pi/180
        omega = direction*omega1
        a1 = numpy.deg2rad(t1)
        a2 = numpy.deg2rad(t2)
        virtual_vel = 0
        ini_states = numpy.array([0,0, a2, a2, a2, a2 ,virtual_vel*direction,0, omega, 0,0,0])
        system1 = System()
        final1,states1,y1 = Cal_robot(direction,given_l,omega,a1,a2,ini_states,'robot_p1.mp4',system1,video_on,x1)
        
        direction  = -1
        omega = direction*omega1
        
        final = final1
        final[8]  = omega
        final[9::] = 0
    
        system2 = System()
        final2,states2,y2 = Cal_robot(direction,given_l,omega,a2,a1,final,'robot_p2.mp4',system2,video_on,x1)
        
        dis1 = states1[:,0]
        dis2 = states2[:,0]
        dis = numpy.append(dis1,dis2)
        real_dis = abs(dis[0]-dis[-1])
        forward_dis = abs(dis1[0]-dis1[-1])
        backward_dis = abs(dis2[0]-dis2[-1])
        ieta = 1-real_dis/abs(dis2[0]-dis2[-1])
    else:
        ieta = 1
    if video_on ==1:
        plt.figure()
        plt.plot(dis*1000)
        plt.title("Robot distance over time")
        plt.ylabel("Distance (mm)")
    else:
        pass
    total_eta = ieta+error1+error2
    return total_eta,forward_dis,backward_dis

def cal_eff1(x):
    print(x)
    x1 = [1.63587104, 1.25042469]
    total_eta,forward_dis,backward_dis = cal_eff(x,x1,0)
    return total_eta


x = [85,-10]    
# theta12 = [88,-10]
cal_eff(x,[1.63587104, 1.25042469],1)