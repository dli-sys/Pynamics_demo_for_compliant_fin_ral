import sys
import os

import pynamics
from pynamics.frame import Frame
from pynamics.variable_types import Differentiable,Constant
from pynamics.system import System
from pynamics.body import Body
from pynamics.dyadic import Dyadic
from pynamics.output import Output,PointsOutput
from pynamics.particle import Particle, PseudoParticle
import pynamics.integration
from pynamics.constraint import AccelerationConstraint
from pynamics import tanh
import logging
import time

from numpy import pi
import matplotlib
from matplotlib import animation
import matplotlib.pyplot as plt
import sympy
import numpy
import scipy
import scipy.optimize
from sympy import sin, cos

# import cma

plt.close('all')
plt.ion()

def Cal_robot(system,direction, angular_vel, ini_states,force_coeff,sim_time=False,video_on=True,video_name="swimming.gif",sim_mode="gm"):

  [body_drag_lg,body_drag_wd,arm_force_prep,arm_force_par] = force_coeff

  # time_a = time.time()
  pynamics.set_system(__name__, system)

  if sim_time is False:
    tfinal = 5
  else:
    tfinal = sim_time
  tinitial = 0
  tstep = 1 / 30
  t = numpy.r_[tinitial:tfinal:tstep]

  # Add some tol so that nothing is divided by zero
  tol_1 = 1e-8
  tol_2 = 1e-8
  # Length of the body
  lO = Constant(30 / 1000, 'lO', system)
  lApater = Constant(10 / 1000, 'lA', system)
  lR = Constant(60 / 1000, 'lR', system)

  # mass of the bodies
  mO = Constant(300 / 1000, 'mO', system)
  mR = Constant(100 / 1000, 'mR', system)
  mL = Constant(100 / 1000, 'mL', system)
  # k = Constant(given_k, 'k', system)

  # Adding friction forces to the body
  friction_arm_perp = Constant(arm_force_prep, 'fr_perp', system)
  friction_arm_par = Constant(arm_force_par, 'fr_par', system)
  # b_damping = Constant(given_b, 'b_damping', system)


  Ixx_O = Constant(1, 'Ixx_O', system)
  Iyy_O = Constant(1, 'Iyy_O', system)
  Izz_O = Constant(1, 'Izz_O', system)

  # This might need toi be adjusted based on the real location of CoM
  Ixx_R = Constant(1, 'Ixx_R', system)
  Iyy_R = Constant(1, 'Iyy_R', system)
  Izz_R = Constant(1, 'Izz_R', system)

  Ixx_L = Constant(1, 'Ixx_L', system)
  Iyy_L = Constant(1, 'Iyy_L', system)
  Izz_L = Constant(1, 'Izz_L', system)

  x, x_d, x_dd = Differentiable('x', system)
  y, y_d, y_dd = Differentiable('y', system)
  qO, qO_d, qO_dd = Differentiable('qO', system)
  qR, qR_d, qR_dd = Differentiable('qR', system)
  qL, qL_d, qL_dd = Differentiable('qL', system)


  # Set initial values. Importsnat
  initialvalues = {}
  initialvalues[x]  = ini_states[0] + tol_1
  initialvalues[y]  = ini_states[1] + tol_1
  initialvalues[qO] = ini_states[2] + tol_1
  initialvalues[qR] = ini_states[3] + tol_1
  initialvalues[qL] = ini_states[4] + tol_1

  initialvalues[x_d]  = ini_states[5] + tol_1
  initialvalues[y_d]  = ini_states[6] + tol_1
  initialvalues[qO_d] = ini_states[7] + tol_1
  initialvalues[qR_d] = ini_states[8] + tol_1
  initialvalues[qL_d] = ini_states[9] + tol_1

  statevariables = system.get_state_variables()
  ini = [initialvalues[item] for item in statevariables]

  N = Frame('N',system)
  O = Frame('O',system)
  R = Frame('R',system)
  L = Frame('L',system)

  system.set_newtonian(N)

  global_q = False

  if not global_q:
    O.rotate_fixed_axis(N, [0, 0, 1], qO, system)
    R.rotate_fixed_axis(O, [0, 0, 1], qR, system)
    L.rotate_fixed_axis(O, [0, 0, -1], qL, system)

  else:
    O.rotate_fixed_axis(N, [0, 0, 1], qO, system)
    R.rotate_fixed_axis(N, [0, 0, 1], qR, system)
    L.rotate_fixed_axis(N, [0, 0, 1], qL, system)


  # pNO =  + y * N.y # For simulation on the rail
  pNO = x * O.x + y * O.y # Free swimming exp
  pOR = pNO + 0.5*lO * O.x
  # Just adding a dummry adapter that does not chaning anything so that I can calculate the location of the real force applied
  pAdapterR = pOR + lApater*R.x
  pRA = pOR + lR * R.x

  pOL = pNO - 0.5*lO * O.x
  # Just adding a dummry adapter that does not chaning anything so that I can calculate the location of the real force applied
  pAdapterL = pOL - lApater*L.x
  pLA = pOL - lR * L.x

  pOcm = pNO + lO / 2 * N.x

  pRcm = pOR + lApater * R.x + (lR-lApater) /2 * R.x
  pLcm = pOL - lApater * L.x - (lR-lApater) / 2 * L.x

  # points = [pLA, pLcm, pAdapterL, pOL, pNO, pOR, pAdapterL, pRcm, pRA]
  # pRcm = pOR + lR / 2 * R.x

  # w is the rotation velocity and R is the angle
  wNO = N.get_w_to(O)
  wOR = N.get_w_to(R)
  wOL = N.get_w_to(L)

  IO = Dyadic.build(O, Ixx_O, Iyy_O, Izz_O)
  IR = Dyadic.build(R, Ixx_R, Iyy_R, Izz_R)
  IL = Dyadic.build(L, Ixx_L, Iyy_L, Izz_L)

  BodyO = Body('BodyO', O, pOcm, mO, IO, system)
  BodyR = Body('BodyR', R, pRcm, mR, IR, system)
  BodyL = Body('BodyL', L, pLcm, mL, IL, system)

  # j_tol = 3 * pi / 180
  # inv_k = 10


  # vOcm = y_d * N.y # Fixed simulatyion
  vOcm = y_d * N.y + x_d*N.x
  vRcm = pRcm.time_derivative()
  vLcm = pLcm.time_derivative()

  nvRcm = 1 / (vRcm.length() + tol_1) * vRcm
  nvLcm = 1 / (vLcm.length() + tol_1) * vLcm
  nvOcm = 1 / (vOcm.length() + tol_1) * vOcm


  if sim_mode == "gm":
    # vSoil = -direction * 1 * N.y
    # nSoil = -1 / vSoil.length()
    # foperp = body_drag_lg * nSoil
    body_vel_factor = tanh.gen_static_friction(vOcm.length(),plot=False)
    # body_vel_factor = 1
    # Is this part real? Since RFT only applies to granular flow, so it's also state related
    # In the paper, do you want to assume that everything is moving -- basically a flow
    # Or you want to do something like velocity or state based programmming?
    system.addforce(-body_drag_lg*nvOcm*body_vel_factor, vOcm)

    fin_vel_factor = tanh.gen_static_friction(vRcm.length(), plot=False)
    # fin_vel_factor = 1
    frperp = friction_arm_perp * nvRcm.dot(R.y) * R.y
    frpar = friction_arm_par * nvRcm.dot(R.x) * R.x

    fin_friction = frperp+frpar

    system.addforce(-fin_vel_factor*fin_friction, vRcm)

    fin_vel_factor_L = tanh.gen_static_friction(vLcm.length(), plot=False)
    # fin_vel_factor = 1
    frperp = friction_arm_perp * nvLcm.dot(L.y) * L.y
    frpar = friction_arm_par * nvLcm.dot(L.x) * L.x
    fin_friction = frperp+frpar
    system.addforce(-fin_vel_factor_L*fin_friction, vLcm)

  if sim_mode =="water":
    rho = Constant(1000,'rho',system)
    Area = Constant(0.01,'Area',system)
    # f_aero_C2 = rho * vAcm.length() * (vAcm.dot(A.y)) * Area * A.y

    # this part is assumoing the force is always perep to the body, right?
    # f_aero_body  = rho*vOcm.length()*(vOcm.dot(O.y))*Area*O.y*lO
    # f_aero_fin_L = rho*vLcm.length()*(vLcm.dot(O.y))*Area*L.y*lR
    # f_aero_fin_R = rho*vRcm.length()*(vRcm.dot(O.y))*Area*R.y*lR

    f_aero_body  = rho*vOcm.dot(O.y)*(vOcm.dot(O.y))*Area*O.y*lO * 1
    f_aero_fin_L = rho*vLcm.dot(L.y)*(vLcm.dot(L.y))*Area*L.y*lR * 0.6
    f_aero_fin_R = rho*vRcm.dot(R.y)*(vRcm.dot(R.y))*Area*R.y*lR * 0.6

    system.addforce(-f_aero_body, vOcm)
    system.addforce(-f_aero_fin_L, vLcm)
    system.addforce(-f_aero_fin_R, vRcm)


  # system.addforce(0, vLcm)

  # This is using forced vel already but vel is still chanigng -- looks like eq_d - constant =0 is the old pynamics way if you wnat to using new pynamics, you want to do eq_dd
  eq = []
  eq_d = [(system.derivative(item)) for item in eq]
  eq_d.append(qR_d - ini_states[8]) # to make qR_d is zero
  eq_d.append(qL_d - ini_states[9])

  eq_dd = [(system.derivative(item)) for item in eq_d]
  # eq_dd_con = eq_dd.dot(R.x)
  eq_dd_scalar = [qR_dd,qL_dd]
  # eq_dd_scalar.append([qR_dd],[qL_dd])
  system.add_constraint(AccelerationConstraint(eq_dd_scalar))

  f, ma = system.getdynamics()


  func1 = system.state_space_post_invert(f,ma,constants = system.constant_values)

  points = [pLA, pOL, pNO, pOR, pRA]

  # points = [pLA,pLcm,pAdapterL , pOL,pNO, pOR, pAdapterL, pRcm, pRA]
  # points = [pNO, pOR,pAdapterL,pRcm,pRA]
  # points = [pNO, pOR, pRA, pAB, pBC, pCtip]

  alpha = 1e2
  beta = 1e2
  error = 1e-3
  error_tol = 1e-3

  states = pynamics.integration.integrate_odeint(func1, ini, t, rtol=1e-4, atol=1e-4,
                                                 args=({'constants': {}, 'alpha': alpha, 'beta': beta},))

  # constants = system.constant_values
  # states = pynamics.integration.integrate_odeint(func1, ini, t,rtol = error, atol = error,  args=({'alpha':alpha,'beta':beta, 'constants':system.constant_values}),full_output = 1,mxstep = int(1e5))
  # plt.plot(states[:,0])
  # final = numpy.asarray(states[-1, :])


  # logger1 = logging.getLogger('pynamics.system')
  # logger2 = logging.getLogger('pynamics.integration')
  # logger3 = logging.getLogger('pynamics.output')
  # logger1.disabled = True
  # logger2.disabled = True
  # logger3.disabled = True

  # Here is how to use points to calculatethe video
  points_output = PointsOutput(points, system, constant_values=system.constant_values)
  y1 = points_output.calc(states,t)

  if video_on:
    # plt.figure()
    # plt.show()
    # print(forward_limits)
    # plt.axis("equal")
    # plt.ion()

    fig,_axs= plt.subplots(2,2)
    [[x_dis_axis, y_dis_axis], [ang_axis, ang_axis_L]] = _axs
    x_dis_axis.plot(states[:,0]*1e3)
    x_dis_axis.plot(states[:, 5]*1e3,'--')
    x_dis_axis.set_xlabel("Time (s)")
    x_dis_axis.set_ylabel("dis/vel (mm)")
    x_dis_axis.set_title("X distance/vel vs. time")
    x_dis_axis.legend({"x", "xd"})
    x_dis_axis.grid('on')

    y_dis_axis.plot(states[:,1]*1e3)
    y_dis_axis.plot(states[:, 6]*1e3,'--')
    y_dis_axis.set_xlabel("Time (s)")
    y_dis_axis.set_ylabel("dis/vel (mm)")
    y_dis_axis.set_title("Y distance/vel vs. time")
    y_dis_axis.legend({"y", "yd"})
    y_dis_axis.grid('on')


    ang_axis.plot(numpy.rad2deg(states[:,3]))
    ang_axis.plot(numpy.rad2deg(states[:, 8]),'--')
    ang_axis.set_title("Fin angle/vel vs. time")
    ang_axis.legend({"qR","qRd"})
    ang_axis.set_xlabel("Time (s)")
    ang_axis.set_ylabel("angle/vel (deg)")
    ang_axis.grid('on')

    ang_axis_L.plot(numpy.rad2deg(states[:,4]))
    ang_axis_L.plot(numpy.rad2deg(states[:, 9]),'--')
    ang_axis_L.set_title("Fin angle/vel vs. time")
    ang_axis_L.legend({"qL","qLd"})
    ang_axis_L.set_xlabel("Time (s)")
    ang_axis_L.set_ylabel("angle/vel (deg)")
    ang_axis_L.grid('on')

    # rot_axis.plot((numpy.rad2deg(states[:,2]))
    # rot_axis.plot((numpy.rad2degstates[:, 7]),'--')
    # rot_axis.set_xlabel("Time (s)")
    # rot_axis.set_ylabel("dis/vel (mm)")
    # rot_axis.set_title("Body Rotation/vel vs. time")
    # rot_axis.legend({"R_O", "w_O"})
    # rot_axis.grid('on')


    from matplotlib import colormaps

    plt.figure()
    cmap = plt.get_cmap('bwr')
    colors = cmap(numpy.linspace(0, 1, len(y1[::int(len(y1) / 20)])))  # Generate colors from the colormap
    y1_plot = y1*1e3
    for i, line in enumerate(y1_plot[::int(len(y1_plot) / 20)]):
      plt.plot(*(line.T) , color=colors[i])
    # plt.arrow(y1[0,0],y1[-1,0])
    plt.arrow(y1_plot[0, 0][0],y1_plot[0, 0][1], y1_plot[-1, 0][0]-y1_plot[0, 0][0],y1_plot[-1, 0][1]-y1_plot[0, 0][1], head_width=2, head_length=5, fc='k', ec='k')

    # plt.plot(*(y1[::int(len(y1) / 20)].T) * 1000,cmap='bwr')
    # plt.axis('equal')
    # plt.axis('equal')
    plt.title("Plate Configuration vs time (blue - red)")
    # plt.xlabel("Configuration")
    plt.ylabel("Distance (mm)")
    plt.axis('equal')
    # plt.show()

    # points_output.animate(fps=1 / tstep, movie_name=video_name,scale=1e3, lw=2, marker='o', color=(1, 0, 0, 1), linestyle='-')
    # plt.close()
    # points_output.animate(fps=1 / tstep, movie_name=video_name, lw=2, marker='o', color=(1, 0, 0, 1), linestyle='-')
    # plt.show()
  else:
    pass
  return states, y1, points


# def cal_eff(video_flag):
#
#
#   return total_eta, forward_dis, backward_dis


if __name__ == "__main__":
  plt.close('all')
  plt.ion()

  # if error1 + error2 == 0:
  # logger1 = logging.getLogger('pynamics.system')
  # logger2 = logging.getLogger('pynamics.integration')
  # logger3 = logging.getLogger('pynamics.output')
  #
  # logger1.disabled = True
  # logger2.disabled = True
  # logger3.disabled = True
  direction = 1

  # Looks like a2 is the initial angle of the arm, vir_vel is the initial velocity, omega is some angular velocity or angle
  # DEfination
  # [y,qO,qR,y_d,qO_d,qR_d]
  # End def
  # above zero -- positive below_zero -- negtive -- counterclockwise

  servo_speed   = pi/180*10
  ini_angle     = pi/180*-60
  ini_states = numpy.array([0, 0, 0, ini_angle, ini_angle, 0, 0, 0, servo_speed,servo_speed])
  # Just add amplitude the direction is handlled inside
  fin_drag_reduction_coef   = 0.3
  body_drag_reduction_coef  = 0.6
  fin_perp    = 10
  fin_par     = -2
  body_drag_lg   = 20
  body_drag_wd = 0

  force_coeff_p = [body_drag_lg,body_drag_wd,fin_perp*fin_drag_reduction_coef,fin_par*fin_drag_reduction_coef]
  force_coeff_r = [body_drag_lg*body_drag_reduction_coef,body_drag_wd*body_drag_reduction_coef, fin_perp, fin_par]

  mode='gm'
  sim_time = 12


  system1 = System()
  states1, y1,forward_points = Cal_robot(system1,direction, servo_speed, ini_states,force_coeff_p,video_on=True,video_name='robot_p1.gif',sim_time=sim_time,sim_mode=mode)
  # plt.figure()
  # plt.plot(numpy.rad2deg(states1[:,2]))
  # plt.show()
  # plt.close()
  final =  numpy.asarray(states1[-1, :])
  # clear the velocity
  final[5::] = 0
  # DEfine the velocity
  final[-2] = -servo_speed * 1
  final[-1] = -servo_speed * 0.5
  sim_time = 12

  system2 = System()
  states2, y2,recovery_points = Cal_robot(system2,-direction, servo_speed, final, force_coeff_r,video_on=True,video_name='robot_p2.gif',sim_time=sim_time,sim_mode=mode)

  full_out_y = numpy.vstack((y1,y2))

  y  = full_out_y
  movie_name = "swimming.mp4"

  PointsOutput.point_anim(full_out_y,movie_name=movie_name,lw=2, marker='o', color=(1, 0, 0, 1), linestyle='-',title="Full swimming, unit (mm)")

  plt.figure()
  x1 = states1[:, 0] * 1e3
  x2 = states2[:, 0] * 1e3
  y1 = states1[:, 1] * 1e3
  y2 = states2[:, 1] * 1e3
  dis_x = numpy.append(x1, x2)
  dis_y = numpy.append(y1, y2)
  trajectories = []
  for x, y in zip(dis_x, dis_y):
    trajectories.append([x, y])

  cmap = plt.get_cmap('bwr')
  colors = cmap(numpy.linspace(0, 1, len(dis_x)) ) # Generate colors from the colormap_

  vis_step = abs(dis_y[-1]-dis_y[0])/len(dis_y) if abs(dis_x[-1]-dis_x[0]) < 0.1 else 0
  print(vis_step)
  for idx,trajectory in enumerate(trajectories):
    plt.plot(trajectory[0]+vis_step*idx, trajectory[1],color=colors[idx],marker='.')
  plt.axis("equal")
  #
  # real_dis = abs(dis[0] - dis[-1])
  # forward_dis = abs(dis1[0] - dis1[-1])
  # backward_dis = abs(dis2[0] - dis2[-1])
  # ieta = 1 - real_dis / abs(dis2[0] - dis2[-1])
  # print(ieta)

  # plt.figure()
  # plt.plot([dis_x,dis_y])
  # plt.title("Robot distance over time")
  # plt.ylabel("Distance (mm)")
  plt.show(block=True)
