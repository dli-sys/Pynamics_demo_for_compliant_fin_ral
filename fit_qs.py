# [t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3]  = fit_10()

import numpy

def fit_10(time_step):
    #fit for 10mm/s
    tinitial = 0
    tfinal = 20
    tstep = time_step
    t = numpy.r_[tinitial:tfinal:tstep]
    qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3 = [59.39,0.09915,1.046,49.89,0.1134,4.763,4.963,0.2814,6.505]
    qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3 = [407.5,0.007948,0.01161,17,0.2033,-0.2733,1.725,0.4184,6.076]    
    qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3 = [65.4,0.08329,-0.0815,8.104,0.2372,6.337,0.7624,0.5186,4.748]
    return t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3

def fit_20(time_step):
#fit    for 20mm/s
    tinitial = 0
    tfinal = 10
    tstep = time_step
    t = numpy.r_[tinitial:tfinal:tstep]
    qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3 = [54.05,0.2079,0.8322,37.66,0.2473,4.563,3.95,0.6362,6.234]
    qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3 = [73.64,0.247,-0.175,357.2,0.5185,2.136,345,0.5222,5.283]
    qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3 = [58.86,0.2063,-0.05415,17.15,1.042,-2.215,15.12,1.081,0.7287]
    return t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3

def fit_30(time_step):
    #fit for 30mm/s
    tinitial = 0
    tfinal = 200/30
    tstep = time_step
    t = numpy.r_[tinitial:tfinal:tstep]
    qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3 = [88.94,0.4628,-0.1873,128.6,0.6458,2.73,78.07,0.7099,5.841]
    qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3 = [64.41,0.251,0.03774,4.598,1.008,-0.6726,0.4017,1.85,-0.3148]
    qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3 = [59.27,0.2861,-0.04419,12,1.432,-1.887,9.104,1.546,0.9067]
    return t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3

def fit_40(time_step):
    #fit for 40mm/s
    tinitial = 0
    tfinal = 5
    tstep = time_step
    t = numpy.r_[tinitial:tfinal:tstep]
    qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3 = [41.65,0.4353,0.2967,394.5,0.9561,3.474,386.8,0.9674,6.597]
    qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3 = [55.94,0.3319,0.269,76.83,0.9285,3.845,72.52,0.9913,6.798]
    qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3 = [54.53,0.401,0.4147,104.6,0.8511,3.842,88.3,0.9431,6.758]
    return t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3

def fit_50(time_step):
    #fit for 50mm/s
    tinitial = 0
    tfinal = 4
    tstep = time_step
    t = numpy.r_[tinitial:tfinal:tstep]
    qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3 = [76.13,0.77,-0.1489,184.6,1.185,2.624,145.9,1.248,5.722]
    qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3 = [77.78,0.6091,-0.1461,1003,1.261,2.232,986.9,1.266,5.373]
    qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3 = [59.02,0.4802,-0.006005,5.151,2.465,-2.068,2.91,2.919,0.1847]
    return t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3

def fit_0(time_step):
    #fit for 50mm/s
    tinitial = 0
    tfinal = 10
    tstep = time_step
    t = numpy.r_[tinitial:tfinal:tstep]
    qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3 = [0,0,0,0,0,0,0,0,0]
    qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3 = [0,0,0,0,0,0,0,0,0]
    qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3 = [0,0,0,0,0,0,0,0,0]
    return t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3


def fit_10_amount(time_step):
    #fit for 10mm/s
    tinitial = 0
    tfinal = 20
    tstep = (tfinal-tinitial)/ time_step
    t = numpy.r_[tinitial:tfinal:tstep]
    qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3 = [59.39,0.09915,1.046,49.89,0.1134,4.763,4.963,0.2814,6.505]
    qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3 = [407.5,0.007948,0.01161,17,0.2033,-0.2733,1.725,0.4184,6.076]    
    qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3 = [65.4,0.08329,-0.0815,8.104,0.2372,6.337,0.7624,0.5186,4.748]
    return t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3

def fit_20_amount(time_step):
#fit    for 20mm/s
    tinitial = 0
    tfinal = 10
    tstep = (tfinal-tinitial)/ time_step
    t = numpy.r_[tinitial:tfinal:tstep]
    qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3 = [54.05,0.2079,0.8322,37.66,0.2473,4.563,3.95,0.6362,6.234]
    qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3 = [73.64,0.247,-0.175,357.2,0.5185,2.136,345,0.5222,5.283]
    qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3 = [58.86,0.2063,-0.05415,17.15,1.042,-2.215,15.12,1.081,0.7287]
    return t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3

def fit_30_amount(time_step):
    #fit for 30mm/s
    tinitial = 0
    tfinal = 200/30
    tstep = (tfinal-tinitial)/ time_step
    t = numpy.r_[tinitial:tfinal:tstep]
    qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3 = [88.94,0.4628,-0.1873,128.6,0.6458,2.73,78.07,0.7099,5.841]
    qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3 = [64.41,0.251,0.03774,4.598,1.008,-0.6726,0.4017,1.85,-0.3148]
    qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3 = [59.27,0.2861,-0.04419,12,1.432,-1.887,9.104,1.546,0.9067]
    return t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3

def fit_40_amount(time_step):
    #fit for 40mm/s
    tinitial = 0
    tfinal = 5
    tstep = (tfinal-tinitial)/ time_step
    t = numpy.r_[tinitial:tfinal:tstep]
    qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3 = [41.65,0.4353,0.2967,394.5,0.9561,3.474,386.8,0.9674,6.597]
    qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3 = [55.94,0.3319,0.269,76.83,0.9285,3.845,72.52,0.9913,6.798]
    qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3 = [54.53,0.401,0.4147,104.6,0.8511,3.842,88.3,0.9431,6.758]
    return t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3

def fit_50_amount(time_step):
    #fit for 50mm/s
    tinitial = 0
    tfinal = 4
    tstep = (tfinal-tinitial)/ time_step
    t = numpy.r_[tinitial:tfinal:tstep]
    qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3 = [76.13,0.77,-0.1489,184.6,1.185,2.624,145.9,1.248,5.722]
    qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3 = [77.78,0.6091,-0.1461,1003,1.261,2.232,986.9,1.266,5.373]
    qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3 = [59.02,0.4802,-0.006005,5.151,2.465,-2.068,2.91,2.919,0.1847]
    return t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3

def fit_0_amount(time_step):
    #fit for 50mm/s
    tinitial = 0
    tfinal = 10
    tstep = (tfinal-tinitial)/ time_step
    t = numpy.r_[tinitial:tfinal:tstep]
    qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3 = [0,0,0,0,0,0,0,0,0]
    qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3 = [0,0,0,0,0,0,0,0,0]
    qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3 = [0,0,0,0,0,0,0,0,0]
    return t,tinitial,tfinal,tstep,qAa1,qAb1,qAc1,qAa2,qAb2,qAc2,qAa3,qAb3,qAc3,qBa1,qBb1,qBc1,qBa2,qBb2,qBc2,qBa3,qBb3,qBc3,qCa1,qCb1,qCc1,qCa2,qCb2,qCc2,qCa3,qCb3,qCc3


def exp_fit(x,a1, b1, c1,a2,b2,c2,a3,b3,c3):
    from scipy import sin
    return a1*sin(b1*x+c1)+a2*sin(b2*x+c2)+a3*sin(b3*x+c3)

def exp_fit2(x,a1, b1, c1,a2,b2,c2,a3,b3,c3,extra_time):
    from scipy import sin
    fit_1 = a1*sin(b1*x+c1)+a2*sin(b2*x+c2)+a3*sin(b3*x+c3)
    return 