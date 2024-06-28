import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import leastsq

#Version 1.2, used since Feb 7 2024 at 21:24

def temp_surface_func(b, c, d, e, f):
    return lambda x, y: (b*x+c)*(y)**(d*x**2+e*x+f) + 27

def inverse_temp_surface_func(b, c, d, e, f):
    return lambda dw, t: ((t-27)/(b*dw+c))**(1/(d*dw**2 + e*dw +f))

def twod_surface(base, a, b, c, d, e):
    return lambda x, y: base + a*x + b*y + c*x**2 + d*y**2 + e*x*y

def cubic_surface(base, a, b, c, d, e, f, g):
    return lambda x, y: base + a*x + b*y + c*x**2 + d*y**2 + e*x*y + f*x**4 + g*y**4

########### Year 2023 #################
def left_right_width_2023(tau, Temp):
    left_width_fit = [ 7.93197988e+02,
                   2.34346130e+01,
                   -4.29544003e+01,
                   4.68763599e+01,
                   9.64778097e-01,
                   -4.95937074e+00,
                   -1.68873411e+00,
                   -9.11029166e-05]
    right_width_fit = [ 4.39339896e+02,
                       -4.51657906e+01,
                       -9.35886967e+00,
                       -3.51825565e+00,
                       -2.20075111e-02,
                       2.44098710e+00]
    left_width = cubic_surface(*left_width_fit) 
    right_width = twod_surface(*right_width_fit)
    log10_velocity = np.log10(88200./tau)
    power = LaserPowerMing_Spring2023(tau, Temp, temp_fit = None)
    return left_width(log10_velocity, power), right_width(log10_velocity, power)

def LaserPowerMing_Spring2023(dwell, Tpeak, temp_fit = None):
    if temp_fit is None:
        temp_fit = [-0.06688143,
                     0.22353149,
                     -0.0556677,
                     0.20606041,
                     2.43702215] 

    velo = 88200/dwell
    log10vel = np.log10(velo)
    get_power = inverse_temp_surface_func(*temp_fit) 
    power = get_power(log10vel, Tpeak)
    return power
########### End of Year 2023 #################

########### Year 2024 #################

def left_right_width_2024(tau, Temp):
    left_width_fit = [ 9.23610385e+02,
                       2.31591101e+02,
                      -2.33724136e+01,
                      -6.63593631e+01,
                       2.38094452e-01,
                      -2.37653526e+00, 
                       5.05393010e+00,
                      -4.32250685e-06]
    right_width_fit =  [ 7.38857381e+02,
                        -2.47679319e+02,
                        -4.12579362e+00,
                        -2.64455714e+01,
                        -1.55386703e-01,
                         7.51253690e+00]


    left_width = cubic_surface(*left_width_fit) 
    right_width = twod_surface(*right_width_fit)
    log10_velocity = np.log10(88200./tau)
    power = LaserPowerMing_Spring2024(tau, Temp, temp_fit = None)
    return left_width(log10_velocity, power), right_width(log10_velocity, power)


def LaserPowerMing_Spring2024(dwell, Tpeak, temp_fit = None):
    if temp_fit is None:
        temp_fit = [-1.71420431e-04,  6.37694255e-04, -8.72721878e-02,  2.21288491e-01, 3.64085959e+00]

    velo = 88200/dwell
    log10vel = np.log10(velo)
    get_power = inverse_temp_surface_func(*temp_fit) 
    power = get_power(log10vel, Tpeak)
    return power

########### End of Year 2024 #################

if __name__ == '__main__':
    pass
    # print("Latest power", LaserPowerMing_Spring2024(10000, 1000))
    # left_width, right_width = left_right_width(10000, 1000)
    # print("Left width", left_width)
    # print("Right width", right_width)
