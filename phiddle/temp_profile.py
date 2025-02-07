import numpy as np
from util import two_lorentz, oned_gaussian_func

#Version 1.2, used since Feb 7 2024 at 21:24

temperature_profile_func_dict = {
        "2025": lambda dwell, tpeak: two_lorentz(tpeak, 0., *left_right_width_2025(dwell, tpeak)),
        "2024": lambda dwell, tpeak: two_lorentz(tpeak, 0., *left_right_width_2024(dwell, tpeak)), 
        "2023": lambda dwell, tpeak: two_lorentz(tpeak, 0., *left_right_width_2023(dwell, tpeak)),
        "2021": lambda dwell, tpeak: oned_gaussian_func(tpeak, 0.,
                                      sigma_Fall2021(dwell, LaserPowerMing_Fall2021(dwell, tpeak)))
        }

width_func_dict = {
        "2025": lambda dwell, tpeak: left_right_width_2025(dwell, tpeak), 
        "2024": lambda dwell, tpeak: left_right_width_2024(dwell, tpeak), 
        "2023": lambda dwell, tpeak: left_right_width_2023(dwell, tpeak),
        "2021": lambda dwell, tpeak: sigma_Fall2021(dwell, LaserPowerMing_Fall2021(dwell, tpeak))
        }

# Surface functions
def temp_surface_func(b, c, d, e, f):
    return lambda x, y: (b*x+c)*(y)**(d*x**2+e*x+f) + 27

def inverse_temp_surface_func(b, c, d, e, f):
    return lambda dw, t: ((t-27)/(b*dw+c))**(1/(d*dw**2 + e*dw +f))

def twod_surface(base, a, b, c, d, e):
    return lambda x, y: base + a*x + b*y + c*x**2 + d*y**2 + e*x*y

def cubic_surface(base, a, b, c, d, e, f, g):
    return lambda x, y: base + a*x + b*y + c*x**2 + d*y**2 + e*x*y + f*x**4 + g*y**4

############## Year 2025 ##################
def left_right_width_2025(tau, Temp):
    left_width_fit =  [ 5.02794524e+02,
                       2.27809447e+02,
                       -1.33599129e+01,
                       -8.12411230e+01,
                       1.76867367e-01,
                       -1.65872055e+00,
                       4.51417474e+00,
                       -5.67899883e-06]
    right_width_fit = [ 9.00132242e+02,
                       -1.60851951e+02,
                       -9.78194951e+00,
                       3.82337277e+01,
                       -4.75956003e-02,
                       2.49843621e+00]

    left_width = cubic_surface(*left_width_fit) 
    right_width = twod_surface(*right_width_fit)
    log10_velocity = np.log10(88200./tau)
    power = LaserPowerMing_Spring2025(tau, Temp, temp_fit = None)
    return left_width(log10_velocity, power), right_width(log10_velocity, power)


def LaserPowerMing_Spring2025(dwell, Tpeak, temp_fit = None):
    print("USING POWER PROFILE MING SPRING 2024")
    if temp_fit is None:
        temp_fit = [-0.01824834, 0.05924233, -0.01709909, 0.01145558, 2.74505353]

    velo = 88200/dwell
    log10vel = np.log10(velo)
    get_power = inverse_temp_surface_func(*temp_fit) 
    power = get_power(log10vel, Tpeak)
    return power

########### Year 2023 #################
def left_right_width_2023(tau, Temp):
    """
    Temperature profile for Spring 2023 CHESS run
    The temperature is asymmetric so left and right width of the temperature
    distribution is fitted separately.

    Args:
        dwell: float = Dwell time of the anneal
        Tpeak: float = Peak temperature of the anneal
    Return:
        left width: float = left width parameter
        right width: float = right width parameter
    """
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
    """
    Laser power for Spring 2023 CHESS run
    Given the required dweel and Tpeak, return the power that gives that Tpeak
    Args:
        dwell: float = Dwell time of the anneal
        Tpeak: float = Peak temperature of the anneal
    Return:
        power: float = Power that reaches that peak temperautre given the dwell time
    """
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
    """
    Temperature profile for Spring 2024 CHESS run
    The temperature is asymmetric so left and right width of the temperature
    distribution is fitted separately.

    Args:
        dwell: float = Dwell time of the anneal
        Tpeak: float = Peak temperature of the anneal
    Return:
        left width: float = left width parameter
        right width: float = right width parameter
    """
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
    """
    Laser power for Spring 2024 CHESS run
    Given the required dweel and Tpeak, return the power that gives that Tpeak
    Args:
        dwell: float = Dwell time of the anneal
        Tpeak: float = Peak temperature of the anneal
    Return:
        power: float = Power that reaches that peak temperautre given the dwell time
    """
    if temp_fit is None:
        temp_fit = [-1.71420431e-04,  6.37694255e-04, -8.72721878e-02,  2.21288491e-01, 3.64085959e+00]

    velo = 88200/dwell
    log10vel = np.log10(velo)
    get_power = inverse_temp_surface_func(*temp_fit) 
    power = get_power(log10vel, Tpeak)
    return power

########### End of Year 2024 #################

########### Year 2021 #############
def LaserPowerMing_Fall2021(tau, Temp):
    """
    Laser power for Spring 2021 CHESS run
    Given the required dweel and Tpeak, return the power that gives that Tpeak
    Args:
        dwell: float = Dwell time of the anneal
        Tpeak: float = Peak temperature of the anneal
    Return:
        power: float = Power that reaches that peak temperautre given the dwell time
    """
    tpeak = Temp
    dwell = np.log10(tau)
    n_sqrt = 0.0190715357103267*dwell**4 - 0.180649351400508*dwell**3 + 0.651605603621984*dwell**2 + 5.21414120166499e-5*dwell*tpeak - dwell - 0.000206391126900198*tpeak + 0.569129395966527
    if n_sqrt >= 0.:
        P = (-6.40883896736885e+15*dwell**2 + 3.78789311741561e+16*dwell + 4.64073240375504e+16*np.sqrt(n_sqrt) - 7.48733837342515e+16)/(112293816201517.0*dwell - 444492129640718.0)
    else:
        print("Requested power cannot be reached!")
        P = 1.e10
    if np.isnan(P) or P < 0.:
        print("Requested power cannot be reached! Exiting now!")
        P = 1.e10
    return P

def sigma_Fall2021(tau, current):
    """ 
    Will return the sigma given a dwell in log10(tau) and current in Amps
    """ 
    dwell = np.log10(tau)
    return  [(6.70640332e+02 
            - 2.22189018e+02*dwell
            - 2.27743561e+00*current
            + 3.02469331e+01*dwell**2
            + 7.70978201e-03*current**2
            +1.70903683e-01*dwell*current)]

########## End of Year 2021 #########

if __name__ == '__main__':
    left_width_fit = [ 9.23610385e+02,
                       2.31591101e+02,
                      -2.33724136e+01,
                      -6.63593631e+01,
                       2.38094452e-01,
                      -2.37653526e+00, 
                       5.05393010e+00,
                      -4.32250685e-06]

    right_width_fit = [ 4.39339896e+02,
                       -4.51657906e+01,
                       -9.35886967e+00,
                       -3.51825565e+00,
                       -2.20075111e-02,
                       2.44098710e+00]

    left_width = cubic_surface(*left_width_fit) 
    right_width = twod_surface(*right_width_fit)
    print("Latest power", LaserPowerMing_Spring2024(10000, 1000))
    left_width, right_width = left_right_width_2024(10000, 1000)
    print("Left width", left_width)
    print("Right width", right_width)
