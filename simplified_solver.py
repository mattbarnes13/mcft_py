import math
import numpy as np
from scipy.optimize import fsolve
import pandas as pd
from matplotlib import pyplot as plt


def simplified_solver(s_x, f_c, agg_size, rho_sy, rho_sx, f_yy, f_yx, f_x, E_s=200000):
    """simplified_solver uses the Simplified MCFT Process to iteratively solve for the maximum shear capacity of a reinforced conrete member
    INPUTS:
    s_x - longitudinal reinforcement spacing (mm)
    f_c - compressive strength of concete (MPa)
    agg_size - diameter of aggregate used in concrete (mm) A value of 6mm is recommended if this is not known
    rho_sy - the ratio of area of transverse steel to total area of a concrete section
    rho_sx - the ratio of area of longitudinal steel to total area of a concrete section
    f_yy - yield strength of transverse steel (MPa)
    f_yx - yield strength of longitidinal steel (MPa)
    f_x - applied horizontal axial load (MPa), tension is positive
    E_s - Optional: Default value of 200,000 (MPa)
    OUTPUTS:
    v - maximum shear stress (MPa)
    """
    tolerance = 10 ** -5
    z = np.array([math.pi / 6, 0.001])
    iters = 0
    eps_x = 0.001
    increment = 10 ** (-5)
    f_sxcr_check = False
    iter_limit = 100000
    run = True
    if f_c > 70:
        agg_size = 0

    s_xe = 35 * s_x / (agg_size + 16)
    while run is True and iters < iter_limit:

        # theta = math.atan((0.568 + 12.58 * s_xe * eps_1 / math.sin(theta)) / (1 + math.sqrt(500 * eps_1)))
        #
        # eps_1 = eps_x * (1 + 1 / math.tan(theta) ** 2) + 1 / math.tan(theta) ** 4 / (15000 * (1 + math.sqrt(500 * eps_1)))
        if iters >= 1:
            eps_x = next_eps_x

        iters += 1

        def my_function(z):
            theta = z[0]
            eps_1 = z[1]
            F = np.empty(2)
            F[0] = math.tan(theta) - (0.568 + 1.258 * s_xe * eps_1 / math.sin(theta)) / (1 + math.sqrt(500 * eps_1))
            F[1] = eps_1 - (eps_x * (1 + 1 / math.tan(theta) ** 2) + 1 / math.tan(theta) ** 4 / (15000 * (1 + math.sqrt(500 * eps_1))))
            F[2] = f_x - rho_sx * min(eps_x * eps_y , f_yy) + 0.33*math.sqrt(f_c) / (1 + math.sqrt(500*eps_1)) - v/math.tan(theta)
            return F

        try:
            z = fsolve(my_function, z)
            theta = min(z[0], 75 * math.pi / 180)
            eps_1 = z[1]

            beta = 0.18 / (0.31 + 0.686 * s_xe * eps_1 / math.sin(theta))

        except Exception:
            beta = 0.4/(1+1500*eps_x) * 1300/(1000+s_xe)
            theta = min((29+7000*eps_x) * (0.88 + s_xe/2500), 75)
            theta = theta * math.pi/180

        v_c = beta * math.sqrt(f_c)
        v_s = rho_sy * f_yy * 1 / math.tan(theta)
        v = v_c + v_s

        if eps_x <= 1.9 * 10 ** (-3):
            eps_x_calc = (v / math.tan(theta) - v_c * math.tan(theta)) / (E_s * rho_sx)
            if abs(eps_x - eps_x_calc) < tolerance:
                run = False

        next_eps_x = (eps_x + eps_x_calc)/2

        f_sxcr = (v + v_c) / math.tan(theta) / rho_sx + f_x / v * v / rho_sx

        if f_sxcr > f_yx:
            if f_sxcr_check is False:
                eps_x = 3 * 10 ** (-3)
            f_sxcr_check = True
            next_eps_x = eps_x + increment



        if f_sxcr_check is True and f_sxcr <= f_yx:
            run = False

        if rho_sy * f_yy /f_c > 0.25:
            v = 0.25 * f_c
            break

        f_x = min(rho_sx*E_s*eps_x - v/math.tan(theta) + v_c * math.tan(theta), rho_sx*f_yx - (v + v_c)/math.tan(theta))

    if iters >= iter_limit:
        v = np.nan

    return v
