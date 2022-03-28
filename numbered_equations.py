# Note, all numbered equations are from (Bentz, Vecchio, and Collins 2006).

import math

# EQUILIBRIUM:


# Average Stresses:
# list of input vars: rho_x, rho_z, f_sx, f_sx, f_1, f_2, v, theta
#        output vars: f_x, f_z, v

def eq_1(rho_x,f_sx, f_1, v, theta):
    # returns f_x
    return rho_x * f_sx + f_1 - v / math.tan(theta)


def eq_2(rho_z,f_sz, f_1, v, theta):
    # returns f_z
    return rho_z * f_sz + f_1 - v * math.tan(theta)


def eq_3(f_1, f_2, theta):
    # returns v
    return (f_1 + f_2)/(math.tan(theta) + 1 / math.tan(theta))


# Stress at Cracks:
# list of input vars: f_x, f_z, v, theta, v_ci, rho_x, rho_z
#        output vars: f_sxcr, f_szcr

def eq_4(f_x, v, theta, v_ci, rho_x):
    # returns f_sxcr
    return (f_x + (v + v_ci) / math.tan(theta)) / rho_x


def eq_5(f_z, v, theta, v_ci, rho_z):
    # returns f_szcr
    return (f_z + (v - v_ci) * math.tan(theta)) / rho_z


# Geometric Conditions


# Average Strains:
# list of input vars: eps_x, eps_z, eps_2
#        output vars: theta, eps_1, gamma_xz

def eq_6(eps_x, eps_z, eps_2):
    # returns tan^2 theta
    return (eps_x + eps_2) / (eps_z + eps_2)


def eq_7(eps_x, eps_z, eps_2):
    # returns eps_1
    return eps_x + eps_z + eps_2


def eq_8(eps_x, eps_2, theta):
    # returns gamma_xz
    return 2 / math.tan(theta) * (eps_x + eps_2)


# Stress- Strain Relationships


# Crack Widths
# list of input vars: rho_x, rho_z, f_sx, f_sx, f_1, f_2, v, theta
#        output vars: f_x, f_z, v

def eq_9(s_theta, eps_1):
    # returns w
    return s_theta * eps_1


def eq_10(theta, s_x, s_z):
    # returns s_theta
    return 1 / (math.sin(theta)/s_x + math.cos(theta)/s_z)


# Stress- Strain relationships


# Reinforcement
# list of input vars: rho_x, rho_z, f_sx, f_sx, f_1, f_2, v, theta
#        output vars: f_x, f_z, v

def eq_11(E_s, eps_x, f_yx):
    # returns f_sx
    return min( E_s * eps_x, f_yx)


def eq_12(E_s, eps_z, f_yz):
    # returns f_sz
    return min( E_s * eps_z, f_yz)


# Concrete
# list of input vars: rho_x, rho_z, f_sx, f_sx, f_1, f_2, v, theta
#        output vars: f_x, f_z, v

def eq_13(f_c, eps_1, eps_2, eps_c):
    # returns f_2
    return f_c/(0.8 + 170 * eps_1)*(2 * eps_2 / eps_c - (eps_2 / eps_c) ** 2)


def eq_14(f_c, eps_1):
    # returns f_1
    return 0.33 * f_c ** 0.5 / (1 + (500 * eps_1) ** 0.5)


def eq_15(f_c, w, a_g):
    # returns v_ci upper limit
    return 0.18 * f_c ** 0.5 / (0.31 + 24 * w / (a_g + 16))
