import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def mcft_solver(s_x, s_y, f_c, agg_size, rho_sy, rho_sx, f_yy, f_yx, f_y, f_x, eps_c=0.002, E_s=200000, shear_area=1000, plot=True):
    """mcft_solver uses the Full MCFT Process to iteratively solve for the maximum shear capacity of a reinforced conrete member
    INPUTS:
    s_x - longitudinal reinforcement spacing (mm)
    s_y - transverse reinforcement spacing (mm)
    f_c - compressive strength of concete (MPa)
    agg_size - diameter of aggregate used in concrete (mm) A value of 6mm is recommended if this is not known
    rho_sy - the ratio of area of transverse steel to total area of a concrete section
    rho_sx - the ratio of area of longitudinal steel to total area of a concrete section
    f_yy - yield strength of transverse steel (MPa)
    f_yx - yield strength of longitidinal steel (MPa)
    f_y - applied vertical axial load (MPa), tension is positive
    f_x - applied horizontal axial load (MPa), tension is positive
    eps_c - Optional: Default value of 0.002, the cracking strain of concrete
    E_s - Optional: Default value of 200,000 (MPa)
    shear area - Optional: Default value of 1000 is equivalent to calculating shear stress (mm^2). The width is to be
    taken as the spalled width (width to outside of the transverse reinforcement). Depth is the total Depth less cover
    plot - Optional: Default value of True will display a plot the predicted member shear strain vs shear stress
    OUTPUTS:
    v - maximum shear stress (MPa) or V - Shear force (kN) if a shear area is provided
    """

    tolerance_f_sy = 1
    tolerance_f_x = 1
    f_c = -f_c
    eps_c = -eps_c
    deltafc1 = 0.001
    E_c = 2 * f_c / eps_c
    f_cr = 0.36 * (-f_c) ** 0.5
    eps_cr = f_cr / E_c

    # step 1
    sm_x = 1.5 * s_x
    sm_y = 1.5 * s_y

    # step 2
    Pi = math.pi

    # step 3
    delta_theta = 0.1

    # step 5
    delta_fsy = 0.1

    v_xy_list = []
    gamma_xy_list =[]
    no_output_count = 0

    for eps_1 in range(200):
        #print(f' iter: {eps_1 + 1:n}/200')

        # resetting vars
        over_counter_limit = False
        fc1_iters = 0
        theta = 1
        f_sy_est = 0.01
        counter = 0
        eps_1 *= 10 ** -4

        run = True
        iter_theta = False
        iter_fsy = False
        f_sy_repeat = 0
        while run:
            counter = counter + 1
            if counter == 10000:
                over_counter_limit = True
                break

            theta *= Pi / 180
            # step 4
            s_theta = 1 / (math.sin(theta) / sm_x + math.cos(theta) / sm_y)
            w = eps_1 * s_theta # maybe multiply by 0.75 as in vicroads paper

            # step 6
            if eps_1 <= eps_cr:  #this condition is never true
                f_c1 = E_c * eps_1
            else:
                f_c1 = f_cr / (1 + (500 * eps_1) ** 0.5) # 200 changed to 500

            k = max(1.64 - 1 / math.tan(theta), 0)
            v_ci_max = (-f_c) ** 0.5 / (0.31 + 24 * w / (agg_size + 16))
            if __name__ == '__main__':
                f_c1 = min(v_ci_max * (0.18 + 0.3 * k ** 2) * math.tan(theta) + rho_sy * (f_yy - f_sy_est), f_c1) \
                    - fc1_iters * deltafc1

            # step 7
            f_cy = f_y - rho_sy * f_sy_est
            v_xy = (f_c1 - f_cy) / math.tan(theta)

            # step 8
            f_c2 = f_c1 - v_xy * (math.tan(theta) + 1 / math.tan(theta))  # assuming v_xy = v_cxy due to theta = theta_c
            # step 9
            f_c2_max = f_c / (0.8 - 0.34 * eps_1 / eps_c)

            f_c2_max = max(f_c, f_c2_max)

            # step 10
            if f_c2 > f_c2_max:   # iterate theta if false

                # step 11
                eps_2 = eps_c * (1 - (1 - f_c2 / f_c2_max) ** 0.5)

                # step 12
                eps_y = (eps_1 + eps_2 * math.tan(theta) ** 2) / (1 + math.tan(theta) ** 2)

                # step 13
                f_sy = min(E_s * eps_y, f_yy)

                # step 14
                if abs(f_sy - f_sy_est) <= tolerance_f_sy:  # iterate f_sy_est if false

                    # step 15
                    eps_x = eps_1 + eps_2 - eps_y

                    # step 16
                    f_sx = min(E_s * eps_x, f_yx)

                    # step 17
                    f_cx = f_c1 - v_xy / math.tan(theta)
                    calc_f_x = f_cx + rho_sx * f_sx

                    # step 18
                    if abs(calc_f_x - f_x) <= tolerance_f_x:  # iterate theta if false
                        delta_f_c1 = f_c1 - rho_sy * (f_yy - f_sy)
                        # step 19
                        if delta_f_c1 <= 0:
                            v_ci = 0
                            f_ci = 0
                        else:
                            c = delta_f_c1 / math.tan(theta) - 0.18 * v_ci_max
                            if c <= 0:
                                f_ci = 0
                                v_ci = delta_f_c1 / math.tan(theta)
                            else:
                                a = 0.82 / v_ci_max
                                b = 1 / math.tan(theta) - 1.64
                                f_ci = (-b - (b ** 2 - 4 * a * c) ** 0.5) / (2 * a)
                                v_ci = (f_c1 + delta_f_c1) / math.tan(theta)
                        f_sx_cr = f_sx + (f_c1 + f_ci + v_ci / math.tan(theta)) / rho_sx
                        # step 21
                        if f_sx_cr <= f_yx:  # iterate back to step 7 if false (assume lower fc1):
                            # step 22
                            gamma_xy = 2 * (eps_x - eps_2) / math.tan(theta)
                            run = False
                        else:
                            fc1_iters += 1

                    else:
                        iter_theta = True
                else:
                    iter_fsy = True
            else:
                iter_theta = True

            if iter_theta:
                theta = theta * 180 / Pi + delta_theta

            else:
                theta = theta * 180 / Pi

            if iter_fsy:
                f_sy_repeat = f_sy_repeat + 1
                if f_sy_repeat > 50:
                    f_sy_est = f_sy
                    f_sy_repeat = 0
                if f_sy_est > f_sy:
                    f_sy_est = f_sy_est - 2 * delta_fsy

                f_sy_est = f_sy_est + delta_fsy

            iter_theta = False
            iter_fsy = False

        if theta <= 75 and not over_counter_limit:
            v_xy_list.append(v_xy*shear_area/1000)
            gamma_xy_list.append(gamma_xy)
            no_output_count = 0
        else:
            no_output_count += 1

        if no_output_count >= 30:

            break

    if plot:
        plt.plot([x * 10 ** 3 for x in gamma_xy_list], v_xy_list, "r-", label="Predicted Response")
        plt.xlabel(r"$\gamma_{xy} \times 10 ^{3} $")
        plt.ylabel("$v_{xy} $ (MPa)")
        scale_val = 1
        if shear_area != 1000:
            scale_val = 100
            plt.ylabel("$V_{xy} $ (kN)")
        plt.title(f"Python Implementation of Full MCFT")
        plt.xlim(-0.5, 20)
        plt.ylim(0, math.ceil(max(v_xy_list, default=np.nan) + scale_val))
        plt.show()

    return max(v_xy_list, default=np.nan)

