import numpy as np
import math
import matplotlib.pyplot as plt
import pandas


def cross_sec_slice(beam_type, beam_dimensions, trans_reo,  num_slices=20):
    """cross_sec_slice takes an input beam_type, beam dimensions and transverse reinforcement to produce a list
    of slices that can be used for more advanced MCFT calculations
    INPUTS:
    beam_type - values are ibeam, tbeam or rectangle
    beam dimensions - as a list
        for an ibeam
            [b_bf, t_bf, b_tf, t_tf, b_w, D] where all values are in mm and D is the total depth of the section
        for a tbeam
            [b_f, t_f, b_w, D] where all values are in mm and D is the total depth of the section
        for a rectangle
            [b, d] where all values are in mm
    trans_reo - as a list
        [b1, diam_sy, s_y] where b1 is the spalled width (mm), diam_sy is the diameter of transverse reinforcement (mm)
        s_y is the spacing of transverse reinforcement (mm)
    num_slices - Optional: default is 20 slices

    OUTPUTS:
    slice_list - list of length 'num_slices' containing [b1, h_i, rho_syi, y_ci] for each slice, where
        b1 - is the spalled width (mm) of the slice,
        h_i - is the thickness or height of the slice (mm)
        rho_syi - is the ratio of transverse reinforcement
        y_ci - is the distance from slice centroid to the top of the section (mm)

    I -  the second moment of inertia (mm^4)
    y_c - the height of section centroid from the bottom of the section (mm)
    """
    # Input Validation
    # beam type
    if beam_type not in ["ibeam", "tbeam", "rectangle"]:
        raise ValueError(f'"{beam_type}" is not one of the valid beam types (ibeam, tbeam or rectangle)')

    s_y = trans_reo[2]
    b1 = trans_reo[0]
    diam_sy = trans_reo[1]
    a_sy = diam_sy ** 2 * math.pi / 4

    # beam dimensions
    if beam_type == "ibeam":
        if len(beam_dimensions) != 6 \
                or (beam_dimensions[1] + beam_dimensions[3]) >= beam_dimensions[5] \
                or min(beam_dimensions[0], beam_dimensions[2]) <= beam_dimensions[4] \
                or b1 >= beam_dimensions[4]:
            raise ValueError(f'beam_dimensions are not valid!')

    if beam_type == "tbeam":
        if len(beam_dimensions) != 4 \
                or beam_dimensions[1] >= beam_dimensions[3] \
                or beam_dimensions[0] <= beam_dimensions[2]\
                or b1 >= beam_dimensions[2]:
            raise ValueError(f'beam_dimensions are not valid!')

    if beam_type == "rectangle":
        if len(beam_dimensions) != 2 or beam_dimensions[0] < b1:
            raise ValueError(f'beam_dimensions are not valid!')

    slice_dims = []

    if beam_type == "ibeam":
        b_bf = beam_dimensions[0]
        t_bf = beam_dimensions[1]
        b_tf = beam_dimensions[2]
        t_tf = beam_dimensions[3]
        b_w = beam_dimensions[4]
        D = beam_dimensions[5]
        A_tf= b_tf * t_tf
        A_bf = b_bf * t_bf
        A_web = (D - t_tf - t_bf) * b_w
        A_t = A_tf + A_bf + A_web
        A_approx = A_t / num_slices

        n1 = round(A_tf / A_approx)
        A_1 = A_tf / n1
        h1 = A_1/ b_bf
        n2 = round(A_web / A_approx)
        A_2 = A_web / n2
        h2 = A_2/ b_w
        n3 = round(A_bf / A_approx)
        A_3 = A_bf / n3
        h3 = A_3/ b_tf
        t_w = h2*n2
        y_c = (t_tf * b_tf * (t_tf / 2 + t_w+t_bf) + t_w * b_w * (t_w / 2+t_bf) + t_bf*b_bf * (t_bf/2)) / A_t
        print(t_w)
        I = A_web * t_w ** 2 / 12 + A_web * (y_c - t_w / 2 -t_bf) ** 2 + A_tf * t_tf ** 2 / 12 + A_tf * (y_c - t_tf / 2 - t_w- t_bf ) ** 2 \
            + A_bf * t_bf ** 2 / 12 + A_bf * (y_c - t_bf /2) ** 2
        slice_dims = [[A_1,h1,n1],[A_2,h2,n2],[A_3,h3,n3]]

    if beam_type == "tbeam":
        b_f = beam_dimensions[0]
        t_f = beam_dimensions[1]
        b_w = beam_dimensions[2]
        D = beam_dimensions[3]
        A_f = b_f * t_f
        A_web = (D - t_f) * b_w
        A_t = A_f + A_web
        A_approx = A_t / num_slices

        n1 = round(A_f / A_approx)
        A_1 = A_f / n1
        h1 = A_1 / b_f
        n2 = round(A_web / A_approx)
        A_2 = A_web / n2
        h2 = A_2 / b_w
        t_w = h2*n2
        y_c = (t_f*b_f * (t_f/2+t_w) + t_w*b_w*(t_w/2))/A_t
        I = A_web * t_w ** 2 /12 + A_web * (y_c - t_w/2) ** 2  + A_f * t_f ** 2 /12 + A_f * (y_c -t_f/2-t_w)**2
        slice_dims = [[A_1, h1, n1], [A_2, h2, n2]]

    if beam_type == "rectangle":
        b = beam_dimensions[0]
        d = beam_dimensions[1]
        A_t = b*d
        A = A_t / num_slices
        h1 = A/b
        total_slices = num_slices
        y_c = d/2
        I = b*d ** 3 /12
        slice_dims = [[A, h1, total_slices]]

    current_height = 0
    output = []

    for area_section in slice_dims:
        a = area_section[0]
        h_i = area_section[1]
        n = area_section[2]
        b_i = a/h_i
        rho_syi = a_sy / (b_i * s_y )
        for i in range(n):
            y_ci = current_height + h_i / 2
            current_height += h_i
            output.append([b1, h_i, rho_syi, y_ci])

    return output, I, y_c

