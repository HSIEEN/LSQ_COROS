# Data : 2023/4/07 18:14
# Author: Shawn Shi
# Right Reserved By COROS
import os.path
import sys
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import numpy as np
from numpy import NaN, nan
import pymap3d as wgs
from create_KML import create_kml
from observables import PR_observable_parse

flow_control = False
check_checksum = False


# #################################################################################################
def least_squares(sat_pos, pseudoranges, weights=1, x_hat=np.array([0, 0, 0, 0])):
    """
    Args:
      sat_pos: The satellite position (meters) in an ECEF coordinate frame
      pseudoranges: The corrected pseudorange
      (i.e. a closer approximation to the geometric range from the phone to the satellite)
      weights:
      x_hat: the phone's initial/previous estimated position (x, y, z, b) and
             b represent the user clock bias in units of distance = clock bias (t) * li ght speed (c)

    Returns:
      x_hat: the user's estimated position
      norm_dp: residuals
    """
    dx = np.Inf * np.ones(3)
    G = np.ones((pseudoranges.size, 4))
    iterations = 0
    res_sat_pos = []
    res_pseudorange = []
    res_weights = []

    if isinstance(weights, np.ndarray):
        weights = np.diag(weights)
    else:
        weights = weights * np.eye(pseudoranges.size)
    norm2_dx = np.linalg.norm(dx)
    mean_dx = norm2_dx/len(dx)
    # del_flag = False  # Avoid the consecutive removal operation
    while norm2_dx > 1e-2 and iterations < 30:
        norms = np.linalg.norm(sat_pos - x_hat[:3], axis=1)
        dp = pseudoranges - norms - x_hat[3]
        abs_dp = np.absolute(dp)
        norms_dp = np.linalg.norm(abs_dp)
        mean_dp = norms_dp/len(dx)
        max_dp = np.max(abs_dp)
        G[:, 0:3] = -(sat_pos - x_hat[:3]) / norms[:, None]
        if (mean_dx < 1e2) and (len(sat_pos) > 7) and (max_dp > 25 or mean_dp > 10):
            index = np.argmax(abs_dp)
            res_sat_pos.append(sat_pos[index, :])
            sat_pos = np.delete(sat_pos, index, axis=0)
            res_pseudorange.append(pseudoranges[index])
            pseudoranges = np.delete(pseudoranges, index)
            G = np.delete(G, index, axis=0)
            res_weights.append(weights[index, index])
            weights = np.delete(weights, index, axis=0)
            weights = np.delete(weights, index, axis=1)
            dp = np.delete(dp, index)

        dx = np.linalg.pinv(weights @ G) @ weights @ dp
        norm2_dx = np.linalg.norm(dx)
        mean_dx = norm2_dx/len(dx)
        x_hat = x_hat + dx
        iterations += 1
    res = np.linalg.norm(dp)/len(sat_pos)
    used_sat_num = len(sat_pos)
    res_sat_pos = np.array(res_sat_pos)
    res_weights = np.array(res_weights)
    res_pseudorange = np.array(res_pseudorange)
    # secondary calculation
    if (res > 2) and len(res_sat_pos) > 7:
        # recursive call?
        x_hat, res, iterations, used_sat_num = least_squares(res_sat_pos, res_pseudorange, res_weights)
    return x_hat, res, iterations, used_sat_num


if __name__ == '__main__':
    if_plot = True
    root = tk.Tk()
    root.withdraw()
    print('---------COROS LSQ satellites selection optimization v1.2.0-----------')
    print('-----------All rights are reserved by COROS------------')
    print('Please enter a number to select the frequency band, 1 --> l1, 2 --> l5, 3--> l1 and l5: ')
    band = input()
    if band == '1':
        band = 0
        print('You have selected band L1 data to generate a trajectory')
    elif band == '2':
        band = 1
        print('You have selected band L5 data to generate a trajectory')
    elif band == '3':
        band = 2
        print('You have selected bands l1 and l5 data to generate a trajectory')
    else:
        band = 0  # default band is l1
        print('Default band L1 data is selected to generate a trajectory')
    print('-----------Please choose a raw data file (*.txt)------------')
    file = filedialog.askopenfile()
    if not file:
        sys.exit("******未选择任何文件，自动退出程序！！！******")
    file = file.name
    coord = []
    # TOW_pair = []
    lat0, lon0 = -1, -1
    line_dict = {}
    ##################################################################################################
    if if_plot:
        fig, ax = plt.subplots()
        xs = []
        ys = []
    PR_observable_parse(file)
    with open(file[:-4] + '_PR_observables.txt', 'r') as fPR:
        for pline in fPR:
            # line_dict.clear()
            # sat_pos = []
            # pseudorange = []
            # pseudorange_weight = []
            line_dict = eval(pline)
            TOW = line_dict[list(line_dict.keys())[0]][0]
            if TOW != 527322.0:
                continue
            for key in list(line_dict.keys()):
                # if dictionary data are incomplete or pseudorange is 0 or GNSS type is L5, delete the element
                if len(line_dict[key]) != 6 or line_dict[key][1] == 0 or (key[0] % 2 != band):
                    del line_dict[key]
            # sat_pos = np.array[pseudorange_measurement[]]
            sat_pos = np.array(
                [[line_dict[key][3], line_dict[key][4], line_dict[key][5]]
                 for key in line_dict.keys()])
            pseudorange = np.array([line_dict[key][1] for key in line_dict.keys()])
            pseudorange_weight = \
                1 / np.array([line_dict[key][2] for key in line_dict.keys()])
            # TOW = line_dict[list(line_dict.keys())[0]][0]
            # line_dict.clear()
            # print(sat_pos)
            # print(pseudorange)
            print('##############################################################################')
            if len(sat_pos) > 7:
                # if len(TOW_pair) == 2:
                print(f'On TOW: {TOW}')
                print(f'Initially, {len(sat_pos)} observations are used for positioning.')
                # if TOW == 108492.001:
                x, res, iterations, sat_number = least_squares(sat_pos, pseudorange, pseudorange_weight)
                # else:
                #     continue
                # print(np.linalg.norm(x))
                # x, res = pr_residual(sat_pos, pseudorange, weights=1)
                lat, lon, alt = wgs.ecef2geodetic(x[0], x[1], x[2])
                if if_plot and res < 1e3:
                    if lat0 == -1 and lon0 == -1:
                        lat0, lon0, alt0 = lat, lon, alt
                    e, n, __ = wgs.ecef2enu(x[0], x[1], x[2], lat0, lon0, alt0)
                    ax.cla()
                    xs.append(e)
                    ys.append(n)
                    ax.plot(xs, ys, color='r', linewidth=3)
                    plt.xlabel('East/m')
                    plt.ylabel('North/m')
                    plt.title('LSQ Trajectory without velocity')
                    # ax.xlabel()
                    # ax.legend()
                    plt.pause(0.2)
                # ax.plot(x, y, 'r-', linewidth=2)
                if res < 1e4:
                    coord.extend([str(lon), str(lat), '0'])
                print('The position is', end=':')
                print([lat, lon, alt])
                print('The residual is', end=':')
                print(res)
                print(f'Iterations: {iterations}')
                print(f'Finally, {sat_number} observations are used for positioning')
            else:
                print('Observations are not enough!')
    plt.show(block=True)
    coord = ','.join(coord).replace(',0,', ',0 ')
    if band == 0:
        create_kml(file[:-4] + '_COROS_weighted_LSQ_l1_only', coord)
    else:
        create_kml(file[:-4] + '_COROS_weighted_LSQ_l5_only', coord)
