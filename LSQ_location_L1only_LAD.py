# Data : 2023/3/23 18:14
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
import linecache
from create_KML import create_kml
# from track_parser import track_parser, wgs_enu2lla
from raw_parser import Measurement_raw_parse
from raw_parser import SV_raw_parse

# import linecache as lc

flow_control = False
check_checksum = False

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

    if isinstance(weights, np.ndarray):
        weights = np.diag(weights)
    else:
        weights = weights * np.eye(pseudoranges.size)
    increment_x = np.linalg.norm(dx)
    while increment_x > 1e-4:
        # dp = pseudoranges - norms - x_hat[3]
        norms = np.linalg.norm(sat_pos - x_hat[:3], axis=1)
        dp = pseudoranges - norms - x_hat[3]
        # norm_dp = np.linalg.norm(dp)
        abs_dp = np.absolute(dp)
        max_dp = np.max(abs_dp)
        # mean = np.mean(abs_dp)
        # var = np.var(dp)
        G[:, 0:3] = -(sat_pos - x_hat[:3]) / norms[:, None]
        # Outlier detect and removal
        if (increment_x < 1e2) and (len(sat_pos) > 10) and (max_dp > 30):
            index = np.argmax(abs_dp)
            sat_pos = np.delete(sat_pos, index, axis=0)
            pseudoranges = np.delete(pseudoranges, index)
            # decrease a obs
            G = np.delete(G, index, axis=0)
            weights = np.delete(weights, index, axis=0)
            weights = np.delete(weights, index, axis=1)
            # wights = np.delete(weights, index, axis=1)
            dp = np.delete(dp, index)
        # G_T = np.transpose(G)
        # dx = np.linalg.inv(G_T@G) @ G_T @ dp

        dx = np.linalg.pinv(weights @ G) @ weights @ dp
        increment_x = np.linalg.norm(dx)
        x_hat = x_hat + dx
        iterations += 1
    # res = np.linalg.norm(dp)
    # if res > 2e2:
    #     pass
    return x_hat, np.linalg.norm(dp), iterations, len(sat_pos)


def least_absolutes(sat_pos, pseudoranges, weights=1, x_hat=np.array([0, 0, 0, 1e6])):
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
    #
    if isinstance(weights, np.ndarray):
        weights = np.diag(weights)
    else:
        weights = weights * np.eye(pseudoranges.size)
    increment_x = np.linalg.norm(dx)
    while increment_x > 1e-3:
        # dp = pseudoranges - norms - x_hat[3]
        norms = np.linalg.norm(sat_pos - x_hat[:3], axis=1)
        dp = pseudoranges - norms - x_hat[3]
        # i_dp = 1/dp
        # if np.linalg.norm(1/dp) > 1e-3:
        abs_dp = np.absolute(dp)
        abs_dp[abs_dp == 0] = (np.max(abs_dp)+np.max(abs_dp))/2
        weights = (np.diag(np.sqrt(1 / abs_dp)))
        G[:, 0:3] = -(sat_pos - x_hat[:3]) / norms[:, None]
        # Outlier detect and removal
        # dx = np.linalg.inv(G_T@G) @ G_T @ dp
        dx = np.linalg.pinv(weights @ G) @ weights @ dp
        x_hat = x_hat + dx
        iterations += 1
        increment_x = np.linalg.norm(dx)
    # res = np.linalg.norm(dp)
    # if res > 2e2:
    #     pass
    return x_hat, np.linalg.norm(dp), iterations, len(sat_pos)


# ################################################################################################

# Create figure for plotting
# This function is called periodically from FuncAnimation

#############################################################################################
if __name__ == '__main__':
    if_plot = True
    root = tk.Tk()
    root.withdraw()
    print('---------COROS LSQ solving v1.0.0-----------')
    print('-----------All rights are reserved by COROS------------')
    print('-----------Please choose a raw data file (*.txt)------------')
    file = filedialog.askopenfile()
    if not file:
        sys.exit("******未选择任何文件，自动退出程序！！！******")
    file = file.name
    # if not os.path.exists(file[:-4] + '_Parsed Measurement raw data.txt'):
    #     fm = open(file[:-4] + '_Parsed Measurement raw data.txt', 'w')
    # if not os.path.exists(file[:-4] + '_Parsed SV raw data.txt'):
    #     fs = open(file[:-4] + '_Parsed SV raw data.txt', 'w')
    # raw data information abstract
    #############################################################################################
    with open(file, 'rb') as fo:
        data = fo.read()
        # n = n+1
        parsed_measurement_file = file[:-4] + '_Parsed Measurement raw data.txt'
        parsed_sv_file = file[:-4] + '_Parsed SV raw data.txt'
        if not os.path.exists(parsed_measurement_file):
            with open(parsed_measurement_file, 'w') as fm:
                for it in Measurement_raw_parse(data, flow_control, check_checksum):
                    fm.writelines(str(it) + "\n")
        if not os.path.exists(parsed_sv_file):
            with open(parsed_sv_file, 'w') as fs:
                for it in SV_raw_parse(data, flow_control, check_checksum):
                    fs.writelines(str(it) + "\n")

    pseudorange_measurement = {}
    coord = []
    # get total number of cycles in measurement data
    mCycle = sCycle = cIOD = 0
    line_number = 1
    mdict_cycle = {0: 1}
    sdict_cycle = {0: 1}
    with open(file[:-4] + '_Parsed Measurement raw data.txt', 'r') as fmr:
        for mline in fmr:
            line_dict = eval(mline)
            # if the current IOD is less than the previous IOD, cycle increases by 1
            if line_dict['IOD'] < cIOD:
                mCycle += 1
                mdict_cycle[mCycle] = line_number
            cIOD = line_dict['IOD']
            line_number += 1
        mdict_cycle[mCycle + 1] = line_number - 1
    line_number = 1
    cIOD = 0
    # get total number of cycles in SV information
    ln = 0
    with open(file[:-4] + '_Parsed SV raw data.txt', 'r') as fsr:
        for sline in fsr:

            line_dict = eval(sline)
            ln += 1

            #  assign the line number to sdict_cycle dictionary
            if line_dict['IOD'] < cIOD:
                sCycle += 1
                sdict_cycle[sCycle] = line_number
            cIOD = line_dict['IOD']
            line_number += 1
        sdict_cycle[sCycle + 1] = line_number - 1
    if len(mdict_cycle) != len(sdict_cycle):
        sys.exit("******原始观测量与卫星信息不匹配！！！！******")
    m_line_number = 1
    s_line_number = 1
    m = len(mdict_cycle)
    # TOW_pair = []
    lat0, lon0 = -1, -1
    ##################################################################################################
    if if_plot:
        fig, ax = plt.subplots()
        xs = []
        ys = []
    with open(file[:-4] + '_Parsed Measurement raw data.txt', 'r') as fmr:
        for mline in fmr:
            # get current number of cycle
            for i in range(0, m - 1):
                if mdict_cycle[i] <= m_line_number < mdict_cycle[i + 1]:
                    mCycle = i
                    break
            line_dict = eval(mline)
            m_line_number = m_line_number + 1
            # print(line_dict['IOD'])
            mIOD = line_dict['IOD']
            # for debugging
            if line_dict['Receiver TOW'] != 108429.001:
                continue
            if line_dict['Clock status'] >= 2 and line_dict['Number of observation'] > 0:
                m_keys = list(line_dict.keys())
                for n in range(line_dict['Number of observation']):
                    # GNSS and signal type, select L1 only
                    if line_dict[m_keys[7 + n * 12]] % 2 == 0:
                        # tail = key.split('_')[-1]
                        # 8+n*12 --> PRN code, 7+n*12 --> GNSS and signal type
                        RM_key = (line_dict[m_keys[7 + n * 12]], line_dict[m_keys[8 + n * 12]])
                        # 12+i*12 --> raw pseudorange, 14+i*14 --> pseudorange variance index
                        if line_dict[m_keys[14 + n * 12]] <= 32:
                            pseudorange_measurement[RM_key] = \
                                [line_dict[m_keys[12 + n * 12]], line_dict[m_keys[14 + n * 12]]]
                            TOW = line_dict['Receiver TOW']

            else:
                continue
            RM_keys = [key for key in pseudorange_measurement]
            # with open(file[:-4] + '_Parsed SV raw data.txt', 'r') as fsr:
            # Specify the lines to be searched
            specified_lines = [i for i in range(sdict_cycle[mCycle], sdict_cycle[mCycle + 1])]
            for sl in specified_lines:
                sline = linecache.getline(file[:-4] + '_Parsed SV raw data.txt', sl)
                # for pos, sline in enumerate(fsr):
                #     if (pos + 1) in specified_lines:
                line_dict = eval(sline)
                SV_number = line_dict['Number of SV']
                if SV_number > 0 and line_dict['IOD'] == mIOD:
                    # SV_keys = \
                    #     [(line_dict[f'GNSS and signal type_{n}'], line_dict[f'PRN_{n}']) for n in
                    #      range(SV_number)]
                    # if all(item in RM_keys for item in SV_keys):
                    s_keys = list(line_dict.keys())
                    # get the keys that both measurement raw data and SV data contain
                    # intersection_keys = [item for item in RM_keys if item in SV_keys]
                    # intersection_keys = set(RM_keys) & set(SV_keys)
                    for n in range(SV_number):
                        # 6+n*17--> Elevation angle
                        # 4+n*17 -->PRN, 3+n*17 --> GNSS and signal type
                        SV_key = (line_dict[s_keys[3 + n * 17]], line_dict[s_keys[4 + n * 17]])
                        if SV_key in RM_keys:
                            if (line_dict[s_keys[6 + n * 17]] >= 5) and \
                                    (0 < line_dict[s_keys[17 + n * 17]] < 1e2) and \
                                    (0 < line_dict[s_keys[18 + n * 17]] < 1e2):
                                pos = [line_dict[s_keys[m + n * 17]] for m in [9, 10, 11]]
                                # Corrected pseudorange
                                # 15+n*17 ---> clock bias, 17+n*17 --> trop delay, 18+n*17 --> Iono delay
                                pseudorange_measurement[SV_key][0] = \
                                    pseudorange_measurement[SV_key][0] + line_dict[s_keys[15 + n * 17]] - \
                                    line_dict[s_keys[17 + n * 17]] - line_dict[s_keys[18 + n * 17]]
                                pseudorange_measurement[SV_key].extend(pos)
                            else:
                                pseudorange_measurement[SV_key][0] = 0
                            # break
                    # s_line_number = s_line_number + 1
            for key in list(pseudorange_measurement.keys()):
                if len(pseudorange_measurement[key]) != 5 or pseudorange_measurement[key][0] == 0:
                    del pseudorange_measurement[key]
            # sat_pos = np.array[pseudorange_measurement[]]
            sat_pos = np.array(
                [[pseudorange_measurement[key][2], pseudorange_measurement[key][3], pseudorange_measurement[key][4]]
                 for key in pseudorange_measurement.keys()])
            pseudorange = np.array([pseudorange_measurement[key][0] for key in pseudorange_measurement.keys()])
            pseudorange_weight = \
                1 / np.array([pseudorange_measurement[key][1] for key in pseudorange_measurement.keys()])
            # print(sat_pos)
            # print(pseudorange)
            print('##############################################################################')
            if len(sat_pos) > 3:
                # if len(TOW_pair) == 2:
                print(f'On TOW: {TOW}')

                print(f'Initially, {len(sat_pos)} observations are used for positioning.')
                x, res, iterations, sat_number = least_squares(sat_pos, pseudorange, pseudorange_weight)
                # print(np.linalg.norm(x))
                # x, res = pr_residual(sat_pos, pseudorange, weights=1)
                lat, lon, alt = wgs.ecef2geodetic(x[0], x[1], x[2])
                if if_plot:
                    if lat0 == -1 and lon0 == -1:
                        lat0, lon0, alt0 = lat, lon, alt
                    e, n, __ = wgs.ecef2enu(x[0], x[1], x[2], lat0, lon0, alt0)
                    ax.cla()
                    xs.append(e)
                    ys.append(n)
                    ax.plot(xs, ys, color='r', linewidth=3)
                    plt.xlabel('East/m')
                    plt.ylabel('North/m')
                    plt.title('Trajectory')
                    # ax.xlabel()
                    # ax.legend()
                    plt.pause(0.2)
                # ax.plot(x, y, 'r-', linewidth=2)
                if res < 5e2:
                    coord.extend([str(lon), str(lat), '0'])
                print('The position is', end=':')
                print([lat, lon, alt])
                print('The residual is', end=':')
                print(res)
                print(f'Iterations: {iterations}')
                print(f'Finally, {sat_number} observations are used for positioning')
            else:
                print('Observations are not enough!')
            pseudorange_measurement.clear()
    plt.show(block=True)
    coord = ','.join(coord).replace(',0,', ',0 ')
    create_kml(file[:-4] + '_COROS_weighted_LAD_l1_only', coord)
