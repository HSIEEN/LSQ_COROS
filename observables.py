# Data : 2023/4/3 10:03
# Author: Shawn Shi
# Right Reserved By COROS
# Generate PR observables file
import os.path
import sys

import numpy as np
from numpy import NaN, nan
import linecache
from raw_parser import Measurement_raw_parse
from raw_parser import SV_raw_parse

flow_control = False
check_checksum = False

# for placeholder
x = nan


def PR_observable_parse(file):
    # Parse binary files to generate parsed SV and Measurement raw data files
    with open(file, 'rb') as fo:
        data = fo.read()
        # n = n+1
        parsed_measurement_file = file[:-4] + '_Parsed Measurement raw data.txt'
        parsed_sv_file = file[:-4] + '_Parsed SV raw data.txt'
        if not os.path.exists(parsed_measurement_file):
            with open(parsed_measurement_file, 'w') as fm:
                for it in Measurement_raw_parse(data, flow_control, check_checksum):
                    if it['Clock status'] >= 2:
                        fm.writelines(str(it) + "\n")
        if not os.path.exists(parsed_sv_file):
            with open(parsed_sv_file, 'w') as fs:
                for it in SV_raw_parse(data, flow_control, check_checksum):
                    fs.writelines(str(it) + "\n")
    # define pseudorange dictionary
    pseudorange_measurement = {}
    # coord = []
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
            # mdict_cycle[mCycle] = cIOD
            line_number += 1
             # line_dict.clear()
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
            # line_dict.clear()
        sdict_cycle[sCycle + 1] = line_number - 1
    if len(mdict_cycle) != len(sdict_cycle):
        sys.exit("******原始观测量与卫星信息不匹配！！！！******")
    m_line_number = 1
    # s_line_number = 1
    m = len(mdict_cycle)
    # TOW_pair = []
    # lat0, lon0 = -1, -1
    #  define observables file and store pseudorange dictionary data
    observable_file = file[:-4] + '_PR_observables.txt'
    if not os.path.exists(observable_file):
        with open(observable_file, 'w') as fob:
            with open(file[:-4] + '_Parsed Measurement raw data.txt', 'r') as fmr:
                for mline in fmr:
                    # get current number of cycle
                    for i in range(m - 1):
                        if mdict_cycle[i] <= m_line_number < mdict_cycle[i + 1]:
                            mCycle = i
                            break
                    # line_dict.clear()
                    line_dict = eval(mline)
                    m_line_number = m_line_number + 1
                    # print(line_dict['IOD'])
                    mIOD = line_dict['IOD']
                    # for debugging
                    # if line_dict['Receiver TOW'] < 527310.0 or line_dict['Receiver TOW'] > 527330.0:
                    #     continue
                    if line_dict['Number of observation'] > 0:
                        m_keys = list(line_dict.keys())
                        for n in range(line_dict['Number of observation']):
                            # 8+n*12 --> PRN code, 7+n*12 --> GNSS and signal type
                            RM_key = (line_dict[m_keys[7 + n * 12]], line_dict[m_keys[8 + n * 12]])
                            # if line_dict['Receiver TOW'] != 376090.001:
                            #     continue
                            # 12+i*12 --> raw pseudorange, 14+i*14 --> pseudorange variance index
                            if (line_dict[m_keys[14 + n * 12]] <= 32) and\
                                    (1.5e7 < line_dict[m_keys[12 + n * 12]] < 4.5e7):
                                pseudorange_measurement[RM_key] =\
                                    [line_dict['Receiver TOW'], line_dict[m_keys[12 + n * 12]],
                                     line_dict[m_keys[14 + n * 12]]]
                    else:
                        continue
                    RM_keys = [key for key in pseudorange_measurement]
                    if len(RM_keys) > 0:
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
                                s_keys = list(line_dict.keys())
                                for n in range(SV_number):
                                    # 6+n*17--> Elevation angle
                                    # 4+n*17 -->PRN, 3+n*17 --> GNSS and signal type
                                    # SV_key, L1 keys
                                    SV_key = (line_dict[s_keys[3 + n * 17]], line_dict[s_keys[4 + n * 17]])
                                    if SV_key in RM_keys:
                                        pos = [line_dict[s_keys[m + n * 17]] for m in [9, 10, 11]]
                                        if (line_dict[s_keys[6 + n * 17]] >= 5) and \
                                                (abs(line_dict[s_keys[15 + n * 17]]) < 3e5) and \
                                                (0 < line_dict[s_keys[17 + n * 17]] < 1e2) and \
                                                (0 < line_dict[s_keys[18 + n * 17]] < 1e2):

                                            # Corrected pseudorange
                                            # 15+n*17 ---> clock bias, 17+n*17 --> trop delay, 18+n*17 --> Iono delay
                                            pseudorange_measurement[SV_key][1] = \
                                                pseudorange_measurement[SV_key][1] + line_dict[s_keys[15 + n * 17]] - \
                                                line_dict[s_keys[17 + n * 17]] - line_dict[s_keys[18 + n * 17]]
                                        if 1e6 < np.max(np.absolute(np.array(pos))) < 4.5e7:
                                            pseudorange_measurement[SV_key].extend(pos)
                                        else:
                                            pseudorange_measurement[SV_key][1] = 0
                                    # L5 keys
                                    SV_key5 = (SV_key[0] + 1, SV_key[1])
                                    if SV_key5 in RM_keys:
                                        pos = [line_dict[s_keys[m + n * 17]] for m in [9, 10, 11]]
                                        # 6+n*17--> Elevation angle
                                        if (line_dict[s_keys[6 + n * 17]] >= 5) and \
                                                (abs(line_dict[s_keys[15 + n * 17]]) < 3e5) and \
                                                (0 < line_dict[s_keys[17 + n * 17]] < 1e2) and \
                                                (0 < line_dict[s_keys[18 + n * 17]] < 1e2):
                                            # pos = [line_dict[s_keys[m + n * 17]] for m in [9, 10, 11]]
                                            pseudorange_measurement[SV_key5][1] = \
                                                pseudorange_measurement[SV_key5][1] + line_dict[s_keys[15 + n * 17]] - \
                                                line_dict[s_keys[17 + n * 17]] - line_dict[s_keys[18 + n * 17]]
                                        if 1e6 < np.max(np.absolute(np.array(pos))) < 4.5e7:
                                            pseudorange_measurement[SV_key5].extend(pos)
                                        else:
                                            pseudorange_measurement[SV_key5][1] = 0
                                fob.writelines(str(pseudorange_measurement) + '\n')
                                print(f'PR observables file line writing: {str(pseudorange_measurement)}')
                                # clear the dictionary, very important!!!
                                pseudorange_measurement.clear()
        print('PR observables file has been generated!')


if __name__ == '__main__':
    raw_data_file = \
        r'F:\LogData\RawData\B20\国际创新园\Binary_file\log_system_B20_6#_V2.5.0_PATCH2_bd3b50_2023-03-09_14_21_49.txt'
    PR_observable_parse(raw_data_file)
