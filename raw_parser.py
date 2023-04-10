# Data : 2023/3/15 15:15
# Author: Shawn Shi
# Right Reserved By COROS

import struct
import sys
import numpy as np
import math


def calc_checksum(checksum_range: bytes) -> int:
    checksum_calc = checksum_range[0]
    for it in checksum_range[1:]:
        checksum_calc ^= it
    return checksum_calc


def convert_flow_control(data: bytes) -> bytes:
    converted_bytes = bytearray()
    need_convert = False
    for it in data:
        if need_convert:
            converted_bytes.append(0xff - it)
            need_convert = False
            continue
        if it == 0x77:
            need_convert = True
            continue
        converted_bytes.append(it)
    return converted_bytes


def PVT_raw_parse(data: bytes, flow_control=True, check_checksum=True):
    """
        param data: binary data from raw measurement log
        :param flow_control:  if software flow is used
        :param check_checksum: if check
        :return: a dictionary with PVT data
    """
    PVT_preamble = b'\x04\x24'
    # SV_preamble = b''
    # SV_end = b''
    PVT_end = b'\xAA\x44'
    start_pos = 0
    checksum_error_cnt = 0
    length_error = 0
    while True:
        start_pos = data.find(PVT_preamble, start_pos)
        if start_pos == -1:
            break
        end_pos = data.find(PVT_end, start_pos)
        if end_pos == -1:
            break

        message_bytes = data[start_pos:end_pos + 2]
        start_pos = end_pos + 2

        if not message_bytes:
            continue

        if flow_control:
            message_bytes = convert_flow_control(message_bytes)

        if check_checksum:
            checksum_calc = calc_checksum(message_bytes[2:-3])  # message_bytes[2:-3]
            if checksum_calc != message_bytes[-3]:
                checksum_error_cnt += 1
                continue

        message_id = struct.unpack('<H', message_bytes[2:4])[0]
        if message_id != 4005:
            continue
        length = int.from_bytes(message_bytes[4:6], byteorder='little')
        msg_data = message_bytes[6:6 + length]
        if len(msg_data) != 120:
            length_error += 1
            continue

        datadict = {}
        (
            datadict['iod'],  # U1
            datadict['fix_type'],  # U1
            datadict['year'],  # U2
            datadict['month'],  # U1
            datadict['day'],  # U1
            datadict['hour'],  # U1
            datadict['minute'],  # U1
            datadict['ms'],  # U2
            datadict['second'],  # U1
            pvt_status,  # U1,
            datadict['lat'],  # I4
            datadict['lon'],  # I4
            datadict['alt'],  # I4
            datadict['geoid'],  # I4
            datadict['speed'],  # U4
            datadict['heading'],  # U4
            datadict['horizontal_accuracy'],  # R4
            datadict['lsq_lat'],  # I4
            datadict['lsq_lon'],  # I4
            datadict['lsq_alt'],  # I4
            datadict['lsq_horizontal_velocity_N'],  # I4
            datadict['lsq_horizontal_velocity_E'],  # I4
            datadict['lsq_horizontal_velocity_D'],  # I4
            datadict['lsq_clock_bias'],  # R8
            datadict['lsq_clock_drift'],  # R8
            datadict['kf_lat'],  # I4
            datadict['kf_lon'],  # I4
            datadict['kf_alt'],  # I4
            datadict['kf_horizontal_velocity_N'],  # I4
            datadict['kf_horizontal_velocity_E'],  # I4
            datadict['kf_horizontal_velocity_D'],  # I4
            datadict['kf_clock_bias'],  # R8
            datadict['kf_clock_drift'],  # R8
        ) = struct.unpack('<BBHBBBBHBBiiiiIIfiiiiiiddiiiiiidd', msg_data)
        # U1->B, U2->H, I4->i, U4->I, R4->f, R8-> d
        datadict['lla_valid'] = bool(pvt_status & 0b1)
        datadict['lsq_valid'] = bool(pvt_status >> 1 & 0b1)
        datadict['kf_valid'] = bool(pvt_status >> 2 & 0b1)

        yield datadict

    print("Checksum error", checksum_error_cnt, file=sys.stderr)
    print("PVT length error", length_error, file=sys.stderr)


def SV_raw_parse(data: bytes, flow_control=True, check_checksum=True):
    S_preamble = b'\x04\x24'
    # SV_preamble = b''
    # SV_end = b''
    S_end = b'\xAA\x44'
    start_pos = 0
    checksum_error_cnt = 0
    length_error = 0
    while True:
        start_pos = data.find(S_preamble, start_pos)
        if start_pos == -1:
            break
        end_pos = data.find(S_end, start_pos)
        if end_pos == -1:
            break

        message_bytes = data[start_pos:end_pos + 2]
        restart_pos = start_pos
        start_pos = end_pos + 2

        if not message_bytes:
            continue

        if flow_control:
            message_bytes = convert_flow_control(message_bytes)

        if check_checksum:
            checksum_calc = calc_checksum(message_bytes[2:-3])  # message_bytes[2:-3]
            if checksum_calc != message_bytes[-3]:
                checksum_error_cnt += 1
                continue

        message_id = struct.unpack('<H', message_bytes[2:4])[0]
        if message_id != 4003:
            continue
        length = int.from_bytes(message_bytes[4:6], byteorder='little')
        # number_of_obs = int.from_bytes(message_bytes[19], byteorder='little')
        number_of_SVs = message_bytes[7]
        # reassign bytes to message_bytes
        message_bytes = data[restart_pos: restart_pos + length]
        msg_data = message_bytes[6:6 + length]
        # when length error occurs, get maximum number of SVs
        if len(msg_data) != (3 + 52 * number_of_SVs):
            length_error += 1
            number_of_SVs = math.floor((len(msg_data) - 3) / 52)
            if number_of_SVs < 1:
                continue
        datadict = {}
        (  # Part I
            datadict['IOD'],  # U1
            datadict['Number of SV'],  # U1
            datadict['Reserved1'],  # U1

        ) = struct.unpack('<BBB', msg_data[0:3])  # partII BBBBfddBBHHH
        # U1->B, U2->H, I4->i, U4->I, R4->f, R8-> d, I1->b
        datadict['Number of SV'] = number_of_SVs
        if datadict['Number of SV'] != 0:
            for n in range(0, datadict['Number of SV']):
                (
                    datadict['GNSS and signal type_' + str(n)],  # U1
                    datadict['PRN_' + str(n)],  # U1
                    # pvt_status,  # U1,
                    datadict['Azimuth_' + str(n)],  # U2
                    datadict['Elevation_' + str(n)],  # I1
                    datadict['SV status_' + str(n)],  # U1
                    datadict['Reserved_' + str(n)],  # U2
                    datadict['Position X_' + str(n)],  # R4
                    datadict['Position Y_' + str(n)],  # R4
                    datadict['Position Z_' + str(n)],  # R4
                    datadict['Velocity X_' + str(n)],  # R4
                    datadict['Velocity Y_' + str(n)],  # R4
                    datadict['Velocity Z_' + str(n)],  # R4
                    datadict['Clock bias_' + str(n)],  # R4
                    datadict['Clock drift_' + str(n)],  # R4
                    datadict['Trop delay_' + str(n)],  # R4
                    datadict['Iono delay_' + str(n)],  # R4
                    datadict['Position variance_' + str(n)],  # R4
                ) = struct.unpack('<BBHbBHfffffffffff', msg_data[3 + n * 52:3 + (n + 1) * 52])
        yield datadict

    print("Checksum error", checksum_error_cnt, file=sys.stderr)
    print("SV length error", length_error, file=sys.stderr)


def Measurement_raw_parse(data: bytes, flow_control=True, check_checksum=True):
    """
    :param data: binary data from MTK raw log
    :param flow_control: if software decoding is used
    :param check_checksum: if CRC checksum needed
    :return: a dictionary containing decoded text information
    """
    m_preamble = b'\x04\x24'
    # SV_preamble = b''
    # SV_end = b''
    m_end = b'\xAA\x44'
    start_pos = 0
    checksum_error_cnt = 0
    length_error = 0
    while True:
        start_pos = data.find(m_preamble, start_pos)
        if start_pos == -1:
            break
        end_pos = data.find(m_end, start_pos)
        if end_pos == -1:
            break

        message_bytes = data[start_pos:end_pos + 2]
        # restart_pos = start_pos
        start_pos = end_pos + 2

        if not message_bytes:
            continue

        if flow_control:
            message_bytes = convert_flow_control(message_bytes)

        if check_checksum:
            checksum_calc = calc_checksum(message_bytes[2:-3])  # message_bytes[2:-3]
            if checksum_calc != message_bytes[-3]:
                checksum_error_cnt += 1
                continue

        message_id = struct.unpack('<H', message_bytes[2:4])[0]
        if message_id != 4001:
            continue
        length = int.from_bytes(message_bytes[4:6], byteorder='little')
        # number_of_obs = int.from_bytes(message_bytes[4:6], byteorder='little')
        number_of_obs = int(message_bytes[18])
        # number_of_obs = struct.unpack('<B', message_bytes[18])
        # message_bytes = data[restart_pos: restart_pos + length]
        msg_data = message_bytes[6:6 + length]
        # when length error occurs, get maximum number of observations
        if len(msg_data) != (16 + 32 * number_of_obs):
            length_error += 1
            # number_of_obs = math.floor((len(msg_data)-16)/32)
            # if number_of_obs < 1:
            continue
        datadict = {}
        (  # Part I
            datadict['Receiver TOW'],  # R8
            datadict['Receiver week number'],  # U2
            datadict['Clock status'],  # U1
            datadict['Leap second'],  # I1
            datadict['Number of observation'],  # U1
            datadict['IOD'],  # U1
            datadict['Reserved1'],  # U2
        ) = struct.unpack('<dHBbBBH', msg_data[0:16])  # partII BBBBfddBBHHH
        # U1->B, U2->H, I4->i, U4->I, R4->f, R8-> d, I1->b
        datadict['Number of observation'] = number_of_obs
        if datadict['Number of observation'] != 0:
            for n in range(0, datadict['Number of observation']):
                (
                    datadict['GNSS and signal type_' + str(n)],  # U1
                    datadict['PRN_' + str(n)],  # U1
                    # pvt_status,  # U1,
                    datadict['GLONASS frequency_' + str(n)],  # U1
                    datadict['C/N0_' + str(n)],  # U1
                    datadict['Doppler_' + str(n)],  # R4
                    datadict['Pseudorange_' + str(n)],  # R8
                    datadict['Carrier phase_' + str(n)],  # R8
                    datadict['Pseudorange variance index_' + str(n)],  # U1
                    datadict['Doppler variance index_' + str(n)],  # U1
                    datadict['Lock time_' + str(n)],  # U2
                    datadict['Channel status_' + str(n)],  # U2
                    datadict['Reserved2_' + str(n)],  # U2
                ) = struct.unpack('<BBBBfddBBHHH', msg_data[16 + n * 32:16 + (n + 1) * 32])
        yield datadict

    print("Checksum error", checksum_error_cnt, file=sys.stderr)
    print("Measurement length error", length_error, file=sys.stderr)


def PR_raw_parse(data: bytes, flow_control=True, check_checksum=True):
    """
    :param data:
    :param flow_control:
    :param check_checksum:
    :return:
    """
    m_preamble = b'\x04\x24'
    # SV_preamble = b''
    # SV_end = b''
    m_end = b'\xAA\x44'
    start_pos = 0
    checksum_error_cnt = 0
    while True:
        start_pos = data.find(m_preamble, start_pos)
        if start_pos == -1:
            break
        end_pos = data.find(m_end, start_pos)
        if end_pos == -1:
            break

        message_bytes = data[start_pos:end_pos + 2]
        start_pos = end_pos + 2

        if not message_bytes:
            continue

        if flow_control:
            message_bytes = convert_flow_control(message_bytes)

        if check_checksum:
            checksum_calc = calc_checksum(message_bytes[2:-3])  # message_bytes[2:-3]
            if checksum_calc != message_bytes[-3]:
                checksum_error_cnt += 1
                continue

        message_id = struct.unpack('<H', message_bytes[2:4])[0]
        if message_id != 4006:
            continue
        length = int.from_bytes(message_bytes[4:6], byteorder='little')
        # number_of_obs = int.from_bytes(message_bytes[19], byteorder='little')
        number_of_PR = message_bytes[10]
        msg_data = message_bytes[6:6 + length]
        if len(msg_data) != (8 + 8 * number_of_PR):
            continue
        datadict = {}
        (  # Part I
            datadict['Receiver TOW'],  # U4
            datadict['Number of PR residuals'],  # U1
            datadict['Reserved1'],  # U3

        ) = struct.unpack('<fB3s', msg_data[0:8])  # partII BBBBfddBBHHH
        # U1->B, U2->H, I4->i, U4->I, R4->f, R8-> d, I1->b
        if datadict['Number of PR residuals'] != 0:
            for n in range(0, datadict['Number of PR residuals']):
                (
                    datadict['GNSS type_' + str(n)],  # U1
                    datadict['PRN_' + str(n)],  # U1
                    # pvt_status,  # U1,
                    datadict['Signal type_' + str(n)],  # U1
                    datadict['Reserved2_' + str(n)],  # U1
                    datadict['PR residual_' + str(n)],  # I4
                ) = struct.unpack('<BBBBi', msg_data[8 + n * 8:8 + (n + 1) * 8])
        yield datadict

    print("Checksum error", checksum_error_cnt, file=sys.stderr)


def velocity_track_parse(file, flow_control, check_checksum):
    lsq_coord = []
    kf_coord = []
    # calculate velocity, vertical velocity excluded
    # lsq_valid_velocity_vector = np.zeros((1, 2))
    lsq_valid_velocity_vector = [[0, 0]]
    # lsq_valid_velocity_vector[0, :] = 0
    # kf_valid_velocity_vector = np.zeros((1, 2))
    kf_valid_velocity_vector = [[0, 0]]
    # kf_valid_velocity_vector[0, :] = 0
    lsq_horizontal_velocity = []
    kf_horizontal_velocity = []
    # Clock bias and drift
    lsq_clock_bias = []
    lsq_clock_drift = []
    kf_clock_bias = []
    kf_clock_drift = []
    # calculate coordinates from velocity
    Doppler_track = np.zeros((2, 1, 2))
    # Doppler_track[:, 0, :] = 0
    # kf_or_lsq, points, E_or_N
    new_coord = np.zeros((2, 1, 2))
    n = 0
    # fp = open(file[:-4] + '_Parsed PVT raw data.txt', 'w')
    # fm = open(file[:-4] + "_Parsed Raw measurement.txt", 'w')
    # fs = open(file[:-4]+'_Parsed SV raw data.txt', 'w')
    with open(file, 'rb') as f:
        data = f.read()
        # n = n+1
        for it in PVT_raw_parse(data, flow_control, check_checksum):
            # update the track data from velocity data
            if math.sqrt(it['lsq_horizontal_velocity_E'] ** 2 + it['lsq_horizontal_velocity_N'] ** 2) / 1000 > 0:
                new_coord[0, 0, 0] = it['lsq_horizontal_velocity_E'] / 1000 + Doppler_track[0, n, 0]
                new_coord[0, 0, 1] = it['lsq_horizontal_velocity_N'] / 1000 + Doppler_track[0, n, 1]
                # lsq_valid_velocity_vector = np.append(lsq_valid_velocity_vector, [[it['lsq_horizontal_velocity_E']
                # / 1000, it['lsq_horizontal_velocity_N'] / 1000]], axis=0)
            if math.sqrt(it['kf_horizontal_velocity_E'] ** 2 + it['kf_horizontal_velocity_N'] ** 2) / 100 > 0:
                new_coord[1, 0, 0] = it['kf_horizontal_velocity_E'] / 100 + Doppler_track[1, n, 0]
                new_coord[1, 0, 1] = it['kf_horizontal_velocity_N'] / 100 + Doppler_track[1, n, 1]
                # kf_valid_velocity_vector = np.append(kf_valid_velocity_vector, [[it['kf_horizontal_velocity_E'] / 100,
                # it['kf_horizontal_velocity_N'] / 100]], axis=0)

            Doppler_track = np.append(Doppler_track, new_coord, axis=1)
            n = n + 1
            if it['lsq_valid']:
                # m = m+1
                lsq_coord.extend([str(it['lsq_lon'] / 1e7), str(it['lsq_lat'] / 1e7), '0'])
                lsq_horizontal_velocity.append(math.sqrt(it['lsq_horizontal_velocity_E'] ** 2
                                                         + it['lsq_horizontal_velocity_N'] ** 2) / 1000)
                lsq_clock_bias.append(it['lsq_clock_bias'])
                lsq_clock_drift.append(it['lsq_clock_drift'])
                lsq_valid_velocity_vector.append([it['lsq_horizontal_velocity_E'] / 1000,
                                                  it['lsq_horizontal_velocity_N'] / 1000])
            else:
                lsq_horizontal_velocity.append(-0.5)
            if it['kf_valid']:
                kf_coord.extend([str(it['kf_lon'] / 1e7), str(it['kf_lat'] / 1e7), '0'])
                kf_horizontal_velocity.append(math.sqrt(it['kf_horizontal_velocity_E'] ** 2
                                                        + it['kf_horizontal_velocity_N'] ** 2) / 100)
                kf_clock_bias.append(it['kf_clock_bias'])
                kf_clock_drift.append(it['kf_clock_drift'])
                kf_valid_velocity_vector.append([it['kf_horizontal_velocity_E'] / 1000,
                                                 it['kf_horizontal_velocity_N'] / 100])
            else:
                kf_horizontal_velocity.append(-0.5)
        #     fp.writelines(str(it) + "\n")
        # for it in Measurement_raw_parse(data, flow_control, check_checksum):
        #     txt = str(it)
        #     fm.writelines(txt + '\n')
        # for it in SV_raw_parse(data, flow_control, check_checksum):
        #     txt = str(it)
        #     fs.writelines(txt+'\n')
        # print(it)
        # print(m)
        # print(n)
        # fp.close()
        # fm.close()
        # fs.close()
        lsq_coord_data = ','.join(lsq_coord).replace(',0,', ',0 ')
        kf_coord_data = ','.join(kf_coord).replace(',0,', ',0 ')
    lsq_valid_velocity_vector = np.array(lsq_valid_velocity_vector)
    kf_valid_velocity_vector = np.array(kf_valid_velocity_vector)
    return lsq_coord_data, kf_coord_data, lsq_valid_velocity_vector, \
           kf_valid_velocity_vector, lsq_horizontal_velocity, kf_horizontal_velocity, Doppler_track
