# Data : 2023/3/23 16:15
# Author: Shawn Shi
# Right Reserved By COROS
import numpy as np
import matplotlib.pyplot as plt
from track_parser import track_parser, wgs_enu2lla
from raw_parser import PVT_raw_parse
from create_KML import create_kml
from tkinter import filedialog
import tkinter as tk
import sys
root = tk.Tk()
root.withdraw()
print('---------LSQ and KF combine filtering v1.0.0-----------')
print('-----------All rights are reserved by COROS------------')
print('-----------Please choose a raw data file (*.txt)------------')
file = filedialog.askopenfile()
if not file:
    sys.exit("******未选择任何文件，自动退出程序！！！******")
file = file.name

dt = 1
# n = 0
plt.figure()
zx = []
zy = []

# ##########################################################################
lsq_coord = []
kf_coord = []
# calculate velocity, vertical velocity excluded
# lsq_velocity_vector = np.zeros((1, 2))
lsq_velocity_vector = [[0, 0]]
# lsq_velocity_vector[0, :] = 0
# kf_velocity_vector = np.zeros((1, 2))
kf_velocity_vector = [[0, 0]]
# kf_velocity_vector[0, :] = 0

# kf_or_lsq, points, E_or_N
# n = 0
flow_control = False
check_checksum = True
fi = open(file[:-4] + '_Parsed raw data.txt', 'a')
# raw data information abstract
with open(file, 'rb') as fo:
    data = fo.read()
    # n = n+1
    for it in PVT_raw_parse(data, flow_control, check_checksum):
        if it['kf_valid']:
            lsq_coord.extend([str(it['lsq_lon'] / 1e7), str(it['lsq_lat'] / 1e7), '0'])
            lsq_velocity_vector.append([it['lsq_horizontal_velocity_E'] / 1000,
                                        it['lsq_horizontal_velocity_N'] / 1000])
            kf_coord.extend([str(it['kf_lon'] / 1e7), str(it['kf_lat'] / 1e7), '0'])

            kf_velocity_vector.append([it['kf_horizontal_velocity_E'] / 100,
                                       it['kf_horizontal_velocity_N'] / 100])
            fi.writelines(str(it) + "\n")

    fi.close()
    lsq_coord_data = ','.join(lsq_coord).replace(',0,', ',0 ')
    kf_coord_data = ','.join(kf_coord).replace(',0,', ',0 ')

lsq_velocity = np.array(lsq_velocity_vector)
lsq_velocity = np.delete(lsq_velocity, 0, axis=0)
kf_velocity = np.array(kf_velocity_vector)
kf_velocity = np.delete(kf_velocity, 0, axis=0)
create_kml(file[:-4] + '_lsq_track_raw', lsq_coord_data)
create_kml(file[:-4] + '_kf_track_raw', kf_coord_data)
lsq_kml = file[:-4] + '_lsq_track_raw.kml'
kf_kml = file[:-4] + '_kf_track_raw.kml'
kf_track, start_lon, start_lat = track_parser(kf_kml)
lsq_track, _, _ = track_parser(lsq_kml, start_lon, start_lat)

# plt.plot(zx, zy, color='b', label='LSQ raw output')
plt.plot(kf_track[:, 0], kf_track[:, 1], color='g', label='Kalman raw output')
# transform enu coordinates system to lla (WGS84)
# flon, flat = wgs_enu2lla(fx, fy, 0, start_lon, start_lat, 0)
# coord = []
# for n in range(0, len(flon)):
#     coord.extend([str(flon[n]), str(flat[n]), '0'])
# coord = ','.join(coord).replace(',0,', ',0 ')
# create_kml(file[:-4] + '_LSQ_kalman_track_without_velocity', coord)
# fcoord = []
# for n in range(0, len(fvlon)):
#     fcoord.extend([str(fvlon[n]), str(fvlat[n]), '0'])
# fcoord = ','.join(fcoord).replace(',0,', ',0 ')
# create_kml(file[:-4] + '_LSQ_kalman_track_with_velocity', fcoord)
plt.legend()

# do_something_with_estimate (f.x)
plt.show(block=True)
