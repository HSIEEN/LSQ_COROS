# Data : 2023/3/23 18:14
# Author: Shawn Shi
# Right Reserved By COROS
from raw_parser import SV_raw_parse
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
print('-----------Please choose a SV raw data file (*SV*.txt)------------')
file = filedialog.askopenfile()
if not file:
    sys.exit("******未选择任何文件，自动退出程序！！！******")
file = file.name

flow_control = False
check_checksum = True
if __name__ == '__main__':

    with open(file[:-4] + '_Parsed SV raw data.txt', 'r') as fs:
        for line in fs:
            line_dict = eval(line)
            for key in line_dict:
                if 'SV' in key:
                    pass
                    # print(line_dict[key], end=' ')
                if 'GNSS and signal type' in key:
                    # if line_dict[key] % 2 == 0:
                    print(line_dict[key], end=' ')
                #     print('L1', end=' ')
                # else:
                #     print('L5', end=' ')
                if 'PRN' in key:
                    pass
                    # print(line_dict[key], end=' ')
            print(' ')
