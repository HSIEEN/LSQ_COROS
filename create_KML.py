# Data : 2023/3/15 15:15
# Author: Shawn Shi
# Right Reserved By COROS

def create_kml(time_name, cor_str):
    if len(cor_str):
        read_data(cor_str, "gps_ori_data.kml", "%s.kml" % time_name)


# calculate coordinates from the start point and all velocity information

def read_data(gps_data, file_kml, paths):
    """传入经纬度信息写入kml文件，区分原始数据与算法处理过的数据"""
    with open(file_kml, 'r', encoding='utf-8') as fp:
        lines = fp.readlines()
    for i in range(len(lines) - 1):
        if '<coordinates>' in lines[i]:
            lines[i] = '<coordinates>' + gps_data + '</coordinates>\n'
    with open(paths, 'w', encoding='utf-8') as fi:
        fi.writelines(lines)
