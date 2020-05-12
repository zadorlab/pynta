# from datetime import datetime
import datetime
from statistics import mean, median
import os


def getAvTime(path):
    time_list = []
    for struc in sorted(os.listdir(path), key=str):
        strucPath = os.path.join(path, struc)
        if os.path.isdir(strucPath):
            for timelog in os.listdir(strucPath):
                if timelog.endswith('time.log'):
                    timelogPath = os.path.join(strucPath, timelog)
                    with open(timelogPath, 'r') as f:
                        lines = f.readlines()
                        for theTime in lines:
                            theTime.strip()
                            pass
                        time_list.append(theTime[:-1])
                    f.close()

    totSec = []
    print('conf', '\t', 'time')
    print('------------------')
    for num, i in enumerate(time_list):
        print(str(num).zfill(3), i)
        i = i[:-7]
        i = i.split(':')  # Hours:Minutes:Second
        h = float(i[0])
        m = float(i[1])
        s = float(i[2])
        h = h * 3600
        m = m * 60
        total = h + m + s
        totSec.append(total)
    avTimeSec = mean(totSec)
    mdTimeSec = median(totSec)
    avTime = str(datetime.timedelta(seconds=avTimeSec))
    mdTime = str(datetime.timedelta(seconds=mdTimeSec))
    print('------------------')
    txt = 'time of calculations is '
    print('The average {}{}\nThe  median {}{}'.format(
        txt, avTime, txt, mdTime))
