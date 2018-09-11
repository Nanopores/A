import UsefulFunctions as f
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from peakutils.peak import indexes
import pyqtgraph as pg
import os
from tkinter import Tk
from tkinter.filedialog import askopenfilenames

plots=1
lowpass=5e3
#filenames = ['/Users/migraf/Desktop/Capacitance Measurements/FgenModelCellSetUp.dat']

filenames = askopenfilenames()  # show an "Open" dialog box and return the path to the selected file

for file in filenames:
    fi = open(str(os.path.split(file)[0]) + os.sep + 'Capacitances.txt', 'a')

    inputData = f.ImportAxopatchData(file)
    trace = inputData['i2']
    Wn = round(2 * lowpass / (inputData['samplerate']), 4)  # [0,1] nyquist frequency
    b, a = signal.bessel(4, Wn, btype='low', analog=False)
    trace = signal.filtfilt(b, a, trace)

    derivative = np.gradient(trace)
    time = np.float32(np.arange(0, len(inputData['i2'])) / inputData['samplerate'])
    dertime=time[:-1]

    #peakind = signal.find_peaks_cwt(derivative, np.arange(1,10))
    #peakind, dertime[peakind], derivative[peakind]

    print('Detect peaks with minimum height and distance filters.')
    indexesPos = indexes(np.array(derivative), thres=0.75, min_dist=1000)
    indexesNeg = indexes(np.array(np.negative(derivative)), thres=0.75, min_dist=1000)

    #Get Local Minimas
    AllData=np.zeros((4,len(indexesPos)),dtype=np.uint32)
    NewData=np.zeros(len(indexesPos))
    for l,i in enumerate(indexesPos):
        k=i
        while derivative[k]>derivative[k-1] and k>0:
            k-=1
        minimaLeft=k
        k=i
        if k<len(derivative)-2:
            while derivative[k]>derivative[k+1] and k<len(derivative)-2:
                k+=1
        minimaRigth=k
        AllData[0,l]=np.uint32(i)
        AllData[1,l]=np.uint32(minimaLeft)
        AllData[2,l]=np.uint32(minimaRigth)
        NewData[l]=np.abs(trace[np.uint32(minimaLeft)]-trace[np.uint32(minimaRigth)])

    fi.write('{}\nThe mean is: {}F, std is: {}F\n\n'.format(str(os.path.split(file)[1]), pg.siFormat(np.mean(NewData)/24), pg.siFormat(np.std(NewData)/24)))
    #fi.write('{}\nThe mean of dI: {}A, std of dI: {}A\n\n'.format(str(os.path.split(file)[1]), pg.siFormat(np.mean(NewData)), pg.siFormat(np.std(NewData))))

    fi.close

    if plots:
        both=np.concatenate((indexesPos,indexesNeg))
        fig, ax1 = plt.subplots()
        ax1.plot(time, derivative, 'b')
        plt.hold
        ax1.plot(time[AllData[0, :]], derivative[AllData[0, :]], 'r.')
        ax1.plot(time[AllData[1, :]], derivative[AllData[1, :]], 'g.')
        ax1.plot(time[AllData[2, :]], derivative[AllData[2, :]], 'y.')
        fig.tight_layout()

        fig, ax1 = plt.subplots()
        ax1.plot(time, trace, 'b')
        plt.hold
        ax1.plot(time[AllData[0, :]], trace[AllData[0, :]], 'r.')
        ax1.plot(time[AllData[1, :]], trace[AllData[1, :]], 'g.')
        ax1.plot(time[AllData[2, :]], trace[AllData[2, :]], 'y.')
        fig.tight_layout()
        plt.show()

