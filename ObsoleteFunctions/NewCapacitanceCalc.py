from ObsoleteFunctions import UsefulFunctions as f
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib
import os
from scipy.optimize import curve_fit
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib.ticker import EngFormatter

Res = EngFormatter(unit='Ω', places=2)
Cap = EngFormatter(unit='F', places=2)
Amp = EngFormatter(unit='A', places=2)

plots = 1
lowpass = 5000
AppliedVoltageLoweringFactor = 10
filenames = ['/Users/migraf/Desktop/FG_100mV_30Hz.dat']

def current(V, R, C, I0):
    return V/R + C * np.gradient(V) * inputData['samplerate'] + I0

#filenames = askopenfilenames()  # show an "Open" dialog box and return the path to the selected file

for file in filenames:
    fi = open(str(os.path.split(file)[0]) + os.sep + 'Capacitances.txt', 'a')
    inputData = f.ImportAxopatchData(file)
    Vapplied = inputData['i2']/AppliedVoltageLoweringFactor
    Wn = round(2 * lowpass / (inputData['samplerate']), 4)  # [0,1] nyquist frequency
    b, a = signal.bessel(4, Wn, btype = 'low', analog = False)
    traceApplied = signal.filtfilt(b, a, Vapplied)
    traceOut = signal.filtfilt(b, a, inputData['i1'])
    derivative = np.gradient(traceApplied)
    time = np.float32(np.arange(0, len(inputData['i2'])) / inputData['samplerate'])
    popt, pcov = curve_fit(current, traceApplied, traceOut)#, p0=[10e6, 10e-12, 300e-12])#, method='lm', xtol=1e-20, ftol=1e-12, gtol=1e-20, p0=[10e6, 10e-12])
    fig = plt.figure(1, figsize=(5, 4))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.plot(time, traceApplied, 'b')
    ax1.plot(time, derivative, 'r')
    ax2.plot(time, traceOut, 'y')
    ax2.plot(time, current(traceApplied, *popt), 'g')
    ax1.xaxis.set_major_formatter(EngFormatter(unit='s'))
    ax1.yaxis.set_major_formatter(EngFormatter(unit='A'))
    ax2.yaxis.set_major_formatter(EngFormatter(unit='A'))
    ax1.set_ylabel('Current')
    ax2.set_ylabel('Current')
    ax1.set_xlabel('Time')
    print('C = {}, R = {}, I0 = {}'.format(Cap.format_data(popt[1]), Res.format_data(popt[0]), Amp.format_data(popt[2])))
    plt.show()