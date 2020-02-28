
# -*- coding: utf-8 -*-
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")

import numpy as np
import scipy
from scipy.optimize import curve_fit

import scipy.signal as sig
from ObsoleteFunctions import UsefulFunctions as uf
import pyqtgraph as pg
import os
import matplotlib.pyplot as plt
import matplotlib
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
root = Tk()
root.withdraw()
from matplotlib.font_manager import FontProperties
import platform
import csv
fontP = FontProperties()
fontP.set_size('small')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib.ticker import EngFormatter

fig = plt.figure(1, figsize=(10, 8))
ax = fig.add_subplot(111)

#file = '/Volumes/lben/lben-commun/2018 User Data/Michael/Axopatch/20181003/NIPm10_5nm_1MKCl_pH75_Noise_100kHz_0mV_1.dat'
filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file
#colors=['k','r','b']
root.update()

def func(x, alpha, beta):
    return beta * x**(-alpha) #beta is I2A, x is f


for i,file in enumerate(filenames):
    dat = uf.OpenFile(file)
    filename = str(os.path.split(file)[1][:-4])
    os.chdir(os.path.dirname(file))
    directory = (str(os.path.split(file)[0]) + os.sep + 'PSD' + '_SavedImages')
    if not os.path.exists(directory):
        os.makedirs(directory)
    f, Pxx = sig.welch(dat['i1'], dat['samplerate'], nperseg=2**18)
    ax.plot(f, Pxx*1e24, label = filename)
    print(np.std(dat['i1'])*1e12)

    # #fitting
    # num_data = 200
    # x_data = f[2:num_data]
    # y_data = Pxx * 1e24
    # y_data = y_data[2:num_data]
    # popt, pcov = curve_fit(func, x_data, y_data)
    # y_fit = func(x_data, *popt)
    # ax.plot(x_data, y_fit, label = "alpha = {}, beta = {}".format(popt[0], popt[1]))
    # print(popt)
    # #ax.text(popt)

ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel(r'PSD ($\frac{pA^2}{Hz}$)')
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend()
plt.show()

# fig.savefig(directory + os.sep + filename + '_PSD.pdf', transparent=True)
# fig.savefig(directory + os.sep + filename + '_PSD.png', dpi=300)



#########################################
##### Curve fitting #####################
#########################################
# print(np.shape(f))

# x_data = f[2:50]
# y_data = Pxx*1e24
# y_data = y_data[2:50]
# popt, pcov = curve_fit(func, x_data, y_data)
# y_fit = func(x_data, *popt)
# ax.plot(x_data, y_fit)
# print(popt)
# plt.show()

