import numpy as np
import scipy
import scipy.signal as sig
import Functions as uf
import pyqtgraph as pg
import os
import matplotlib.pyplot as plt
import matplotlib
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
from matplotlib.font_manager import FontProperties
import platform
import csv

fontP = FontProperties()
fontP.set_size('small')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib.ticker import EngFormatter

file =  '\\\\sti1arch.epfl.ch\\lben-archives\\2019 - CURRENT\\Mukesh_Archives\\Axopatch200B\\16052019\\GlassChip02_1MKCl_beforeECR_noise_0mV.dat'
dat = uf.OpenFile(file)
fs = dat['samplerate']
endp = np.uint(147.658 * dat['samplerate'])
i1 = dat['i1'][0:endp]
i1Fil = scipy.signal.savgol_filter(i1, 301, 5)
laserOn = [[7.4, 57], [90, 128.4]]
laserOff = [[59, 88.8], [129.9, 147.6]]
t = np.arange(len(i1))/dat['samplerate']
l = laserOn[1]
l_off = laserOff[1]
s_on = np.uint(l[0] * dat['samplerate'])
e_on = np.uint(l[1] * dat['samplerate'])
s_off = np.uint(l_off[0] * dat['samplerate'])
e_off = np.uint(l_off[1] * dat['samplerate'])

fig = plt.figure(1, figsize=(10, 8))
ax = fig.add_subplot(111)

f_on, Pxx_den_on = sig.welch(i1[s_on:e_on], fs, nperseg=20*256)
f_off, Pxx_den_off = sig.welch(i1[s_off:e_off], fs, nperseg=20*256)

ax.plot(f_on, Pxx_den_on*1e18, 'b', label = 'Laser On')
ax.plot(f_off, Pxx_den_off*1e18, 'r', label = 'Laser Off')
#plt.ylim([0.5e-3, 1])
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel(r'PSD ($\frac{nA^2}{Hz}$)')
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend()
fig.savefig('/Users/migraf/Desktop/PSD_Laser2.pdf', transparent=True)
plt.show()