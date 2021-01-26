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

#file = '\\sti1arch.epfl.ch\lben-archives\\2019 - CURRENT\Mukesh_Archives\Axopatch200B\\09052019/gc_2_1MKCl_noise_0mV.dat'
filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file
#colors=['k','r','b']
root.update()
for i,file in enumerate(filenames):
    dat = uf.OpenFile(file)
    filename = str(os.path.split(file)[1][:-4])
    os.chdir(os.path.dirname(file))
    directory = (str(os.path.split(file)[0]) + os.sep + 'PSD' + '_SavedImages')
    if not os.path.exists(directory):
        os.makedirs(directory)
    f, Pxx = sig.welch(dat['i1'], dat['samplerate'], nperseg=2**18)
    ax.plot(f, Pxx*1e24, label = filename)#, color = colors[i])

ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel(r'PSD ($\frac{pA^2}{Hz}$)')
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend()

fig.savefig(directory + os.sep + filename + '_PSD.pdf', transparent=True)
fig.savefig(directory + os.sep + filename + '_PSD.png', dpi=300)

#plt.show()


