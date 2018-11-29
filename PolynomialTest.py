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


def Spec(f, a0, a1, a2, a3):
    return a0 / f + a1 + a2 * f + a3 * f ** 2

fig = plt.figure(1, figsize=(10, 8))
ax = fig.add_subplot(111)

f=np.linspace(1, 1e5, 10000)
ax.plot(f, Spec(f, 1e4, 10, 3.4e-6, 6e-4))
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel(r'PSD ($\frac{pA^2}{Hz}$)')
ax.set_yscale('log')
ax.set_xscale('log')

plt.show()