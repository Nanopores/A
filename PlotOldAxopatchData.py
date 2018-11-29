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

filename = '/Volumes/lben-archives/2017/Michael/Axopatch/20170227/Roche1_somthing.dat'
inp = uf.OpenFile(filename, ChimeraLowPass=50e3)

fig = plt.figure(1, figsize=(10, 4))
ax = fig.add_subplot(111)
ax.plot(inp['i1'])
plt.show()