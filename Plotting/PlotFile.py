import numpy as np
import scipy
from scipy.optimize import curve_fit
from scipy import constants as cst
import Functions as uf
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import h5py
import pandas as pd
from matplotlib.ticker import EngFormatter
from matplotlib import rc
import LoadData
#rc('text', usetex=True)
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
root = Tk()
root.withdraw()
formCurrent = EngFormatter(unit='A', places=3)
formA = EngFormatter(unit='A', places=3)
formB = EngFormatter(unit='s', places=3)


filenames = askopenfilenames()
root.update()

for i,file in enumerate(filenames):
    dat = LoadData.OpenFile(file)
    fig1 = plt.figure(1, figsize=(8, 4))
    ax_part1 = fig1.add_subplot(1, 1, 1)
    t = np.arange(len(dat['i1']))/dat['samplerate']
    ax_part1.plot(t, dat['i1']*1e9)
plt.show()


from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
matplotlib.use("TkAgg")

matplotlib.use("TkAgg")
import numpy as np
import scipy
from scipy.optimize import curve_fit
from scipy import constants as cst
import Functions as uf
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import h5py
import pandas as pd
from matplotlib.ticker import EngFormatter
from matplotlib import rc
import LoadData
#rc('text', usetex=True)
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
import seaborn as sns
from matplotlib import gridspec
root = Tk()
root.withdraw()


formCurrent = EngFormatter(unit='A', places=3)
formA = EngFormatter(unit='A', places=3)
formB = EngFormatter(unit='s', places=3)

filenames = askopenfilenames()
root.update()

# #1  Plot 211 and 212 type graph::: Current-time (top), Voltage-time (below)
#
#
# for i,file in enumerate(filenames):
#     dat = LoadData.OpenFile(file)
#     fig1 = plt.figure(1, figsize=(8, 4))
#     ax_part1 = fig1.add_subplot(2, 1, 1)
#     t = np.arange(len(dat['i1']))/dat['samplerate']
#     ax_part1.plot(t, dat['i1'] * 1e9)
#     ax_part1.set_ylabel('Ionic Current (nA)')
#     current_array = np.array(dat['i1'] * 1e9)
#     # print("mean/std at 0 mV : ", np.mean(current_array[0:47821], axis = 0), np.std(current_array[0:47821], axis = 0))
#     # print("mean/std at 50 mV  : ", np.mean(current_array[49992:20274], axis=0),  np.std(current_array[49992:20274], axis = 0))
#     # print("mean/std at 50 mV  : ", np.mean(current_array[204383:209565], axis=0),np.std(current_array[204383:209565], axis=0))
#
# # ax_part1.set_xlim(0, 50)
# # ax_part1.set_ylim(-5, 5)
#     fig2 = plt.figure(1, figsize=(8, 4))
#     ax_part2 = fig1.add_subplot(2, 1, 2, sharex=ax_part1)
#     t = np.arange(len(dat['v1']))/dat['samplerate']
#     ax_part2.plot(t, dat['v1']*1e3)
# ax_part2.set_xlabel('Time (s)')
# ax_part2.set_ylabel('Voltage (mV)')
# # ax_part2.set_xlim(0, 50)
# # ax_part2.set_ylim(-60, 60)
# plt.show()


# # 2  Plot: Current-time (top), Voltage-time ON SAME graph
#
# for i,file in enumerate(filenames):
#     dat = LoadData.OpenFile(file)
#     fig1 = plt.figure(1, figsize=(8, 4))
#     ax_part1 = fig1.add_subplot(111)
#     t1 = np.arange(len(dat['i1']))/dat['samplerate']
#     color1 = 'tab:green'
#     ax_part1.plot(t1, dat['i1']*1e9, color = color1)
#     ax_part1.set_ylabel('Ionic Current (nA)')
#     ax_part1.set_xlabel('Time (s)')
#     ax_part1.tick_params(labelcolor=color1)
#
# for i,file in enumerate(filenames):
#     dat = LoadData.OpenFile(file)
#     fig1 = plt.figure(1, figsize=(8, 4))
#     ax_part2 = fig1.add_subplot(111, sharex=ax_part1)
#     t2 = np.arange(len(dat['v1']))/dat['samplerate']
#     color2 = 'tab:red'
#     ax_part2 = ax_part1.twinx()
#     ax_part2.plot(t2, dat['v1']*1e3, color=color2)
#     ax_part2.set_ylabel('Voltage (mV)')
#     ax_part2.tick_params( labelcolor=color2)
#
# plt.show()
# # ax_part2.set_xlim(0, 50)
# # ax_part2.set_ylim(-60, 60)


# # 3  Plot: Conductance-time (top), Voltage-time ON SAME graph
# for i,file in enumerate(filenames):
#     dat = LoadData.OpenFile(file)
#     fig1 = plt.figure(1, figsize=(8, 4))
#     ax_part1 = fig1.add_subplot(2, 1, 1)
#     t = np.arange(len(dat['i1']))/dat['samplerate']
#     ax_part1.plot(t, dat['i1']/dat['v1']*1e9)
#     # ax_part1.set_xlabel('Time (s)')
#     ax_part1.set_ylabel('Conductance (nS)')
#     fig2 = plt.figure(1, figsize=(8, 4))
#     ax_part2 = fig1.add_subplot(2, 1, 2, sharex=ax_part1)
#     t = np.arange(len(dat['v1']))/dat['samplerate']
#     ax_part2.plot(t, dat['v1']*1e3)
# ax_part2.set_xlabel('Time (s)')
# ax_part2.set_ylabel('Voltage (mV)')
# # # ax_part2.set_xlim(0, 50)
# # # ax_part2.set_ylim(-60, 60)
# plt.show()


# 4  Plot: Current-time (top) and Histogram-counts (sideways)
for i,file in enumerate(filenames):
    dat = LoadData.OpenFile(file)
    fig1 = plt.figure(1, figsize=(8, 4))
    ax_part1 = fig1.add_subplot(1,1,1)
    t = np.arange(len(dat['i1']))/dat['samplerate']
    ax_part1.plot(t, dat['i1'] * 1e9)


    fig2 = plt.figure(1, figsize=(1, 1))
    ax_part2 = fig1.add_subplot(1,2,2, sharey=ax_part1)
    current_array = np.array(dat['i1']*1e9)
    ax_part2 = sns.distplot(current_array, color="r",  norm_hist = True, bins = 3000, vertical = True) #bins = 3000,
    # ax_part2 = sns.kdeplot(current_array, color="r",  shade=True, vertical=True, kernel='gau') #bins = 3000,


    ax_part1.set_xlabel('Time (s)')
    ax_part1.set_ylabel('Ionic Current (nA)')
    # ax_part1.set_ylim(0, 1)
    # ax_part2.set_xlabel('Counts')
    plt.tight_layout()
    plt.show()
root.update()
