
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
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

#1  Plot 211 and 212 type graph::: Current-time (top), Voltage-time (below)


# for i,file in enumerate(filenames):
#     dat = LoadData.OpenFile(file)
#     fig1 = plt.figure(1, figsize=(8, 4))
#     ax_part1 = fig1.add_subplot(2, 1, 1)
#     t = np.arange(len(dat['i1']))/dat['samplerate']
#     ax_part1.plot(t, dat['i1']*1e12)
#     # ax_part1.set_xlabel('Time (s)')
#     ax_part1.set_ylabel('Ionic Current (pA)')
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


#2  Plot: Current-time (top), Voltage-time ON SAME graph

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

# plt.show()
# ax_part2.set_xlim(0, 50)
# ax_part2.set_ylim(-60, 60)


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
    ax_part1 = fig1.add_subplot(1,2,1)
    t = np.arange(len(dat['i1']))/dat['samplerate']
    ax_part1.plot(t, dat['i1'] * 1e9)


    fig2 = plt.figure(1, figsize=(1, 1))
    ax_part2 = fig1.add_subplot(1,2,2, sharey=ax_part1)
    current_array = np.array(dat['i1']*1e9)
    ax_part2 = sns.distplot(current_array, color="r",  norm_hist = True, vertical = True)


    ax_part1.set_xlabel('Time (s)')
    ax_part1.set_ylabel('Ionic Current (nA)')
    # ax_part1.set_ylim(0, 1)
    ax_part2.set_xlabel('Counts')
    plt.tight_layout()
    plt.show()
#





# import os.path, time
#
# print("last modified: %s" % time.ctime(os.path.getmtime(file)))
#
# print("created: %s" % time.ctime(os.path.getctime(file)))
#
# plt.show()


# file = '/Volumes/lben-archives/2019 - CURRENT/Mukesh_Archives/Axopatch200B/09102019/100bp_1MKCl/Dev_34_1MKCl_100bp_500mV_1.dat'
#
# #file = '/Volumes/lben/lben-commun/2018 User Data/Michael/Axopatch/20180706/NorcadaCh1_100mMKCl_1mM_Cis_trans_ph74_640nm_150mmlens_50mW_IV_Later_Trace.dat'
#
# dat = LoadData.OpenFile(file)
#
# fig1 = plt.figure(1, figsize=(8, 4))
# ax_part1 = fig1.add_subplot(1, 1, 1)
#
# '''
# start2 = int(1.44e7)
# end2 = int(2.04e7)
# start1 = int(0)
# end1 = int(900000)
# i1 = np.array([])
# #i1 = np.append(i1, dat['i1'][start1:end1])
# #i1 = np.append(i1, dat['i1'][start2:end2]+(i1[-1:]-dat['i1'][start2]))
# '''
# t = np.arange(len(dat['i1']))/dat['samplerate']
# ax_part1.plot(t, dat['i1']*1e9)
#
# ax_part1.set(xlim=(0, 5), ylim=(-2, 2))
# #ax_part1.plot(dat['i1'])
#
# plt.show()
# #fig1.savefig('/Volumes/lben-archives/2019 - CURRENT/Mukesh_Archives/Axopatch200B/09102019/100bp_1MKCl/100 mV noise' + '.eps')