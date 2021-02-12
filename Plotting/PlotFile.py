import matplotlib
import platform
if platform.system() == 'Darwin':
    matplotlib.use('MacOSX')
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
formCurrent = EngFormatter(unit='A', places=3)
formA = EngFormatter(unit='A', places=3)
formB = EngFormatter(unit='s', places=3)

filenames = askopenfilenames()


#1 Current (top), Voltage (below)

# for i,file in enumerate(filenames):
#     dat = LoadData.OpenFile(file)
#     fig1 = plt.figure(1, figsize=(8, 4))
#     ax_part1 = fig1.add_subplot(2, 1, 1)
#
#     t = np.arange(len(dat['i1']))/dat['samplerate']
#     # ax_part1.plot(t, dat['i1'] * 1e9)
#
#     ax_part1.set_ylabel('Ionic Current (nA)')
#
#     current_array = np.array(dat['i1'] * 1e9)
#     ax_part1.plot(t, current_array)
#
#     voltage_array = np.array(len(dat['v1']))
#
#     conductance_array = (current_array/voltage_array)
#
#     print(np.mean(conductance_array))
#     print(filenames)
#
# # # ax_part1.set_xlim(0, 50)
# # # ax_part1.set_ylim(-5, 5)
#     fig2 = plt.figure(1, figsize=(8, 4))
#     ax_part2 = fig1.add_subplot(2, 1, 2, sharex=ax_part1)
#     t = np.arange(len(dat['v1']))/dat['samplerate']
#     ax_part2.plot(t, dat['v1']*1e3)
#     ax_part2.set_xlabel('Time (s)')
#     ax_part2.set_ylabel('Voltage (mV)')
# # ax_part2.set_xlim(0, 50)      §
# # ax_part2.set_ylim(-60, 60)
#
# plt.show()

# os.chdir(os.path.dirname(file))
# directory = (str(os.path.split(file)[0]) + os.sep + 'Extracted' + '_traces')
# if not os.path.exists(directory):
#     os.makedirs(directory)
#
# fig.savefig(directory + os.sep + filename)
#
# plt.savefig("test.svg")


# # 2 Current-time (top), Voltage-time (below), Histogram (right)

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
# # # ax_part1.set_xlim(0, 50)
# # # ax_part1.set_ylim(-5, 5)
#     fig2 = plt.figure(1, figsize=(8, 4))
#     ax_part2 = fig1.add_subplot(2, 1, 2, sharex=ax_part1)
#     t = np.arange(len(dat['v1']))/dat['samplerate']
#     ax_part2.plot(t, dat['v1']*1e3)
# # ax_part2.set_xlabel('Time (s)')
# # ax_part2.set_ylabel('Voltage (mV)')
# # # ax_part2.set_xlim(0, 50)      §
# # # ax_part2.set_ylim(-60, 60)
# #
# fig3 = plt.figure(1, figsize = (1,1))
# ax_part3 = fig1.add_subplot(2,2,2, sharey=ax_part1)
# current_array = np.array(dat['i1']*1e9)
# ax_part3 = sns.distplot(current_array, color="r",  norm_hist = True, vertical = True) #bins = 3000,
# plt.show()



# 3  Plot: Current-time (top), Voltage-time ON SAME graph
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
# print(file)
# plt.show()
# # ax_part2.set_xlim(0, 50)
# # ax_part2.set_ylim(-60, 60)


# # 4  Plot: Conductance or Pore Size (top), Voltage-time ON SAME graph

for i,file in enumerate(filenames):
    dat = LoadData.OpenFile(file)

    fig1 = plt.figure(1, figsize=(10, 6))
    ax_part1 = fig1.add_subplot(2, 1, 1)
    t = np.arange(len(dat['i1']))/dat['samplerate']
    # t2 = t/60

    G = dat['i1']/dat['v1']*1e9

    # size = (G+np.sqrt(G*(G+16*10.5/np.pi)))/(2*10.5) #change the bulk conductivity and thickness
    G_mean = np.mean(G)
    G_std = np.std(G)
    print(G_mean, G_std)

    ax_part1.plot(t[::], G[::], label = G_mean)  #dat['i1']/dat['v1']*1e9 or #size
    ax_part1.legend()
    ax_part1.set_xlabel('Time (s)')
    ax_part1.set_ylabel('Conductance (nS)')


    fig2 = plt.figure(1, figsize=(10, 6))
    ax_part2 = fig1.add_subplot(2, 1, 2, sharex=ax_part1)
    ax_part2.plot(t[::], dat['v1'][::])
    ax_part2.set_xlabel('Time (s)')
    ax_part2.set_ylabel('Voltage (V)')
    # # ax_part2.set_xlim(0, 50)
    # ax_part1.set_ylim(0, 5)
    # ax_part2.set_ylim(0, 0.5)

plt.show()
# plt.savefig("Stability_300mV.svg")

#5  Plot: Current-time (top) and Histogram-counts (sideways)
# for i,file in enumerate(filenames):
#     dat = LoadData.OpenFile(file)
#     fig1 = plt.figure(1, figsize=(8, 4))
#     ax_part1 = fig1.add_subplot(1,2,1)
#     t = np.arange(len(dat['i1']))/dat['samplerate']
#     ax_part1.plot(t, dat['i1'] * 1e9)
#
#
#     fig2 = plt.figure(1, figsize=(1, 1))
#     ax_part2 = fig1.add_subplot(1,2,2, sharey=ax_part1)
#     current_array = np.array(dat['i1']*1e9)
#     ax_part2 = sns.distplot(current_array, color="r",bins = 1000,  norm_hist = True, vertical = True) #bins = 3000,
#     # ax_part2 = sns.kdeplot(current_array, color="r",  shade=True, vertical=True, kernel='gau') #bins = 3000,
#
#
#     ax_part1.set_xlabel('Time (s)')
#     ax_part1.set_ylabel('Ionic Current (nA)')
#     ax_part1.set_ylim(0, 1.2)
#     # ax_part1.set_xlim(0, 11.8)
#     ax_part2.set_xlabel('Counts (a.u.)')
#     # plt.tight_layout()
#     plt.show()
# root.update()

# for i,file in enumerate(filenames):
#     dat = LoadData.ImportABF(file)
#     fig1 = plt.figure(1, figsize=(8, 4))
#     ax_part1 = fig1.add_subplot(2, 1, 1)
#
#     t = np.arange(len(dat['i1']))/dat['samplerate']
#     # ax_part1.plot(t, dat['i1'] * 1e9)
#
#     ax_part1.set_ylabel('Ionic Current (nA)')
#
#     current_array = np.array(dat['i1'] * 1e9)
#     ax_part1.plot(t, current_array)
#
#     voltage_array = np.array(len(dat['v1']))
#
#     conductance_array = (current_array/voltage_array)
#     print(np.mean(conductance_array))
#
# # # ax_part1.set_xlim(0, 50)
# # # ax_part1.set_ylim(-5, 5)
#     fig2 = plt.figure(1, figsize=(8, 4))
#     ax_part2 = fig1.add_subplot(2, 1, 2, sharex=ax_part1)
#     t = np.arange(len(dat['v1']))/dat['samplerate']
#     ax_part2.plot(t, dat['v1']*1e3)
#     ax_part2.set_xlabel('Time (s)')
#     ax_part2.set_ylabel('Voltage (mV)')
# # ax_part2.set_xlim(0, 50)      §
# # ax_part2.set_ylim(-60, 60)

# plt.show()

# # 6  Plot: Conductance - Time and Histogram-counts (sideways)
#
# for i,file in enumerate(filenames):
#     dat = LoadData.OpenFile(file)
#     fig1 = plt.figure(1, figsize=(8, 4))
#     ax_part1 = fig1.add_subplot(1,1,1)
#     t = np.arange(len(dat['i1']))/dat['samplerate']
#     v = np.arange(len(dat['v1']))/dat['samplerate']
#     G = (dat['i1'] * 1e9/ v )
#     ax_part1.plot(t, G)
#
#
#     fig2 = plt.figure(1, figsize=(1, 1))
#     ax_part2 = fig1.add_subplot(1,2,2, sharey=ax_part1)
#     current_array = np.array(dat['i1']*1e9)
#     ax_part2 = sns.distplot(current_array, color="r",bins = 2000,  norm_hist = True, vertical = True) #bins = 3000,
#     # ax_part2 = sns.kdeplot(current_array, color="r",  shade=True, vertical=True, kernel='gau') #bins = 3000,
#
#
#     ax_part1.set_xlabel('Time (s)')
#     ax_part1.set_ylabel('Conductance (nS)')
#     # ax_part1.set_ylim(0, 1)
#     # ax_part1.set_xlim(0, 11.8)
#     # ax_part2.set_xlabel('Counts')
#     # plt.tight_layout()
#
# plt.show()


# --------------------------------------------------------------------------------------
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