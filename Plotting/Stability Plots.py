import platform

from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")

import matplotlib
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
import os
import platform
import os, time
# np.seterr(divide='ignore', invalid='ignore')
root = Tk()
root.withdraw()

#Set units formatting
from matplotlib.ticker import EngFormatter

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
formCurrent = EngFormatter(unit='A', places=3)
formA = EngFormatter(unit='A', places=3)
formB = EngFormatter(unit='s', places=3)

filenames = askopenfilenames()
root.update()


time_array_long = []
current_array_long = []
voltage_array_long = []
conductance_array_long = []

last_el = 0

for i, file in enumerate(filenames):
    dat = LoadData.OpenFile(file)

    fig1 = plt.figure(1, figsize=(8, 4))
    ax_part1 = fig1.add_subplot(2, 1, 1)

    #time data
    time_array = np.arange(len(dat['i1']))/dat['samplerate']

    #Ionic current
    current_array = np.array(dat['i1'] * 1e9)
    current_array_long = np.append([current_array_long],[[current_array]])

    # #Conductance (G)
    G = dat['i1']/dat['v1']*1e9
    conductance_array = np.array(G)
    conductance_array_long = np.append([conductance_array_long],[[conductance_array]])

    #Voltage
    voltage_array = np.array(len(dat['v1'] * 1e3))
    voltage_array_long = np.append([voltage_array_long], [[voltage_array]])

   #fix the time after joining
    time_interval = time_array[1] - time_array[0]
    time_array_long =  np.concatenate((time_array_long, time_array+last_el+time_interval))
    last_el = time_array_long[-1]

    fig2 = plt.figure(1, figsize=(8, 4))
    ax_part2 = fig2.add_subplot(2, 1, 2, sharex=ax_part1)



time_array_long -= time_interval

# clean_data = np.where(
#     ((time_array_long < 5) & (time_array_long > 0))|
#     ((time_array_long < 15) & (time_array_long > 10))
#     )

# # time_array_long[clean_data] = np.NaN #specify ranges that are masked and not plotted
# time_array_long [int(3e5):int(5e5)] = np.NaN
# time_array_long [int(7e5):int(9e5)] = np.NaN
# time_array_long [int(11e5):int(13e5)] = np.NaN
# time_array_long [int(15e5):int(17e5)] = np.NaN
# time_array_long [int(19e5):int(21e5)] = np.NaN
# time_array_long [int(23e5):int(25e5)] = np.NaN
# time_array_long [int(27e5):int(29e5)] = np.NaN
# time_array_long [int(31e5):int(33e5)] = np.NaN
# time_array_long [int(35e5):int(37e5)] = np.NaN

print(len(time_array_long),
len(current_array_long),
len(voltage_array_long),
len(conductance_array_long))
print('Done')

ax_part1.plot(time_array_long[::], conductance_array_long [::])
ax_part2.plot(time_array_long[::], current_array_long[::])

plt.xlabel ('Time(s)')
# plt.ylabel ('Ionic current (nA)')

#plt.savefig("joy25p67m_stability_plot.svg")
# plt.savefig( "joy25p67m_stability_plot.png")

plt.show()
root.update()

# print("last modified: %s" % time.ctime(os.path.getmtime(file)))
# print("created: %s" % time.ctime(os.path.getctime(file)))