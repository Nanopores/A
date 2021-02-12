import AnalysisParameters as pm
import numpy as np
import scipy
import scipy.signal as sig
import Functions as f
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
import pandas as pd
import h5py
from tkinter.filedialog import askopenfilenames
pm.init()

#currents = {'i1','i2'}
currents = {'i1'}

filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file
# filenames = {'/Users/migraf/Desktop/Chan/Data/R43_KCl gradient_100mMcis_1mM trans_4_20180405_182358_OriginalDB.hdf5'}
count = 0

dt = {}
dI = {}
count = {}
ind = {}
for k in currents:
    dt[k] = np.array([])
    dI[k] = np.array([])
    count[k] = 0

for k in currents:
    for file in filenames:
        print(file)
        f = h5py.File(file, 'r')
        dt[k] = np.append(dt[k], f['LowPassSegmentation/' + k + '/DwellTime'].value)
        dI[k] = np.append(dI[k], f['LowPassSegmentation/' + k + '/FractionalCurrentDrop'].value)
        ind[k] = f['LowPassSegmentation/' + k + '/CommonIndex'].value
        count[k] += len(dt[k])

fig1, ax = plt.subplots(1)
ax.plot(dt['i1'] * 1e3, dI['i1']*100, 'ob', label = 'Ionic Current')
ax.plot(dt['i2'] * 1e3, dI['i2']*100, 'or', label = 'Transverse Current')

#ax.set_xlim(0, 2)
#ax.set_ylim(10, 40)

ax.set_xlabel('Time [ms]')
ax.set_ylabel('Current Drop dI/I0 [%]')
ax.set_title('MoS2 transverse detection, {} events'.format(count))

fig1.savefig(file[:-4] + 'Scatter.png', dpi=300)
fig1.savefig(file[:-4] + 'Scatter.eps')

plt.show()

#fig1.clear()
#fig1.close()