import MiscParameters as pm
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
from matplotlib import rc
rc('mathtext', default='regular')
pm.init(LoadFiles = 0)

currents = {'i1', 'i2'}
#filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file
filenames = {'/Users/migraf/SWITCHdrive/PhD/BPS/Poster Resources/Data/17B_10mMCis100mMtransKCl_80mer_2_OriginalDB.hdf5'}

count = 0
dt = {}
dI = {}
delay = np.array([], dtype=np.int64)
count = {}
ind = {}
fig2, ax3 = plt.subplots()

for k in currents:
    dt[k] = np.array([])
    dI[k] = np.array([])
    count[k] = 0
for file in filenames:
    for k in currents:
        print(file)
        filen = os.sep + str(os.path.split(file)[1][:-4])
        f = h5py.File(file, 'r')
        dt[k] = np.append(dt[k], f['LowPassSegmentation/' + k + '/DwellTime'].value)
        dI[k] = np.append(dI[k], f['LowPassSegmentation/' + k + '/FractionalCurrentDrop'].value)
        ind[k] = np.int64(f['LowPassSegmentation/' + k + '/' + pm.filter].value)
        count[k] += len(dt[k])
        sr = f['General/Samplerate'].value
    delay = f['XCorrelation/' + pm.filter + '/Lag'].value
    delay = delay[delay < 0.0004]
#    ax3.hist()
    ax3.hist(delay*1e6, 20,
                 label=r'Mean: {:03.2f}$\mu$s, std: {:03.2f}$\mu$s'.format(np.mean(delay*1e6), np.std(delay*1e6)))

## SCATTER PLOTS
fig1, ax = plt.subplots()
ax.plot(dt['i1'][ind['i1']] * 1e3, dI['i1'][ind['i1']]*100, 'ob', label = 'Ionic Current')
ax.set_ylabel(r'Ionic Current Drop $\frac{\Delta I}{I_0}$ [%]', color='b')
ax.tick_params('y', colors='b')
ax.set_xlabel('Time [ms]')
ax.set_title('MoS2 transverse detection, {} events'.format(count))

ax2 = ax.twinx()
ax2.plot(dt['i2'][ind['i2']] * 1e3, dI['i2'][ind['i2']]*100, 'or', label = 'Transverse Current')
ax2.set_ylabel(r'Transverse Current Drop $\frac{\Delta I}{I_0}$ [%]', color='r')
ax2.tick_params('y', colors='r')
#ax2.set_xlim(0, 1)


## Correlation Histograms
ax3.set_xlabel(r'Time [$\mu$s]')
ax3.set_ylabel('Counts')
ax3.set_title('Correlated events: {}, Histogram of delay'.format(len(delay)))
ax3.legend()

fig1.savefig(pm.OutputFolder + filen + 'Scatter.png', dpi=300)
fig1.savefig(pm.OutputFolder + filen + 'Scatter.eps')
fig2.savefig(pm.OutputFolder + filen + 'DelayHist.png', dpi=300)
fig2.savefig(pm.OutputFolder + filen + 'DelayHist.eps')

plt.show()
#fig1.clear()
#fig1.close()