import numpy as np
import h5py
import Functions as uf
import os
import matplotlib.pyplot as plt
import matplotlib
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
import MiscParameters as pm
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
pm.init()

#On the same plot?
sameplot = True
toPlot = 'i1'
baselinetime = 10e-3
timeconvert = 1000 #ms

filenames = [['/Users/migraf/Desktop/Temporary Analysis/Roche_Pdms53_1MKCl_2kb_1_OriginalDB.hdf5', 58],
             ['/Users/migraf/Desktop/Temporary Analysis/Roche_Pdms53_1MKCl_2kb_200mV_1_OriginalDB.hdf5', 131],
             ['/Users/migraf/Desktop/Temporary Analysis/Roche_Pdms53_1MKCl_2kb_400mV_2_OriginalDB.hdf5', 846],
             ['/Users/migraf/Desktop/Temporary Analysis/Roche_Pdms53_1MKCl_2kb_500mV_1_OriginalDB.hdf5', 2788]]

if sameplot:
    fig = plt.figure(1, figsize=(6, 5))
    ax = fig.add_subplot(111)
    ax.set_ylabel('Current')
    ax.set_xlabel('Time')

for file in filenames:
    f = h5py.File(file[0], 'r')
    out = uf.OpenFile(f['General/FileName'].value)
    baseline = np.uint64(baselinetime * out['samplerate'])
    filena = str(os.path.split(file[0])[1][:-5])
    startp_i1 = np.int64(
        f['LowPassSegmentation/' + toPlot + '/' + 'StartPoints'].value[file[1]] - baseline)
    endp_i1 = np.int64(
        f['LowPassSegmentation/' + toPlot + '/' + 'EndPoints'].value[file[1]] + baseline)
    if startp_i1 < 2 * baseline:
        continue
    if endp_i1 > len(out['i1']):
        endp_i1 = len(out['i1'])
    ax.plot(np.arange(endp_i1 - startp_i1) / out['samplerate'] * timeconvert, out['i1'][startp_i1:endp_i1] * 1e9, label = filena + '_' + str(file[1]))

ax.legend()
fig.savefig(pm.OutputFolder + filena + '_IndividualEventPlot.pdf', transparent=True)
plt.show()