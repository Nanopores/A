import numpy as np
import pandas as pd
import scipy.signal as sig
import Functions as uf
import os
import matplotlib.pyplot as plt
import matplotlib
import h5py
import seaborn as sns
import pickle
from tkinter.filedialog import askopenfilenames
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file

#If turned on reloads file and adds some baseline
AddBaseline = 0#10e-3

Name = '2kb, Roche Substrate, MoS2'
cm = plt.cm.get_cmap('RdYlBu')
voltageLimits = [0.11, 0.51]
dwellrange = (0, 10)
dIrange = (0, 7)
FracDIRange = (0, 95)

count = 0
dt = {}
dIF = {}
count = {}

dt = np.array([])
dIF = np.array([])
dI = np.array([])
v = np.array([])
t = np.array([])
t2 = np.array([])
count = 0
events=[]
print('filenames = [')
for file in filenames:
    print('\'{}\','.format(file))
    f = h5py.File(file, 'r')
    if AddBaseline:
        out = uf.OpenFile(f['General/FileName'].value)
    for i in range(0, len(f['LowPassSegmentation/' + 'i1' + '/DwellTime'].value)):
        if AddBaseline:
            s = f['LowPassSegmentation/' + 'i1' + '/StartPoints'].value[i]
            e = f['LowPassSegmentation/' + 'i1' + '/EndPoints'].value[i]
            events.append(out['i1'][int(s-out['samplerate']*AddBaseline):int(e+out['samplerate']*AddBaseline)]*1e9)
        else:
            events.append(f['LowPassSegmentation/' + 'i1' + '/AllEvents/{:09d}'.format(i)].value*1e9)
    dt = np.append(dt, f['LowPassSegmentation/' + 'i1' + '/DwellTime'].value*1e3)
    dIF = np.append(dIF, f['LowPassSegmentation/' + 'i1' + '/FractionalCurrentDrop'].value*100)
    dI = np.append(dI, f['LowPassSegmentation/' + 'i1' + '/DeltaI'].value*1e9)
    t = np.append(t, f['General/TimeFileWritten'].value*np.ones(len(f['LowPassSegmentation/' + 'i1' + '/DeltaI'].value)))
    t2 = np.append(t2, f['General/TimeFileLastModified'].value*np.ones(len(f['LowPassSegmentation/' + 'i1' + '/DeltaI'].value)))
    v = np.append(v, f['LowPassSegmentation/' + 'i1' + '/LocalVoltage'].value)
    count += len(f['LowPassSegmentation/' + 'i1' + '/DeltaI'].value)
print(']')

# Sort the voltages.
availableVoltages = np.unique(v[np.where((v > voltageLimits[0]) & (v < voltageLimits[1]))])
print(availableVoltages)
data = []
datadII0 = []
dataDwell = []
dataRate = []
labels = []
NEvents = []
allevents = []

for vo in availableVoltages:
    ind1 = np.argwhere(v == vo)
    ind2 = np.intersect1d(ind1, np.argwhere((dt > dwellrange[0]) & (dt < dwellrange[1])))
    ind = np.intersect1d(ind2, np.argwhere((dI > dIrange[0]) & (dI < dIrange[1])))
    data.append(dI[ind])
    datadII0.append(dIF[ind])
    dataDwell.append(dt[ind])
    NEvents.append(len(ind))
    all = np.array([])
    for j in ind:
        all = np.append(all, events[j])
    allevents.append(all)
    if f['General/Machine'].value == 'Axopatch':
        print('Modified-Written is applied')
        if (np.max(t2[ind])-np.min(t[ind])) == 0:
            dataRate.append(0)
        else:
            dataRate.append(len(ind)/(np.max(t2[ind])-np.min(t[ind])))
    else:
        if (np.max(t[ind])-np.min(t[ind])) == 0:
            dataRate.append(0)
        else:
            dataRate.append(len(ind)/(np.max(t[ind])-np.min(t[ind])))
    labels.append('{:0.0f}'.format(vo*1000))

cond = np.array([])
for i,dati in enumerate(allevents):
    cond = np.append(cond, dati/availableVoltages[i])

pickle.dump((allevents, availableVoltages), open(file[:-24] + 'AllEventPoints.p', 'wb'))

fig18 = plt.figure(18, figsize=(6, 9))
for i,dati in enumerate(allevents):
    # KDE DwellTime vs Voltages
    ax18 = fig18.add_subplot(len(data), 1, i + 1)
    ax18 = sns.distplot(dati/availableVoltages[i])
    ax18.set_xlabel('Conductance (nS)')
    ax18.set_ylabel('Frequency')
    ax18.set_title('{}mV'.format(int(availableVoltages[i] * 1e3)))
    ax18.set_xlim(left=np.min(cond), right=np.max(cond))
fig18.savefig(file[:-24] + 'KDEConductanceAllEventPoints.pdf')  # , transparent=True)




