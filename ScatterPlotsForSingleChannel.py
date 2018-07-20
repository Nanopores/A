import numpy as np
import scipy
import scipy.signal as sig
import Functions as f
import os
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import EngFormatter
import pandas as pd
import h5py
from tkinter.filedialog import askopenfilenames


filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file
#filenames = {'/Volumes/lben/lben-commun/2018 User Data/Martina/Axonpatch/20180518_13A/Results/13A_1MKCl_events250mV3_OriginalDB.hdf5'}

Name = '13A Crown in 1M KCl + something, Axopatch'
cm = plt.cm.get_cmap('RdYlBu')

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

for file in filenames:
    print(file)
    f = h5py.File(file, 'r')
    dt = np.append(dt, f['LowPassSegmentation/' + 'i1' + '/DwellTime'].value)
    dIF = np.append(dIF, f['LowPassSegmentation/' + 'i1' + '/FractionalCurrentDrop'].value)
    dI = np.append(dI, f['LowPassSegmentation/' + 'i1' + '/DeltaI'].value)
    t = np.append(t, f['General/TimeFileWritten'].value*np.ones(len(f['LowPassSegmentation/' + 'i1' + '/DeltaI'].value)))
    t2 = np.append(t2, f['General/TimeFileLastModified'].value*np.ones(len(f['LowPassSegmentation/' + 'i1' + '/DeltaI'].value)))
    v = np.append(v, f['LowPassSegmentation/' + 'i1' + '/LocalVoltage'].value)
    count += len(f['LowPassSegmentation/' + 'i1' + '/DeltaI'].value)

# Sort the voltages.
availableVoltages = np.unique(v[np.where(v > 0)])
print(availableVoltages)
data = []
datadII0 = []
dataDwell = []
dataRate = []
labels = []
NEvents = []
for vo in availableVoltages:
    ind = np.argwhere(v == vo)
    data.append(dI[ind])
    datadII0.append(dIF[ind])
    dataDwell.append(dt[ind])
    NEvents.append(len(ind))
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
    labels.append('{:0.0f}mV'.format(vo*1000))


# Fractional Current Drop Scatter Plot
fig1 = plt.figure(1, figsize=(9, 6))
ax = fig1.add_subplot(111)
sc = ax.scatter(dt[np.where(v > 0)], dIF[np.where(v > 0)]*100, c=v[np.where(v > 0)], vmin=min(v[np.where(v > 0)]), vmax=max(v[np.where(v > 0)]), s=35, cmap=cm)
cbar = plt.colorbar(sc, ticks=availableVoltages)
cbar.ax.set_yticklabels(labels)  # vertically oriented colorbar
ax.xaxis.set_major_formatter(EngFormatter(unit='s'))
ax.set_xlabel('Time')
ax.set_ylabel('Current Drop dI/I0 [%]')
ax.set_title('{}\nScatter Plot, {} events'.format(Name, count))
fig1.savefig(file[:-4] + 'ScatterFrac.png', dpi=300)
fig1.savefig(file[:-4] + 'ScatterFrac.eps')
ax.set_xlim(0, 50e-3)
ax.set_ylim(0, 30)
fig1.savefig(file[:-4] + 'ScatterFracZoomed.png', dpi=300)
fig1.savefig(file[:-4] + 'ScatterFracZoomed.eps')

#fig1.clear()
#fig1.close()

# Fractional Current Drop Scatter Plot
fig6 = plt.figure(6, figsize=(9, 6))
ax6 = fig6.add_subplot(111)
sc6 = ax6.scatter(dt[np.where(v > 0)], dI[np.where(v > 0)], c=v[np.where(v > 0)], vmin=min(v[np.where(v > 0)]), vmax=max(v[np.where(v > 0)]), s=35, cmap=cm)
cbar6 = plt.colorbar(sc6, ticks=availableVoltages)
cbar6.ax.set_yticklabels(labels)  # vertically oriented colorbar
ax6.xaxis.set_major_formatter(EngFormatter(unit='s'))
ax6.yaxis.set_major_formatter(EngFormatter(unit='A'))
ax6.set_xlabel('Time')
ax6.set_ylabel('Current Drop dI')
ax6.set_title('{}\nScatter Plot, {} events'.format(Name, count))
fig6.savefig(file[:-4] + 'ScatterdI.png', dpi=300)
fig6.savefig(file[:-4] + 'ScatterdI.eps')
ax6.set_xlim(0, 20e-3)
ax6.set_ylim(0, 2.5e-9)
fig6.savefig(file[:-4] + 'ScatterdIZoomed.png', dpi=300)
fig6.savefig(file[:-4] + 'ScatterdIZoomed.eps')

#fig1.clear()
#fig1.close()

# BoxPlot Delta I vs Voltages
fig2 = plt.figure(2, figsize=(9, 6))
ax2 = fig2.add_subplot(111)
bp = ax2.boxplot(data, notch=False, sym='')
ax2.yaxis.set_major_formatter(EngFormatter(unit='A'))
ax2.set_xlabel('Voltage')
ax2.set_ylabel('Current Drop dI')
ax2.set_title('{}\nBoxplot (outliers removed)\n{} events'.format(Name, count))
ax2.set_xticklabels(labels)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
fig2.savefig(file[:-4] + 'BoxplotdI.png', dpi=300)
fig2.savefig(file[:-4] + 'BoxplotdI.eps')

# BoxPlot DwellTime vs Voltages
fig3 = plt.figure(3, figsize=(9, 6))
ax3 = fig3.add_subplot(111)
bp = ax3.boxplot(dataDwell, notch=False, sym='')
ax3.yaxis.set_major_formatter(EngFormatter(unit='s'))
ax3.set_xlabel('Voltage')
ax3.set_ylabel('Dwell Time dt')
ax3.set_title('{}\nBoxplot (outliers removed)\n{} events'.format(Name, count))
ax3.set_xticklabels(labels)
ax3.get_xaxis().tick_bottom()
ax3.get_yaxis().tick_left()
fig3.savefig(file[:-4] + 'BoxplotDwellTime.png', dpi=300)
fig3.savefig(file[:-4] + 'BoxplotDwellTime.eps')

# BoxPlot Event Rate vs Voltages
fig4 = plt.figure(4, figsize=(9, 6))
ax4 = fig4.add_subplot(111)
bp = ax4.bar(np.arange(len(dataRate)), dataRate)
ax4.yaxis.set_major_formatter(EngFormatter(unit=r'$\frac{1}{s}$'))
ax4.set_xlabel('Voltage')
ax4.set_ylabel('Event Rate')
ax4.set_title('{}\Barplot\n{} events'.format(Name, count))
ax4.set_xticks(np.arange(len(dataRate)))
ax4.set_xticklabels(labels)
ax4.get_xaxis().tick_bottom()
ax4.get_yaxis().tick_left()
fig4.savefig(file[:-4] + 'BarplotEventRate.png', dpi=300)
fig4.savefig(file[:-4] + 'BarplotEventRate.eps')

# Bar Plot Event Number vs Voltages
fig5 = plt.figure(5, figsize=(9, 6))
ax5 = fig5.add_subplot(111)
bp = ax5.bar(np.arange(len(NEvents)), NEvents)
ax5.yaxis.set_major_formatter(EngFormatter(unit=''))
ax5.set_xlabel('Voltage')
ax5.set_ylabel('Number Of Events')
ax5.set_title('{}\nBarplot\n{} events'.format(Name, count))
ax5.set_xticks(np.arange(len(NEvents)))
ax5.set_xticklabels(labels)
ax5.get_xaxis().tick_bottom()
ax5.get_yaxis().tick_left()
fig5.savefig(file[:-4] + 'BarplotEventNumber.png', dpi=300)
fig5.savefig(file[:-4] + 'BarplotEventNumber.eps')

# BoxPlot Delta I vs Voltages
fig7 = plt.figure(7, figsize=(9, 6))
ax7 = fig7.add_subplot(111)
bp = ax7.boxplot(datadII0, notch=False, sym='')
ax7.set_xlabel('Voltage')
ax7.set_ylabel('Current Drop dI/I0 [%]')
ax7.set_title('{}\nBoxplot (outliers removed)\n{} events'.format(Name, count))
ax7.set_xticklabels(labels)
ax7.get_xaxis().tick_bottom()
ax7.get_yaxis().tick_left()
fig7.savefig(file[:-4] + 'BoxplotdII0.png', dpi=300)
fig7.savefig(file[:-4] + 'BoxplotdII0.eps')

plt.show()