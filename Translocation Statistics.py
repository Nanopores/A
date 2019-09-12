#Michael's Code modified for for Translocation Kinetics on Jochem's Event Detection
import numpy as np
import pandas as pd
import scipy.signal as sig
import Functions as f
import shelve
from sklearn.mixture import GaussianMixture as GMM
import os
import matplotlib
#matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

import h5py
import seaborn as sns
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
root = Tk()
root.withdraw()
from scipy.stats import gaussian_kde


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

filenames = askopenfilenames(filetypes=[('data files', '*.dat')], parent = root)
root.update()


Name = ''
cm = plt.cm.get_cmap('RdYlBu')
voltageLimits = [0.19, 0.51]
dwellrange = (0, 25)
dIrange = (0, 2.5)
FracDIRange = (0, 100)

count = 0
dt = {}
dIF = {}
count = {}

dt = np.array([])
dIF = np.array([])
dI = np.array([])
localBase = np.array([])
v = np.array([])
t = np.array([])
t2 = np.array([])
count = 0

##This part needs to be adapted to Jochem-style files
for file in filenames:
    print(file)
    fileNamewithourend, end = os.path.splitext(file)
    shelfFile = shelve.open(fileNamewithourend)
    event = shelfFile['TranslocationEvents']
    dIF = np.append(dIF, np.array(event.GetAllIdropsNorm())*100)
    dI = np.append(dI, np.array(event.GetAllIdrops())*1e9)
    tempvoltage = np.array([])
    templocalbaseline = np.array([])
    tempeventlength = np.array([])
    for ev in event.events:
        tempvoltage = np.append(tempvoltage, ev.voltage)
        templocalbaseline = np.append(templocalbaseline, ev.baseline)
        tempeventlength = np.append(tempeventlength, ev.eventLength)
    filepath = event.events[0].filename[19:] #.replace('2019 - CURRENT', '2019')
    st = os.stat('/Volumes/'+filepath.replace('\\','/'))
    dt = np.append(dt, tempeventlength*1e3)
    t = np.append(t, st.st_birthtime * np.ones(len(event.GetAllIdropsNorm())))
    t2 = np.append(t2, st.st_mtime * np.ones(len(event.GetAllIdropsNorm())))
    v = np.append(v, tempvoltage)
    count += len(event.GetAllIdropsNorm())
    localBase = np.append(localBase, templocalbaseline*1e9/tempvoltage)
##

# Sort the voltages.
availableVoltages = np.unique(v[np.where((v > voltageLimits[0]) & (v < voltageLimits[1]))])
print(availableVoltages)
data = []
datadII0 = []
datalBase = []
dataDwell = []
dataRate = []
labels = []
NEvents = []
for vo in availableVoltages:
    ind1 = np.argwhere(v == vo)
    ind2 = np.intersect1d(ind1, np.argwhere((dt > dwellrange[0]) & (dt < dwellrange[1])))
    ind = np.intersect1d(ind2, np.argwhere((dI > dIrange[0]) & (dI < dIrange[1])))
    data.append(dI[ind])
    datadII0.append(dIF[ind])
    dataDwell.append(dt[ind])
    datalBase.append(localBase[ind])
    NEvents.append(len(ind))
    if 1:
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
    labels.append('{:0.0f}'.format(np.round(vo*1000)))

cond = np.array([])
for i,dati in enumerate(data):
    cond = np.append(cond, dati/availableVoltages[i])

# Fractional Current Drop % Scatter Plot
fig1 = plt.figure(1, figsize=(9, 6))
ax = fig1.add_subplot(111)
sc = ax.scatter(dt[np.where(v > 0)], dIF[np.where(v > 0)], c=v[np.where(v > 0)], vmin=min(v[np.where(v > 0)]), vmax=max(v[np.where(v > 0)]), s=35, cmap=cm)
cbar = plt.colorbar(sc, ticks=availableVoltages)
cbar.ax.set_yticklabels(labels)  # vertically oriented colorbar
ax.set_xlabel('Time (s)')
ax.set_ylabel('Current Drop dI/I0 (%)')
ax.set_title('{}\nScatter Plot, {} events'.format(Name, count))
#plt.show()
fig1.savefig(file[:-24] + 'ScatterFrac.pdf', transparent=True)
ax.set_xlim(dwellrange[0], dwellrange[1])
ax.set_ylim(FracDIRange[0], FracDIRange[1])
fig1.savefig(file[:-24] + 'ScatterFracZoomed.pdf', transparent=True)

#fig1.clear()
#fig1.close()

# 1.Fractional Current Drop Scatter Plot...
fig6 = plt.figure(6, figsize=(9, 6))
ax6 = fig6.add_subplot(111)
sc6 = ax6.scatter(dt[np.where(v > 0)], dI[np.where(v > 0)], c=v[np.where(v > 0)], vmin=min(v[np.where(v > 0)]), vmax=max(v[np.where(v > 0)]), s=35, cmap=cm)
cbar6 = plt.colorbar(sc6, ticks=availableVoltages)
cbar6.ax.set_yticklabels(labels)  # vertically oriented colorbar
ax6.set_xlabel('Time (s)')
ax6.set_ylabel('Current Drop dI (nA)')
ax6.set_title('{}\nScatter Plot, {} events'.format(Name, count))
#plt.show()
fig6.savefig(file[:-24] + 'ScatterdI.pdf', transparent=True)
ax6.set_xlim(dwellrange[0], dwellrange[1])
ax6.set_ylim(dIrange[0], dIrange[1])
fig6.savefig(file[:-24] + 'ScatterdIZoomed.pdf', transparent=True)

#plt.show()

#fig1.clear()
#fig1.close()

# 2.BoxPlot Delta I vs Voltages
fig2 = plt.figure(2, figsize=(9, 6))
ax2 = fig2.add_subplot(111)
bp = ax2.boxplot(data, notch=True, autorange=True, sym='', whis=1)
ax2.set_xlabel('Voltage (mV)')
ax2.set_ylabel('Current Drop dI (nA)')
ax2.set_title('{}\nBoxplot (outliers removed)\n{} events'.format(Name, count))
ax2.set_xticklabels(labels)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
fig2.savefig(file[:-24] + 'BoxplotdI.pdf', transparent=True)

# 3.BoxPlot DwellTime vs Voltages
fig3 = plt.figure(3, figsize=(9, 6))
ax3 = fig3.add_subplot(111)
bp = ax3.boxplot(dataDwell, notch=True, autorange=True, sym='', whis=1)
ax3.set_xlabel('Voltage (mV)')
ax3.set_ylabel('Dwell Time dt (ms)')
ax3.set_title('{}\nBoxplot (outliers removed)\n{} events'.format(Name, count))
ax3.set_xticklabels(labels)
ax3.get_xaxis().tick_bottom()
ax3.get_yaxis().tick_left()
fig3.savefig(file[:-24] + 'BoxplotDwellTime.pdf', transparent=True)

# 4.BoxPlot Event Rate vs Voltages
fig4 = plt.figure(4, figsize=(9, 6))
ax4 = fig4.add_subplot(111)
bp = ax4.bar(np.arange(len(dataRate)), dataRate)
ax4.set_xlabel('Voltage (mV)')
ax4.set_ylabel('Event Rate (Hz)')
ax4.set_title('{}\Barplot\n{} events'.format(Name, count))
ax4.set_xticks(np.arange(len(dataRate)))
ax4.set_xticklabels(labels)
ax4.get_xaxis().tick_bottom()
ax4.get_yaxis().tick_left()
fig4.savefig(file[:-24] + 'BarplotEventRate.pdf', transparent=True)

# 5.Bar Plot Event Number vs Voltages
fig5 = plt.figure(5, figsize=(9, 6))
ax5 = fig5.add_subplot(111)
print(NEvents)
bp = ax5.bar(np.arange(len(NEvents)), NEvents)
ax5.set_xlabel('Voltage (mV)')
ax5.set_ylabel('Number Of Events')
ax5.set_title('{}\nBarplot\n{} events'.format(Name, count))
ax5.set_xticks(np.arange(len(NEvents)))
ax5.set_xticklabels(labels)
ax5.get_xaxis().tick_bottom()
ax5.get_yaxis().tick_left()
fig5.savefig(file[:-24] + 'BarplotEventNumber.pdf', transparent=True)

# 7.BoxPlot Delta I vs Voltages
fig7 = plt.figure(7, figsize=(9, 6))
ax7 = fig7.add_subplot(111)
bp = ax7.boxplot(datadII0, notch=True, autorange=True, sym='', whis=1)
ax7.set_xlabel('Voltage (mV)')
ax7.set_ylabel('Current Drop dI/I0 (%)')
ax7.set_title('{}\nBoxplot (outliers removed)\n{} events'.format(Name, count))
ax7.set_xticklabels(labels)
ax7.get_xaxis().tick_bottom()
ax7.get_yaxis().tick_left()
fig7.savefig(file[:-24] + 'BoxplotdII0.pdf', transparent=True)


# 8.Histogram DwellTime vs Voltages
fig8 = plt.figure(8, figsize=(6, 12))
for i,dati in enumerate(dataDwell):
    ax8 = fig8.add_subplot(len(dataDwell), 1, i+1)
    n, bins, patches = ax8.hist(x=dati, bins = 100, alpha=0.7, rwidth=0.85)
    ax8.grid(axis='y', alpha=0.75)
    ax8.set_xlabel('Dwell Time (ms)')
    ax8.set_ylabel('Frequency')
    ax8.set_title('{}mV'.format(int(availableVoltages[i]*1e3)))
    #ax8.set_xlim(left=np.min(np.concatenate(dataDwell).ravel()), right=np.max(np.concatenate(dataDwell).ravel()))
    ax8.set_xlim(0,2)
fig8.savefig(file[:-24] + 'HistogramDwellTime.pdf', transparent=True)

# 9.Histogram DwellTime vs Voltages
fig9 = plt.figure(9, figsize=(6, 12))
for i,dati in enumerate(data):
    ax9 = fig9.add_subplot(len(data), 1, i+1)
    n, bins, patches = ax9.hist(x=dati, bins = 100, alpha=0.7, rwidth=0.85)
    ax9.grid(axis='y', alpha=0.75)
    ax9.set_xlabel('Current Drop (nA)')
    ax9.set_ylabel('Frequency')
    ax9.set_title('{}mV'.format(int(availableVoltages[i]*1e3)))
    #ax9.set_xlim(left=np.min(np.concatenate(data).ravel()), right=np.max(np.concatenate(data).ravel()))
    ax9.set_xlim(0, 2)
fig9.savefig(file[:-24] + 'HistogramCurrentDrop.pdf', transparent=True)
plt.style.use("seaborn-darkgrid")

# 10 KDE DwellTime vs Voltages
fig10 = plt.figure(10, figsize=(6, 12))
for i, dati in enumerate(dataDwell):
    ax10 = fig10.add_subplot(len(dataDwell), 1, i+1)
    ax10 = sns.distplot(dati, rug= False)
    ax10.set_xlabel('Dwell Time (ms)')
    ax10.set_ylabel('Frequency')
    ax10.set_title('{}mV'.format(int(availableVoltages[i]*1e3)))
    ax10.set_xlim(left=np.min(np.concatenate(dataDwell).ravel()), right=np.max(np.concatenate(dataDwell).ravel()))
fig10.savefig(file[:-24] + 'KDEDwellTime.pdf')#, transparent=True)

# 11 KDE DwellTime vs Voltages
fig11 = plt.figure(11, figsize=(6, 12))
for i,dati in enumerate(data):
    ax11 = fig11.add_subplot(len(data), 1, i+1)
    ax11 = sns.distplot(dati)
    ax11.set_xlabel('Current Drop (nA)')
    ax11.set_ylabel('Frequency')
    # ax11.set_xlim(0, 2)
    ax11.set_title('{}mV'.format(int(availableVoltages[i]*1e3)))
    ax11.set_xlim(left=np.min(np.concatenate(data).ravel()), right=np.max(np.concatenate(data).ravel()))
    ax8.set_xlim(0, 2)
fig11.savefig(file[:-24] + 'KDECurrentDrop.pdf')#, transparent=True)

# 17 KDE DwellTime vs Voltages
fig17 = plt.figure(17, figsize=(6, 12))
for i,dati in enumerate(datadII0):
    ax17 = fig17.add_subplot(len(datadII0), 1, i+1)
    ax17 = sns.distplot(dati)
    ax17.set_xlabel('Current Drop (nA)')
    ax17.set_ylabel('Frequency')
    ax17.set_title('{}mV'.format(int(availableVoltages[i]*1e3)))
    ax17.set_xlim(left=0, right=100)
fig17.savefig(file[:-24] + 'KDECurrentDropdII0.pdf')#, transparent=True)

# 18 KDE DwellTime vs Voltages
fig18 = plt.figure(18, figsize=(6, 12))
for i, dati in enumerate(data):
    ax18 = fig18.add_subplot(len(data), 1, i+1)
    ax18 = sns.distplot(dati/availableVoltages[i])
    ax18.set_xlabel('Conductance Drop (nS)')
    ax18.set_ylabel('Frequency')
    ax18.set_title('{}mV'.format(int(availableVoltages[i]*1e3)))
    ax18.set_xlim(left=np.min(cond), right=np.max(cond))
fig18.savefig(file[:-24] + 'KDEConductanceDrop.pdf')#, transparent=True)

# 15 BiVariate KDE plot Individual
fig15 = plt.figure(15, figsize=(9, 9))
for i,dati in enumerate(data):
    ax15 = fig15.add_subplot(2, 2, i + 1)
    # cmap = sns.cubehelix_palette(as_cmap=True, dark=0, light=1, reverse=True)
    ax15 = sns.kdeplot(dataDwell[i]*1e-3, dati*1e-9, shade=True, cbar=True, cut=0, cmap='rainbow', n_levels=60, reverse=False, shade_lowest=True, label='{}mV'.format(str(int(availableVoltages[i]*1e3))))
    ax15.set_ylabel('Current Drop (nA)')
    ax15.set_xlabel('Dwell Time (ms)')
    ax15.set_title('{}mV'.format(int(availableVoltages[i]*1e3)))
    ax15.set_xlim(left=np.min(np.concatenate(dataDwell).ravel())-1, right=np.max(np.concatenate(dataDwell).ravel())+1)
    ax15.set_ylim(bottom=np.min(np.concatenate(data).ravel())-1, top=np.max(np.concatenate(data).ravel()))
fig15.savefig(file[:-24] + 'KDEScatterIndividualdI.pdf')#, transparent=True)

# BiVariate KDE plot Individual Fractional
fig16 = plt.figure(16, figsize=(9, 9))
for i,dati in enumerate(datadII0):
    ax16 = fig16.add_subplot(2, 2, i + 1)
    ax16 = sns.kdeplot(dataDwell[i], dati, shade=True, cbar=True, cmap='gist_heat', shade_lowest=False, label='{}mV'.format(str(int(availableVoltages[i]*1e3))))
    ax16.set_ylabel('Fractional Current Drop (dI/I0)')
    ax16.set_xlabel('Dwell Time (ms)')
    ax16.set_title('{}mV'.format(int(availableVoltages[i]*1e3)))
    ax16.set_xlim(left=np.min(np.concatenate(dataDwell).ravel())-1, right=np.max(np.concatenate(dataDwell).ravel())+1)
    #ax16.set_ylim(bottom=0, top=100)
    # Set graph limits manually
    # ax16.set_xlim(0, 8)
    # ax16.set_ylim(0, 8)
fig16.savefig(file[:-24] + 'KDEScatterIndividualdIdI0.pdf')#, transparent=True)

# BiVariate KDE plot Individual
fig19 = plt.figure(19, figsize=(9, 9))
for i,dati in enumerate(data):
    ax19 = fig19.add_subplot(2, 2, i + 1)
    ax19 = sns.kdeplot(dataDwell[i], dati/availableVoltages[i], cbar=True, cmap='gist_heat', shade=True, shade_lowest=False, label='{}mV'.format(str(int(availableVoltages[i]*1e3))))
    ax19.set_ylabel('Conductance Drop (nS)')
    ax19.set_xlabel('Dwell Time (ms)')
    ax19.set_title('{}mV'.format(int(availableVoltages[i]*1e3)))
    ax19.set_xlim(left=np.min(np.concatenate(dataDwell).ravel())-1, right=np.max(np.concatenate(dataDwell).ravel())+1)
    ax19.set_ylim(bottom=np.min(cond)-1, top=np.max(cond))
    # Set graph limits manually
    # ax19.set_xlim(0, )
    # ax19.set_ylim(0, )
fig19.savefig(file[:-24] + 'KDEScatterConductanceDrop.pdf')#, transparent=True)

matplotlib.rcParams['font.size'] = 22

# KDE DwellTime vs Voltages ALL on One
fig12 = plt.figure(12, figsize=(16, 9))
ax12 = fig12.add_subplot(111)
for i, dati in enumerate(dataDwell):
    ax12 = sns.distplot(dati, label='{}mV'.format(str(int(availableVoltages[i]*1e3))))
ax12.legend()
ax12.set_xlabel('Dwell Time (ms)')
ax12.set_ylabel('Frequency')
ax12.set_title('Dwell Time Histograms')
#ax12.set_xlim(left=np.min(np.concatenate(dataDwell).ravel()), right=np.max(np.concatenate(dataDwell).ravel()))
ax12.set_xlim(0, 6)
fig12.savefig(file[:-24] + 'KDESingleDwellTime.pdf')#, transparent=True)

# KDE DwellTime vs Voltages ALL on One
fig13 = plt.figure(13, figsize=(16, 9))
ax13 = fig13.add_subplot(111)
for i,dati in enumerate(data):
    ax13 = sns.distplot(dati, label='{}mV'.format(str(int(availableVoltages[i]*1e3)), norm_hist=True, hist=True))
ax13.legend()
ax13.set_xlabel('Current Drop (nA)')
ax13.set_ylabel('Frequency')
ax13.set_title('Current Drop Histograms')
ax13.set_xlim(left=np.min(np.concatenate(data).ravel()), right=np.max(np.concatenate(data).ravel()))
fig13.savefig(file[:-24] + 'KDESingleCurrentDrop.pdf')

# KDE DwellTime vs Voltages ALL on One
fig20 = plt.figure(20, figsize=(16, 9))
ax20 = fig20.add_subplot(111)
for i,dati in enumerate(data):
    ax20 = sns.distplot(dati/availableVoltages[i], label='{}mV'.format(str(int(availableVoltages[i]*1e3))))
ax20.legend()
ax20.set_xlabel('Conductance Drop (nS)')
ax20.set_ylabel('Frequency')
ax20.set_title('Conductance Drop Histograms')
ax20.set_xlim(left=np.min(cond), right=np.max(cond))
fig20.savefig(file[:-24] + 'KDESingleConductanceDrop.pdf')#, transparent=True)

# BiVariate KDE plot All In One
fig14 = plt.figure(14, figsize=(9, 9))
ax14 = fig14.add_subplot(111)
for i,dati in enumerate(data):
    ax14 = sns.kdeplot(dataDwell[i], dati, shade=True, shade_lowest=False, label='{}mV'.format(str(int(availableVoltages[i]*1e3))))
ax14.legend()
ax14.set_ylabel('Current Drop (nA)')
ax14.set_xlabel('Dwell Time (ms)')
ax14.set_title('KDE Scatter')
# Set graph limits
# ax14.set_xlim(0, 8)
# ax14.set_ylim(0, 5)
#ax14.set_xlim(left=np.min(np.concatenate(data).ravel()), right=np.max(np.concatenate(data).ravel()))
fig14.savefig(file[:-24] + 'KDEScatterdI.pdf')#, transparent=True)

# BiVariate KDE plot All In One
fig21 = plt.figure(21, figsize=(9, 9))
ax21 = fig21.add_subplot(111)
for i,dati in enumerate(data):
    ax21 = sns.kdeplot(dataDwell[i], dati/availableVoltages[i], shade=True, shade_lowest=False, label='{}mV'.format(str(int(availableVoltages[i]*1e3))))
ax21.legend()
ax21.set_ylabel('Conductance Drop (nA)')
ax21.set_xlabel('Dwell Time (ms)')
ax21.set_title('KDE Scatter')
# ax21.set_xlim(0, 8)
ax21.set_ylim(0, 5)
#ax21.set_xlim(left=np.min(np.concatenate(data).ravel()), right=np.max(np.concatenate(data).ravel()))
fig21.savefig(file[:-24] + 'KDEScatterConductance.pdf')#, transparent=True)

# KDE DwellTime vs Voltages
fig22 = plt.figure(22, figsize=(6, 12))
for i, dati in enumerate(datalBase):
    ax22 = fig22.add_subplot(len(data), 1, i+1)
    ax22 = sns.distplot(dati)
    ax22.set_xlabel('Conductance (nS)')
    ax22.set_ylabel('Frequency')
    # ax22.set_xlim(0, 50)
    ax22.set_title('{}mV'.format(int(availableVoltages[i]*1e3)))
    ax22.set_xlim(left= 0, right=50)#np.min(np.concatenate(datalBase).ravel()), right=np.max(np.concatenate(datalBase).ravel()))
fig22.savefig(file[:-24] + 'LocalBaseline.pdf')#, transparent=True)

