## Cross-correlation of two channels
import AnalysisParameters as pm
import numpy as np
import scipy
import scipy.signal as sig
import Functions as fu
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
import pandas as pd
#import MiscParameters as pm
from matplotlib.ticker import EngFormatter
import h5py
from tkinter.filedialog import askopenfilenames
from matplotlib import rc
rc('mathtext', default='regular')
pm.init(LoadFiles = 0)

baselinetime = 10e-3
fitOn = 1
percentage = 1
plotPoints=0
AllInOneFile = 0

#filenames = {'/Volumes/lben/lben-commun/2018 User Data/Martina/Axonpatch/20180518_13A/Results/13A_1MKCl_events250mV3_OriginalDB.hdf5'}
filenames = askopenfilenames()

if AllInOneFile:
    pp = PdfPages(pm.OutputFolder + 'All_EventPlots.pdf')
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax2.set_title('Ionic Voltage')
    ax1.set_ylabel('Current')
    ax2.set_ylabel('Voltage')
    ax1.set_xlabel('Time')
    ax2.set_xlabel('Time')
    ax1.yaxis.set_major_formatter(EngFormatter(unit='A', places=1))
    ax2.yaxis.set_major_formatter(EngFormatter(unit='V', places=1))
    ax1.xaxis.set_major_formatter(EngFormatter(unit='s', places=1))
    ax2.xaxis.set_major_formatter(EngFormatter(unit='s', places=1))
    line1, = ax1.plot([])
    line2, = ax2.plot([], 'c')
    if fitOn:
        line1f, = ax1.plot([], 'r', ls='dashed')

for filename in filenames:
    print(filename)
    f = h5py.File(filename, 'r')
    out = fu.OpenFile(f['General/FileName'].value)
    baseline = np.uint64(baselinetime * out['samplerate'])
    file = os.sep + str(os.path.split(filename)[1][:-5])
    toPlot = ['i1']
    for k in toPlot:
        NumberOfEvents=f['LowPassSegmentation/' + k + '/NumberOfEvents'].value
        if NumberOfEvents is not 0:
            if not AllInOneFile:
                pp = PdfPages(pm.OutputFolder + file + '_Only_' + k + '_EventPlots.pdf')
                fig = plt.figure()
                ax1 = fig.add_subplot(211)
                ax2 = fig.add_subplot(212, sharex=ax1)
                ax2.set_title('Ionic Voltage')
                ax1.set_ylabel('Current')
                ax2.set_ylabel('Voltage')
                ax1.set_xlabel('Time')
                ax2.set_xlabel('Time')
                ax1.yaxis.set_major_formatter(EngFormatter(unit='A', places=1))
                ax2.yaxis.set_major_formatter(EngFormatter(unit='V', places=1))
                ax1.xaxis.set_major_formatter(EngFormatter(unit='s', places=1))
                ax2.xaxis.set_major_formatter(EngFormatter(unit='s', places=1))
                line1, = ax1.plot([])
                line2, = ax2.plot([], 'c')
                if fitOn:
                    line1f, = ax1.plot([], 'r', ls='dashed')
            for i in np.linspace(0, NumberOfEvents-1, NumberOfEvents*percentage, dtype=np.uint64):
                startp_i1 = np.int64(
                    f['LowPassSegmentation/' + k + '/' + 'StartPoints'].value[i] - baseline)
                endp_i1 = np.int64(
                    f['LowPassSegmentation/' + k + '/' + 'EndPoints'].value[i] + baseline)
                if startp_i1 < 2*baseline:
                    continue
                if endp_i1 > len(out['i1']):
                    endp_i1 = len(out['i1'])
                if fitOn:
                    if plotPoints:
                        fit1 = np.concatenate(([np.ones(baseline) *
                                               f['LowPassSegmentation/' + k + '/LocalBaseline'].value[i],
                                               f['LowPassSegmentation/' + k + '/AllEvents/' + '{:09d}'.format(i)].value[i],
                                               np.ones(baseline) *
                                               f['LowPassSegmentation/' + k + '/LocalBaseline'].value[i]]))
                    else:
                        fit1 = np.concatenate([np.ones(baseline) *
                                               f['LowPassSegmentation/' + k + '/LocalBaseline'].value[i],
                                               np.ones(np.uint64(endp_i1 - startp_i1 - 2 * baseline)) *
                                               f['LowPassSegmentation/' + k + '/FitLevel'].value[i],
                                               np.ones(baseline) *
                                               f['LowPassSegmentation/' + k + '/LocalBaseline'].value[i]])
                    line1f.set_data(np.arange(endp_i1-startp_i1)/out['samplerate'], fit1)
                line1.set_data(np.arange(endp_i1-startp_i1)/out['samplerate'], out['i1'][startp_i1:endp_i1])
                line2.set_data(np.arange(endp_i1-startp_i1)/out['samplerate'], out['v1'][startp_i1:endp_i1])

                ax1.set_title('Event {}\n {}'.format(i, file[1:-11]))
                ax1.relim()
                ax2.relim()
                ax1.autoscale_view(True, True, True)
                ax2.autoscale_view(True, True, True)
                fig.canvas.draw()
                pp.savefig(fig)
                print('{}: {} out of {} saved!'.format(toPlot, str(i), NumberOfEvents-1))
            if not AllInOneFile:
                pp.close()
if AllInOneFile:
    pp.close()