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
import MiscParameters as pm
from matplotlib.ticker import EngFormatter
import h5py
from tkinter.filedialog import askopenfilenames
from matplotlib import rc
rc('mathtext', default='regular')
pm.init(LoadFiles = 0)

baseline = 100
fitOn = 1
percentage = 1
plotPoints=0

filenames = {'/Users/migraf/SWITCHdrive/PhD/BPS/Poster Resources/Data/17B_1MKCl_IV_80mer_both_Ch25_8_OriginalDB.hdf5'}
#filenames = askopenfilenames()

for filename in filenames:
    print(filename)
    f = h5py.File(filename, 'r')
    out = fu.OpenFile(f['General/FileName'].value)
    file = os.sep + str(os.path.split(filename)[1][:-5])

    # Graphene Or Not??
    if not out['graphene']:
        ## Only need to plot i1 and i1_Up and v1
        toPlot = ['i1', 'i1_Up']
        for k in toPlot:
            NumberOfEvents=f['LowPassSegmentation/' + k + '/NumberOfEvents'].value
            if NumberOfEvents is not 0:
                pp = PdfPages(pm.OutputFolder + file + '_Only_' + k + '_EventPlots.pdf')
                fig = plt.figure()
                ax1 = fig.add_subplot(211)
                ax2 = fig.add_subplot(212, sharex=ax1)
                ax1.set_title('Ionic Current Trace "' + k + '"')
                ax2.set_title('Ionic Voltage')
                ax1.set_ylabel('Current [nA]')
                ax2.set_ylabel('Voltage [V]')
                ax1.set_xlabel('Time [s]')
                ax2.set_xlabel('Time [s]')
                ax1.yaxis.set_major_formatter(EngFormatter(unit='A', places=1))
                ax2.yaxis.set_major_formatter(EngFormatter(unit='V', places=1))
                ax1.xaxis.set_major_formatter(EngFormatter(unit='s', places=1))
                ax2.xaxis.set_major_formatter(EngFormatter(unit='s', places=1))
                line1, = ax1.plot([])
                line2, = ax2.plot([], 'c')
                if fitOn:
                    line1f, = ax1.plot([], 'r')
                for i in np.linspace(0, NumberOfEvents-1, NumberOfEvents*percentage, dtype=np.uint64):
                    startp_i1 = np.int64(
                        f['LowPassSegmentation/' + k + '/' + 'StartPoints'].value[i] - baseline)
                    endp_i1 = np.int64(
                        f['LowPassSegmentation/' + k + '/' + 'EndPoints'].value[i] + baseline)
                    if startp_i1 < 0:
                        startp_i1 = 0
                    if endp_i1 > len(out['i1']):
                        endp_i1 = len(out['i1'])
                    if fitOn:
                        if plotPoints:
                            fit1 = np.concatenate([np.ones(baseline) *
                                                   f['LowPassSegmentation/' + k + '/LocalBaseline'].value[i],
                                                   f['LowPassSegmentation/' + k + '/AllEvents/' + '{:09d}'.format(i)].value[i],
                                                   np.ones(baseline) *
                                                   f['LowPassSegmentation/' + k + '/LocalBaseline'].value[i]])
                        else:
                            fit1 = np.concatenate([np.ones(baseline) *
                                                   f['LowPassSegmentation/' + k + '/LocalBaseline'].value[i],
                                                   np.ones(endp_i1 - startp_i1 - 2 * baseline) *
                                                   f['LowPassSegmentation/' + k + '/FitLevel'].value[i],
                                                   np.ones(baseline) *
                                                   f['LowPassSegmentation/' + k + '/LocalBaseline'].value[i]])
                        line1f.set_data(np.arange(startp_i1, endp_i1)/out['samplerate'], fit1)
                    line1.set_data(np.arange(startp_i1, endp_i1)/out['samplerate'], out['i1'][startp_i1:endp_i1])
                    line2.set_data(np.arange(startp_i1, endp_i1)/out['samplerate'], out['v1'][startp_i1:endp_i1])

                    ax1.set_title('{} Event {}\n Ionic Current Trace'.format(filter, i))
                    ax1.relim()
                    ax2.relim()
                    ax1.autoscale_view(True, True, True)
                    ax2.autoscale_view(True, True, True)
                    fig.canvas.draw()
                    pp.savefig(fig)
                    print('{}: {} out of {} saved!'.format(toPlot, str(i), NumberOfEvents-1))
                pp.close()
    else:
        toPlot = []
        #Initialize all the cases
        toPlot.append(['i1', '', np.uint64(f['LowPassSegmentation/i1/OnlyIndexes'].value), np.uint64(f['LowPassSegmentation/i1/OnlyIndexes'].value)])
        #toPlot.append(['i1_Up', '', np.uint64(f['LowPassSegmentation/i1_Up/OnlyIndexes'].value), np.uint64(f['LowPassSegmentation/i1_Up/OnlyIndexes'].value)])
        #toPlot.append(['', 'i2', np.uint64(f['LowPassSegmentation/i2/OnlyIndexes'].value), np.uint64(f['LowPassSegmentation/i2/OnlyIndexes'].value)])
        #toPlot.append(['', 'i2_Up', np.uint64(f['LowPassSegmentation/i2_Up/OnlyIndexes'].value), np.uint64(f['LowPassSegmentation/i2_Up/OnlyIndexes'].value)])
        #toPlot.append(['i1', 'i2', np.uint64(f['LowPassSegmentation/i1/CommonIndexesWithi2'].value), np.uint64(f['LowPassSegmentation/i2/CommonIndexesWithi1'].value)])
        #toPlot.append(['i1_Up', 'i2', np.uint64(f['LowPassSegmentation/i1_Up/CommonIndexesWithi2'].value), np.uint64(f['LowPassSegmentation/i2/CommonIndexesWithi1_Up'].value)])
        #toPlot.append(['i1_Up', 'i2_Up', np.uint64(f['LowPassSegmentation/i1_Up/CommonIndexesWithi2_Up'].value), np.uint64(f['LowPassSegmentation/i2_Up/CommonIndexesWithi1_Up'].value)])
        #toPlot.append(['i1', 'i2_Up', np.uint64(f['LowPassSegmentation/i1/CommonIndexesWithi2_Up'].value), np.uint64(f['LowPassSegmentation/i2_Up/CommonIndexesWithi1'].value)])


        line1data = np.array([])
        line2data = np.array([])
        line3data = np.array([])
        line4data = np.array([])
        line1fdata = np.array([])
        line2fdata = np.array([])
        for cases in toPlot:
            NumberOfEvents = len(cases[2])
            for i in np.linspace(0, NumberOfEvents-1, NumberOfEvents*percentage, dtype=np.uint64):
                if not cases[0] == '':
                    startp_i1 = np.int64(
                        f['LowPassSegmentation/' + cases[0] + '/' + 'StartPoints'].value[cases[2][i]] - baseline)
                    endp_i1 = np.int64(
                        f['LowPassSegmentation/' + cases[0] + '/' + 'EndPoints'].value[cases[2][i]] + baseline)
                else:
                    startp_i1 = np.int64(
                        f['LowPassSegmentation/' + cases[1] + '/' + 'StartPoints'].value[cases[3][i]] - baseline)
                    endp_i1 = np.int64(
                        f['LowPassSegmentation/' + cases[1] + '/' + 'EndPoints'].value[cases[3][i]] + baseline)
                if not cases[1] == '':
                    startp_i2 = np.int64(
                        f['LowPassSegmentation/' + cases[1] + '/' + 'StartPoints'].value[cases[3][i]] - baseline)
                    endp_i2 = np.int64(
                        f['LowPassSegmentation/' + cases[1] + '/' + 'EndPoints'].value[cases[3][i]] + baseline)
                else:
                    startp_i2 = np.int64(
                        f['LowPassSegmentation/' + cases[0] + '/' + 'StartPoints'].value[cases[2][i]] - baseline)
                    endp_i2 = np.int64(
                        f['LowPassSegmentation/' + cases[0] + '/' + 'EndPoints'].value[cases[2][i]] + baseline)

                if startp_i1 < 0:
                    startp_i1 = 0
                if startp_i2 < 0:
                    startp_i2 = 0
                if endp_i1 > len(out['i1']):
                    endp_i1 = len(out['i1'])
                if endp_i2 > len(out['i2']):
                    endp_i2 = len(out['i2'])

                i1 = out['i1'][startp_i1:endp_i1]
                v1 = out['v1'][startp_i1:endp_i1]
                i2 = out['i2'][startp_i2:endp_i2]
                v2 = out['v2'][startp_i2:endp_i2]

                diff = len(i1)-len(i2)
                if diff>0: #i1 is longer
                    i1 = i1[0:len(i2)]
                    v1 = v1[0:len(i2)]
                if diff < 0:  # i2 is longer
                    i2 = i2[0:len(i1)]
                    v2 = v2[0:len(i1)]

                line1data = np.append(line1data, i1)
                line2data = np.append(line2data, i2)
                line3data = np.append(line3data, v1)
                line4data = np.append(line4data, v2)


                if fitOn:
                    if not cases[0] == '':
                        if plotPoints:
                            fit1 = np.concatenate([np.ones(baseline) *
                                                   f['LowPassSegmentation/' + cases[0] + '/LocalBaseline'].value[cases[2][i]],
                                                   f['LowPassSegmentation/' + cases[0] + '/AllEvents/' + '{:09d}'.format(cases[2][i])].value,
                                                   np.ones(baseline) *
                                                   f['LowPassSegmentation/' + cases[0] + '/LocalBaseline'].value[cases[2][i]]])
                        else:
                            fit1 = np.concatenate([np.ones(baseline) *
                                                   f['LowPassSegmentation/' + cases[0] + '/LocalBaseline'].value[cases[2][i]],
                                                   np.ones(endp_i1 - startp_i1 - 2 * baseline) * f[
                                                       'LowPassSegmentation/' + cases[0] + '/FitLevel'].value[cases[2][i]],
                                                   np.ones(baseline) *
                                                   f['LowPassSegmentation/' + cases[0] + '/LocalBaseline'].value[cases[2][i]]])
                            if diff > 0:  # i1 is longer
                                fit1 = fit1[0:len(i2)]
                            line1fdata = np.append(line1fdata, fit1)
                    if not cases[1] == '':
                        if plotPoints:
                            fit2 = np.concatenate([np.ones(baseline) *
                                                   f['LowPassSegmentation/' + cases[1] + '/LocalBaseline'].value[cases[3][i]],
                                                   f['LowPassSegmentation/' + cases[1] + '/AllEvents/' + '{:09d}'.format(
                                                       cases[3][i])].value,
                                                   np.ones(baseline) *
                                                   f['LowPassSegmentation/' + cases[1] + '/LocalBaseline'].value[cases[3][i]]])
                        else:
                            fit2 = np.concatenate([np.ones(baseline) *
                                                   f['LowPassSegmentation/' + cases[1] + '/LocalBaseline'].value[cases[3][i]],
                                                   np.ones(endp_i2 - startp_i2 - 2 * baseline) * f['LowPassSegmentation/' + cases[1] + '/FitLevel'].value[cases[3][i]],
                                                   np.ones(baseline) *
                                                   f['LowPassSegmentation/' + cases[1] + '/LocalBaseline'].value[cases[3][i]]])
                            if diff < 0:  # i2 is longer
                                fit2 = fit2[0:len(i1)]
                            line2fdata = np.append(line2fdata, fit2)


            fig = plt.figure()
            ax1 = fig.add_subplot(311)
            ax2 = fig.add_subplot(312, sharex=ax1)
            ax3 = fig.add_subplot(313, sharex=ax1)
            ax4 = ax3.twinx()
            ax4.set_title('Voltage Trace')
            ax1.set_ylabel('Current [nA]')
            ax2.set_ylabel('Current [nA]')
            ax3.set_ylabel('Ionic Voltage [V]')
            ax3.set_xlabel('Time [s]')
            ax1.set_xlabel('Time [s]')
            ax2.set_xlabel('Time [s]')
            ax1.yaxis.set_major_formatter(EngFormatter(unit='A', places=1))
            ax2.yaxis.set_major_formatter(EngFormatter(unit='A', places=1))
            ax3.yaxis.set_major_formatter(EngFormatter(unit='V', places=1))
            ax4.yaxis.set_major_formatter(EngFormatter(unit='V', places=1))
            ax3.xaxis.set_major_formatter(EngFormatter(unit='s', places=1))
            ax1.xaxis.set_major_formatter(EngFormatter(unit='s', places=1))
            ax2.xaxis.set_major_formatter(EngFormatter(unit='s', places=1))
            ax4.set_ylabel('Tr. Voltage [V]')

            t = np.arange(len(line1data))/out['samplerate']
            t2 = np.arange(len(line2data))/out['samplerate']

            ax1.plot(t, line1data)
            ax2.plot(t2, line2data)
            ax3.plot(t, line3data, 'c')
            ax4.plot(t2, line4data, 'y')
            if fitOn:
                ax1.plot(t, line1fdata, 'r')
                #ax2.plot(t2, line2fdata, 'r')
            ax4.tick_params(axis='y', colors='y')
            ax4.yaxis.label.set_color('y')
            ax3.tick_params(axis='y', colors='c')
            ax3.yaxis.label.set_color('c')
            ax1.set_title('Event {} of {}\n Ionic Current Trace "{}"'.format(i, NumberOfEvents, cases[0]))
            ax1.relim()
            ax2.relim()
            ax3.relim()
            ax4.relim()
            ax1.autoscale_view(True, True, True)
            ax2.autoscale_view(True, True, True)
            # voltage plots
            ax3.autoscale_view(True, True, False)
            ax4.autoscale_view(True, True, False)
            ax3.set_ylim(np.min(line3data) - np.mean(line3data) * 5 / 10, np.max(line3data) + np.mean(line3data) * 1 / 10)
            ax4.set_ylim([np.min(line4data) - np.mean(line4data) * 2 / 10, np.max(line4data) + np.mean(line4data) * 5 / 10])
            fig.tight_layout()
            fig.canvas.draw()
            plt.show()


