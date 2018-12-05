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
import h5py
from tkinter.filedialog import askopenfilenames
from matplotlib import rc
rc('mathtext', default='regular')
pm.init(LoadFiles = 0)

baseline = 250
fitOn = 1

filenames = {'/Users/migraf/Desktop/Chan/Data/R43_KCl gradient_100mMcis_1mM trans_4_20180405_182358_OriginalDB.hdf5'}
#filenames = askopenfilenames()

for filename in filenames:
    print(filename)
    f = h5py.File(filename, 'r')
    out = fu.OpenFile(f['General/FileName'].value)
    file = os.sep + str(os.path.split(filename)[1][:-5])

    Combinations={}


    categories1 = ['i1', 'i1_Up']



    if out['graphene']:
        filters = ['CommonIndexes', 'OnlyIndexesi1', 'OnlyIndexesi2']
        categories2 = ['i2', 'i2_Up']
    else:
        categories2 = ['']
        filters = ['OnlyIndexesi1']

    for cat1 in categories1:
        for cat2 in categories2:
            for filter in filters:
                print(filter)
                pp = PdfPages(pm.OutputFolder + file + '_' + cat1 + '_vs_' + cat2 + '_' + filter + '_EventPlots.pdf')
                fig = plt.figure()
                ax1 = fig.add_subplot(311)
                ax2 = fig.add_subplot(312, sharex = ax1)
                ax3 = fig.add_subplot(313, sharex = ax1)
                ax4 = ax3.twinx()
                ax1.set_title('Ionic Current Trace "' + cat1 + '"')
                ax2.set_title('Transverse Current Trace "' + cat2 + '"')
                ax4.set_title('Voltage Trace')
                ax1.set_ylabel('Current [nA]')
                ax2.set_ylabel('Current [nA]')
                ax3.set_ylabel('Ionic Voltage [V]')
                ax3.set_xlabel('Time [ms]')
                ax1.set_xlabel('Time [ms]')
                ax2.set_xlabel('Time [ms]')
                ax4.set_ylabel('Tr. Voltage [V]')
                line1, = ax1.plot([])
                line2, = ax2.plot([])
                line3, = ax3.plot([], 'c')
                line4, = ax4.plot([], 'y')
                ax4.tick_params(axis='y', colors = 'y')
                ax4.yaxis.label.set_color('y')
                ax3.tick_params(axis='y', colors = 'c')
                ax3.yaxis.label.set_color('c')

                if fitOn:
                    line1f, = ax1.plot([], 'r')
                    line2f, = ax2.plot([], 'r')

                # Generate Event List
                if filter is 'CommonIndexes':
                    i1_ActualEventList = np.uint64(f['LowPassSegmentation/' + cat1 + '/CommonIndexesWith' + cat2].value)
                    i2_ActualEventList = np.uint64(f['LowPassSegmentation/' + cat2 + '/CommonIndexesWith' + cat1].value)
                    NumberOfEvents = len(i2_ActualEventList)
                if filter is 'OnlyIndexesi1' and out['graphene']:
                    i1_ActualEventList = np.uint64(f['LowPassSegmentation/' + cat1 + '/OnlyIndexesWith' + cat2].value)
                    i2_ActualEventList = np.uint64(f['LowPassSegmentation/' + cat1 + '/OnlyIndexesWith' + cat2].value)
                    NumberOfEvents = len(i1_ActualEventList)
                if filter is 'OnlyIndexesi1' and not out['graphene']:
                    i1_ActualEventList = np.arange(f['LowPassSegmentation/' + cat1 + '/NumberOfEvents'].value)
                    NumberOfEvents = len(i1_ActualEventList)
                if filter is 'OnlyIndexesi2':
                    i1_ActualEventList = np.uint64(f['LowPassSegmentation/' + cat2 + '/OnlyIndexesWith' + cat1].value)
                    i2_ActualEventList = np.uint64(f['LowPassSegmentation/' + cat2 + '/OnlyIndexesWith' + cat1].value)
                    NumberOfEvents = len(i1_ActualEventList)

                for i in np.arange(0, NumberOfEvents):
                    if filter is 'OnlyIndexesi2':
                        startp_i2 = np.int64(
                            f['LowPassSegmentation/'+ cat2 + '/' + 'StartPoints'].value[i2_ActualEventList[i]] - baseline)
                        endp_i2 = np.int64(f['LowPassSegmentation/' + cat2 + '/' + 'EndPoints'].value[i2_ActualEventList[i]] + baseline)
                        startp_i1 = startp_i2
                        endp_i1 = endp_i2
                    if filter is 'OnlyIndexesi1':
                        startp_i1 = np.int64(
                            f['LowPassSegmentation/' + cat1 + '/' + 'StartPoints'].value[i1_ActualEventList[i]] - baseline)
                        endp_i1 = np.int64(f['LowPassSegmentation/' + cat1 + '/' + 'EndPoints'].value[i1_ActualEventList[i]] + baseline)
                        startp_i2 = startp_i1
                        endp_i2 = endp_i1
                    if filter is 'CommonIndexes':
                        startp_i1 = np.int64(f['LowPassSegmentation/' + cat1 + '/' + 'StartPoints'].value[i1_ActualEventList[i]] - baseline)
                        endp_i1 = np.int64(f['LowPassSegmentation/' + cat1 + '/' + 'EndPoints'].value[i1_ActualEventList[i]] + baseline)
                        startp_i2 = np.int64(f['LowPassSegmentation/' + cat2 + '/' + 'StartPoints'].value[i2_ActualEventList[i]] - baseline)
                        endp_i2 = np.int64(f['LowPassSegmentation/' + cat2 + '/' + 'EndPoints'].value[i2_ActualEventList[i]] + baseline)
                    if startp_i1 < 0:
                        startp_i1 = 0
                    if startp_i2 < 0:
                        startp_i2 = 0
                    if endp_i1 > len(out['i1']):
                        endp_i1 = len(out['i1'])
                    if out['graphene']:
                        if endp_i2 > len(out['i2']):
                            endp_i2 = len(out['i2'])
                        i2 = out['i2'][startp_i2:endp_i2]
                        v2 = out['v2'][startp_i1:endp_i1]
                        t2 = np.arange(startp_i2, endp_i2) / out['samplerate']

                    i1 = out['i1'][startp_i1:endp_i1]
                    v1 = out['v1'][startp_i1:endp_i1]
                    t1 = np.arange(startp_i1, endp_i1)/out['samplerate']

                    if fitOn:
                        if filter is 'OnlyIndexesi1' and (endp_i2 - startp_i2 - 2 * baseline)>0:
                            if len(cat1) > 2:
                                fit1 = np.concatenate([np.ones(baseline) * f['LowPassSegmentation/' + cat1 + '/LocalBaseline'].value[i1_ActualEventList[i]],
                                                   np.ones(endp_i1-startp_i1-2*baseline) * (
                                                       f['LowPassSegmentation/' + cat1 + '/LocalBaseline'].value[i1_ActualEventList[i]] +
                                                       f['LowPassSegmentation/' + cat1 + '/DeltaI'].value[i1_ActualEventList[i]]),
                                                   np.ones(baseline) * f['LowPassSegmentation/' + cat1 + '/LocalBaseline'].value[i1_ActualEventList[i]]])
                            else:
                                fit1 = np.concatenate([np.ones(baseline) *
                                                       f['LowPassSegmentation/' + cat1 + '/LocalBaseline'].value[
                                                           i1_ActualEventList[i]],
                                                       np.ones(endp_i1 - startp_i1 - 2 * baseline) * (
                                                           f['LowPassSegmentation/' + cat1 + '/LocalBaseline'].value[
                                                               i1_ActualEventList[i]] -
                                                           f['LowPassSegmentation/' + cat1 + '/DeltaI'].value[
                                                               i1_ActualEventList[i]]),
                                                       np.ones(baseline) *
                                                       f['LowPassSegmentation/' + cat1 + '/LocalBaseline'].value[
                                                           i1_ActualEventList[i]]])

                            fit2=0
                            line1f.set_data(t1, fit1)
                            line2f.set_data([], [])
                        elif filter is 'OnlyIndexesi2' and (endp_i2 - startp_i2 - 2 * baseline)>0:
                            if len(cat1) > 2:
                                fit2 = np.concatenate([np.ones(baseline) * f['LowPassSegmentation/' + cat2 + '/LocalBaseline'].value[i2_ActualEventList[i]],
                                               np.ones(endp_i2-startp_i2-2*baseline) * (
                                                   f['LowPassSegmentation/' + cat2 + '/LocalBaseline'].value[i2_ActualEventList[i]] +
                                                   f['LowPassSegmentation/' + cat2 + '/DeltaI'].value[i2_ActualEventList[i]]),
                                               np.ones(baseline) * f['LowPassSegmentation/' + cat2 + '/LocalBaseline'].value[i2_ActualEventList[i]]])
                            else:
                                fit2 = np.concatenate([np.ones(baseline) *
                                                       f['LowPassSegmentation/' + cat2 + '/LocalBaseline'].value[
                                                           i2_ActualEventList[i]],
                                                       np.ones(endp_i2 - startp_i2 - 2 * baseline) * (
                                                           f['LowPassSegmentation/' + cat2 + '/LocalBaseline'].value[
                                                               i2_ActualEventList[i]] -
                                                           f['LowPassSegmentation/' + cat2 + '/DeltaI'].value[
                                                               i2_ActualEventList[i]]),
                                                       np.ones(baseline) *
                                                       f['LowPassSegmentation/' + cat2 + '/LocalBaseline'].value[
                                                           i2_ActualEventList[i]]])
                            fit1=0
                            line1f.set_data([], [])
                            line2f.set_data(t2, fit2)
                        elif filter is 'CommonIndexes' and (endp_i2 - startp_i2 - 2 * baseline)>0:
                            if len(cat2) > 2:
                                fit2 = np.concatenate(
                                    [np.ones(baseline) * f['LowPassSegmentation/' + cat2 + '/LocalBaseline'].value[i2_ActualEventList[i]],
                                     np.ones(endp_i2 - startp_i2 - 2 * baseline) * (
                                         f['LowPassSegmentation/' + cat2 + '/LocalBaseline'].value[i2_ActualEventList[i]] +
                                         f['LowPassSegmentation/' + cat2 + '/DeltaI'].value[i2_ActualEventList[i]]),
                                     np.ones(baseline) * f['LowPassSegmentation/' + cat2 + '/LocalBaseline'].value[i2_ActualEventList[i]]])
                            else:
                                fit2 = np.concatenate(
                                    [np.ones(baseline) * f['LowPassSegmentation/' + cat2 + '/LocalBaseline'].value[
                                        i2_ActualEventList[i]],
                                     np.ones(endp_i2 - startp_i2 - 2 * baseline) * (
                                         f['LowPassSegmentation/' + cat2 + '/LocalBaseline'].value[i2_ActualEventList[i]] -
                                         f['LowPassSegmentation/' + cat2 + '/DeltaI'].value[i2_ActualEventList[i]]),
                                     np.ones(baseline) * f['LowPassSegmentation/' + cat2 + '/LocalBaseline'].value[
                                         i2_ActualEventList[i]]])
                            if len(cat1) > 2:
                                fit1 = np.concatenate(
                                    [np.ones(baseline) * f['LowPassSegmentation/' + cat1 + '/LocalBaseline'].value[i1_ActualEventList[i]],
                                     np.ones(endp_i1 - startp_i1 - 2 * baseline) * (
                                         f['LowPassSegmentation/' + cat1 + '/LocalBaseline'].value[i1_ActualEventList[i]] +
                                         f['LowPassSegmentation/' + cat1 + '/DeltaI'].value[i1_ActualEventList[i]]),
                                     np.ones(baseline) * f['LowPassSegmentation/' + cat1 + '/LocalBaseline'].value[i1_ActualEventList[i]]])
                            else:
                                fit1 = np.concatenate(
                                    [np.ones(baseline) * f['LowPassSegmentation/' + cat1 + '/LocalBaseline'].value[i1_ActualEventList[i]],
                                     np.ones(endp_i1 - startp_i1 - 2 * baseline) * (
                                         f['LowPassSegmentation/' + cat1 + '/LocalBaseline'].value[i1_ActualEventList[i]] -
                                         f['LowPassSegmentation/' + cat1 + '/DeltaI'].value[i1_ActualEventList[i]]),
                                     np.ones(baseline) * f['LowPassSegmentation/' + cat1 + '/LocalBaseline'].value[i1_ActualEventList[i]]])

                            line1f.set_data(t1, fit1)
                            line2f.set_data(t2, fit2)
                        else:
                            line1f.set_data([], [])
                            line2f.set_data([], [])


                    #print('Values\t i1: {}, t1:{}, i2:{}, t2:{}\nFit1:{}, Fit2:{}'.format(len(i1),len(t1), len(i2),len(t2), len(fit1),len(fit2)))

                    line1.set_data(t1, i1)
                    line3.set_data(t1, v1)
                    if out['graphene']:
                        line4.set_data(t1, v2)
                        line2.set_data(t2, i2)

                    ax1.set_title('{} Event {}\n Ionic Current Trace'.format(filter, i))
                    ax1.relim()
                    ax2.relim()
                    ax3.relim()
                    ax4.relim()
                    ax1.autoscale_view(True, True, True)
                    ax2.autoscale_view(True, True, True)
                    #voltage plots
                    ax3.autoscale_view(True, True, False)
                    ax4.autoscale_view(True, True, False)
                    ax3.set_ylim(np.min(v1)-np.mean(v1)*5/10, np.max(v1)+np.mean(v1)*1/10)
                    if out['graphene']:
                        ax4.set_ylim([np.min(v2)-np.mean(v2)*2/10, np.max(v2)+np.mean(v2)*5/10])
                    fig.tight_layout()
                    fig.canvas.draw()
                    pp.savefig(fig)
                    print('{} out of {} saved!'.format(str(i), NumberOfEvents-1))

                pp.close()
                if NumberOfEvents == 0:
                    print('Empty PDF deleted: {}'.format(pm.OutputFolder + file + '_' + cat1 + '_vs_' + cat2 + '_' + filter + '_EventPlots.pdf'))
                    os.remove(pm.OutputFolder + file + '_' + cat1 + '_vs_' + cat2 + '_' + filter + '_EventPlots.pdf')
                plt.close()