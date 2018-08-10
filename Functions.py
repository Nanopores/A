import matplotlib
#matplotlib.use('Qt5Agg')
import numpy as np
import math
import scipy
import scipy.signal as sig
import os
import pickle as pkl
from scipy import constants as cst
from scipy import io
from scipy import signal
from PyQt5 import QtGui, QtWidgets
from numpy import linalg as lin
import pyqtgraph as pg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import h5py
from timeit import default_timer as timer
import platform
from scipy.optimize import curve_fit
from sys import platform
#import pyabf
from matplotlib.ticker import EngFormatter

if not 'verbose' in globals():
    verbose = False

verboseprint = print if verbose else lambda *a, **k: None


def GetKClConductivity(Conc, Temp):
    p = pkl.load(open('KCl_ConductivityValues.p', 'rb'))
    if Conc==1.0:
        Conc = np.uint(1)
    return np.polyval(p[str(Conc)], Temp)

def GetTempFromKClConductivity(Conc, Cond):
    p = pkl.load(open('KCl_ConductivityValues.p', 'rb'))
    if Conc==1.0:
        Conc = np.uint(1)
    return (Cond-p[str(Conc)][1])/p[str(Conc)][0]



def RecursiveLowPassFast(signal, coeff, samplerate):
    padlen = np.uint64(samplerate)
    prepadded = np.ones(padlen)*np.mean(signal[0:1000])
    signaltofilter = np.concatenate((prepadded, signal))

    mltemp = scipy.signal.lfilter([1 - coeff['a'], 0], [1, -coeff['a']], signaltofilter)
    vltemp = scipy.signal.lfilter([1 - coeff['a'], 0], [1, -coeff['a']], np.square(signaltofilter - mltemp))

    ml = np.delete(mltemp, np.arange(padlen))
    vl = np.delete(vltemp, np.arange(padlen))

    sl = ml - coeff['S'] * np.sqrt(vl)
    Ni = len(signal)
    points = np.array(np.where(signal <= sl)[0])
    to_pop = np.array([])
    for i in range(1, len(points)):
        if points[i] - points[i - 1] == 1:
            to_pop = np.append(to_pop, i)
    points = np.unique(np.delete(points, to_pop))
    RoughEventLocations = []
    NumberOfEvents = 0

    for i in points:
        if NumberOfEvents is not 0:
            if i >= RoughEventLocations[NumberOfEvents-1][0] and i <= RoughEventLocations[NumberOfEvents-1][1]:
                continue
        NumberOfEvents += 1
        start = i
        El = ml[i] + coeff['E'] * np.sqrt(vl[i])
        Mm = ml[i]
        Vv = vl[i]
        duration = 0
        endp = start
        if (endp + 1) < len(signal):
            while signal[endp + 1] < El and endp < (Ni - 2):# and duration < coeff['eventlengthLimit']*samplerate:
                duration += 1
                endp += 1
        if duration >= coeff['eventlengthLimit'] * samplerate or endp > (Ni - 10): # or duration <= coeff['minEventLength'] * samplerate:
            NumberOfEvents -= 1
            continue
        else:
            k = start
            while signal[k] < Mm and k > 1:
                k -= 1
            start = k - 1
            k2 = i + 1
            #while signal[k2] > Mm:
            #    k2 -= 1
            #endp = k2
            if start < 0:
                start = 0
            RoughEventLocations.append((start, endp, ml[start], vl[start]))

    return np.array(RoughEventLocations)#, ml, vl, sl

def RecursiveLowPassFastUp(signal, coeff, samplerate):
    ml = scipy.signal.lfilter([1 - coeff['a'], 0], [1, -coeff['a']], signal)
    vl = scipy.signal.lfilter([1 - coeff['a'], 0], [1, -coeff['a']], np.square(signal - ml))
    sl = ml + coeff['S'] * np.sqrt(vl)
    Ni = len(signal)
    points = np.array(np.where(signal>=sl)[0])
    to_pop=np.array([])
    for i in range(1,len(points)):
        if points[i] - points[i - 1] == 1:
            to_pop=np.append(to_pop, i)
    points = np.delete(points, to_pop)

    points =np.delete(points, np.array(np.where(points == 0)[0]))

    RoughEventLocations = []
    NumberOfEvents=0
    for i in points:
        if NumberOfEvents is not 0:
            if i >= RoughEventLocations[NumberOfEvents-1][0] and i <= RoughEventLocations[NumberOfEvents-1][1]:
                continue
        NumberOfEvents += 1
        start = i
        El = ml[i] + coeff['E'] * np.sqrt(vl[i])
        Mm = ml[i]
        duration = 0
        while signal[i + 1] > El and i < (Ni - 2) and duration < coeff['eventlengthLimit']*samplerate:
            duration += 1
            i += 1
        if duration >= coeff['eventlengthLimit']*samplerate or i > (Ni - 10):
            NumberOfEvents -= 1
        else:
            k = start
            while signal[k] > Mm and k > 2:
                k -= 1
            start = k - 1
            k2 = i + 1
            while signal[k2] > Mm:
                k2 -= 1
            endp = k2
            RoughEventLocations.append((start, endp, ml[start], vl[start]))

    return np.array(RoughEventLocations)

def ImportABF(datafilename):
    abf = pyabf.ABF(datafilename)
    abf.info()  # shows what is available
    output={'type': 'Clampfit', 'graphene': 0, 'samplerate': abf.pointsPerSec, 'i1': -20000./65536 * abf.dataY, 'v1': abf.dataC, 'filename': datafilename}
    return output

def ImportAxopatchData(datafilename):
    x=np.fromfile(datafilename, np.dtype('>f4'))
    f=open(datafilename, 'rb')
    graphene=0
    for i in range(0, 10):
        a=str(f.readline())
        #print(a)
        if 'Acquisition' in a or 'Sample Rate' in a:
            samplerate=int(''.join(i for i in a if i.isdigit()))/1000
        if 'FEMTO preamp Bandwidth' in a:
            femtoLP=int(''.join(i for i in a if i.isdigit()))
        if 'I_Graphene' in a:
            graphene=1
            print('This File Has a Graphene Channel!')
    end = len(x)
    if graphene:
        #pore current
        i1 = x[250:end-3:4]
        #graphene current
        i2 = x[251:end-2:4]
        #pore voltage
        v1 = x[252:end-1:4]
        #graphene voltage
        v2 = x[253:end:4]
        print('The femto was set to : {} Hz, if this value was correctly entered in the LabView!'.format(str(femtoLP)))
        output={'FemtoLowPass': femtoLP, 'type': 'Axopatch', 'graphene': 1, 'samplerate': samplerate, 'i1': i1, 'v1': v1, 'i2': i2, 'v2': v2, 'filename': datafilename}
    else:
        i1 = np.array(x[250:end-1:2])
        v1 = np.array(x[251:end:2])
        output={'type': 'Axopatch', 'graphene': 0, 'samplerate': samplerate, 'i1': i1, 'v1': v1, 'filename': datafilename}
    return output

def ImportChimeraRaw(datafilename):
    matfile=io.loadmat(str(os.path.splitext(datafilename)[0]))
    #buffersize=matfile['DisplayBuffer']
    data = np.fromfile(datafilename, np.dtype('<u2'))
    samplerate = np.float64(matfile['ADCSAMPLERATE'])
    TIAgain = np.int32(matfile['SETUP_TIAgain'])
    preADCgain = np.float64(matfile['SETUP_preADCgain'])
    currentoffset = np.float64(matfile['SETUP_pAoffset'])
    ADCvref = np.float64(matfile['SETUP_ADCVREF'])
    ADCbits = np.int32(matfile['SETUP_ADCBITS'])

    closedloop_gain = TIAgain * preADCgain
    bitmask = (2 ** 16 - 1) - (2 ** (16 - ADCbits) - 1)
    data = -ADCvref + (2 * ADCvref) * (data & bitmask) / 2 ** 16
    data = (data / closedloop_gain + currentoffset)
    data.shape = [data.shape[1], ]
    output = {'matfilename': str(os.path.splitext(datafilename)[0]),'i1raw': data, 'v1': np.float64(matfile['SETUP_biasvoltage']), 'samplerateRaw': np.int64(samplerate), 'type': 'ChimeraRaw', 'filename': datafilename, 'graphene': 0}
    return output

def ImportChimeraData(datafilename):
    matfile = io.loadmat(str(os.path.splitext(datafilename)[0]))
    samplerate = matfile['ADCSAMPLERATE']
    if samplerate<4e6:
        data = np.fromfile(datafilename, np.dtype('float64'))
        buffersize = matfile['DisplayBuffer']
        out = Reshape1DTo2D(data, buffersize)
        output = {'i1': out['i1'], 'v1': out['v1'], 'samplerate':float(samplerate), 'type': 'ChimeraNotRaw', 'filename': datafilename}
    else:
        output = ImportChimeraRaw(datafilename)
    return output

def OpenFile(filename = '', ChimeraLowPass = 10e3):
    if filename == '':
        datafilename = QtGui.QFileDialog.getOpenFileName()
        datafilename=datafilename[0]
        print(datafilename)
    else:
        datafilename=filename

    verboseprint('Loading file... ' +filename)

    if datafilename[-3::] == 'dat':
        isdat = 1
        output = ImportAxopatchData(datafilename)
    elif datafilename[-3::] == 'log':
        isdat = 0
        output = ImportChimeraData(datafilename)
        verboseprint('length: ' + str(len(output['i1raw'])))
        if output['type'] is 'ChimeraRaw':  # Lowpass and downsample
            Wn = round(2 * ChimeraLowPass / output['samplerateRaw'], 4)  # [0,1] nyquist frequency
            b, a = signal.bessel(4, Wn, btype='low', analog=False)  # 4-th order digital filter
            Filt_sig = signal.filtfilt(b, a, output['i1raw'], method = 'gust')
            ds_factor = output['samplerateRaw'] / (5 * ChimeraLowPass)
            output['i1'] = scipy.signal.resample(Filt_sig, int(len(output['i1raw']) / ds_factor))
            output['samplerate'] = output['samplerateRaw'] / ds_factor
            output['v1'] = output['v1']*np.ones(len(output['i1']))
            verboseprint('Samplerate after filtering:' + str(output['samplerate']))
            verboseprint(len(output['i1']))
    elif datafilename[-3::] == 'abf':
        output = ImportABF(datafilename)
    st = os.stat(datafilename)
    if platform == 'darwin':
        verboseprint('Platform Is Darwin')
        output['TimeFileWritten'] = st.st_birthtime
        output['TimeFileLastModified'] = st.st_mtime
        output['ExperimentDuration'] = st.st_mtime - st.st_birthtime
    else:
        verboseprint('Platform Is WinShit')
        output['TimeFileWritten'] = st.st_ctime
        output['TimeFileLastModified'] = st.st_mtime
        output['ExperimentDuration'] = st.st_mtime - st.st_ctime
    return output

def MakeIVData(output, approach = 'mean', delay = 0.7, UseVoltageRange=0):
    (ChangePoints, sweepedChannel) = CutDataIntoVoltageSegments(output)
    if ChangePoints is 0:
        return 0

    if output['graphene']:
        currents = ['i1', 'i2']
    else:
        currents = ['i1']

    Values = output[sweepedChannel][ChangePoints]
    Values = np.append(Values, output[sweepedChannel][::-1][0])
    delayinpoints = np.int64(delay * output['samplerate'])
    #   Store All Data
    All={}
    for current in currents:
        Item = {}
        ItemExp = {}
        l = len(Values)
        Item['Voltage'] = np.zeros(l)
        Item['StartPoint'] = np.zeros(l,dtype=np.uint64)
        Item['EndPoint'] = np.zeros(l, dtype=np.uint64)
        Item['Mean'] = np.zeros(l)
        Item['STD'] = np.zeros(l)
        Item['SweepedChannel'] = sweepedChannel
        ItemExp['SweepedChannel'] = sweepedChannel
        Item['YorkFitValues'] = {}
        ItemExp['YorkFitValues'] = {}
        ### MEAN METHOD###
        # First item
        Item['Voltage'][0] = Values[0]
        trace = output[current][0 + delayinpoints:ChangePoints[0]]
        Item['StartPoint'][0] = 0 + delayinpoints
        Item['EndPoint'][0] = ChangePoints[0]
        Item['Mean'][0] = np.mean(trace)
        isProblematic = math.isclose(Values[0], 0, rel_tol=1e-3) or math.isclose(Values[0], np.min(Values), rel_tol=1e-3) or math.isclose(Values[0], np.max(Values), rel_tol=1e-3)
        if sweepedChannel is 'v2' and isProblematic:
            print('In Problematic If for Values: {}, index {}'.format(Values[0],0))
            Item['STD'][0] = np.std(trace[-np.int64(output['samplerate']):])
        else:
            Item['STD'][0] = np.std(trace)
        for i in range(1, len(Values) - 1):
            trace = output[current][ChangePoints[i - 1] + delayinpoints : ChangePoints[i]]
            Item['Voltage'][i] = Values[i]
            Item['StartPoint'][i] = ChangePoints[i - 1]+ delayinpoints
            Item['EndPoint'][i] = ChangePoints[i]
            Item['Mean'][i] = np.mean(trace)
            isProblematic = math.isclose(Values[i], 0, rel_tol=1e-3) or math.isclose(Values[i], np.min(Values),
                                                                                     rel_tol=1e-3) or math.isclose(Values[i], np.max(Values), rel_tol=1e-3)
            if sweepedChannel is 'v2' and isProblematic:
                print('In Problematic If for Values: {}, index {}'.format(Values[i], i))
                Item['STD'][i] = np.std(trace[-np.int64(output['samplerate']):])
            else:
                Item['STD'][i] = np.std(trace)
            print('{}, {},{}'.format(i, ChangePoints[i - 1] + delayinpoints, ChangePoints[i]))
        # Last
        if 1:
            trace = output[current][ChangePoints[len(ChangePoints) - 1] + delayinpoints : len(output[current]) - 1]
            Item['Voltage'][-1:] = Values[len(Values) - 1]
            Item['StartPoint'][-1:] = ChangePoints[len(ChangePoints) - 1]+delayinpoints
            Item['EndPoint'][-1:] = len(output[current]) - 1
            Item['Mean'][-1:] = np.mean(trace)
            isProblematic = math.isclose(Values[-1:], 0, rel_tol=1e-3) or math.isclose(Values[-1:], np.min(Values),
                                                                                     rel_tol=1e-3) or math.isclose(
                Values[-1:], np.max(Values), rel_tol=1e-3)
            if sweepedChannel is 'v2' and isProblematic:
                print('In Problematic If for Values: {}, index {}'.format(Values[-1:], i))
                Item['STD'][-1:] = np.std(trace[-np.int64(output['samplerate']):])
            else:
                Item['STD'][-1:] = np.std(trace)

        ## GET RECTIFICATION FACTOR
        parts = {'pos': Item['Voltage'] > 0, 'neg': Item['Voltage'] < 0}
        a = {}
        b = {}
        sigma_a = {}
        sigma_b = {}
        b_save = {}
        x_values = {}

        for part in parts:
            (a[part], b[part], sigma_a[part], sigma_b[part], b_save[part]) = YorkFit(
                Item['Voltage'][parts[part]].flatten(), Item['Mean'][parts[part]].flatten(),
                1e-12 * np.ones(len(Item['Voltage'][parts[part]].flatten())),
                Item['STD'][parts[part]].flatten())
        factor = b['neg'] / b['pos']
        Item['Rectification'] = factor
        ItemExp['Rectification'] = factor

        ###EXPONENTIAL METHOD###
        if approach is not 'mean':
            l = 1000
            ItemExp['StartPoint'] = np.zeros(l, dtype=np.uint64)
            ItemExp['EndPoint'] = np.zeros(l, dtype=np.uint64)
            baselinestd = np.std(output[current][0:np.int64(1 * output['samplerate'])])
            movingstd = pd.rolling_std(output[current], 10)

            # Extract The Desired Parts
            ind = (np.abs(movingstd - baselinestd) < 1e-9) & (np.abs(output[current]) < np.max(output[current]) / 2)

            # How Many Parts?
            restart = True
            counter = 0
            for i, value in enumerate(ind):
                if value and restart:
                    ItemExp['StartPoint'][counter] = i
                    print('Start: {}\t'.format(i))
                    restart = False
                elif not value and not restart:
                    if ((i - 1) - ItemExp['StartPoint'][counter]) > 1 * output['samplerate']:
                        ItemExp['EndPoint'][counter] = i - 1
                        print('End: {}'.format(i - 1))
                        restart = True
                        counter += 1
                    else:
                        restart = True

            if not restart:
                counter += 1
                ItemExp['EndPoint'][counter] = len(ind)

            ItemExp['StartPoint'] = ItemExp['StartPoint'][np.where(ItemExp['StartPoint'])]
            ItemExp['EndPoint'] = ItemExp['EndPoint'][np.where(ItemExp['EndPoint'])]
            NumberOfPieces = len(ItemExp['StartPoint'])
            ItemExp['Voltage'] = np.zeros(NumberOfPieces)
            ItemExp['ExpValues'] = np.zeros((3, NumberOfPieces))
            ItemExp['STD'] = np.zeros(NumberOfPieces)
            ItemExp['Mean'] = np.zeros(NumberOfPieces)

            # Make It Piecewise
            for i in np.arange(0, NumberOfPieces):
                trace = output[current][ItemExp['StartPoint'][i]:ItemExp['EndPoint'][i]]
                popt, pcov = MakeExponentialFit(np.arange(len(trace)) / output['samplerate'], trace)

                ItemExp['Voltage'][i] = output[sweepedChannel][ItemExp['StartPoint'][i]]
                if popt[0]:
                    ItemExp['ExpValues'][:, i] = popt
                    ItemExp['Mean'][i] = popt[2]
                    try:
                        ItemExp['STD'][i] = np.sqrt(np.diag(pcov))[2]
                    except RuntimeWarning:
                        ItemExp['STD'][i] = np.std(trace)/np.sqrt(len(trace))
                        print('ExpFit: STD of voltage {} failed calculating...'.format(ItemExp['Voltage'][i]))
                else:
                    print('Exponential Fit on for ' + current + ' failed at V=' + str(ItemExp['Voltage'][0]))
                    ItemExp['Mean'][i] = np.mean(trace)
                    ItemExp['STD'][i] = np.std(trace)

                ## FIT THE EXP CALCULATED VALUES ###
                (a, b, sigma_a, sigma_b, b_save) = YorkFit(ItemExp['Voltage'], ItemExp['Mean'],
                                                           1e-12 * np.ones(len(ItemExp['Voltage'])), ItemExp['STD'])
                x_fit = np.linspace(min(Item['Voltage']), max(Item['Voltage']), 1000)
                y_fit = scipy.polyval([b, a], x_fit)
                ItemExp['YorkFitValues'] = {'x_fit': x_fit, 'y_fit': y_fit, 'Yintercept': a, 'Slope': b,
                                            'Sigma_Yintercept': sigma_a,
                                            'Sigma_Slope': sigma_b, 'Parameter': b_save}
                All[current + '_Exp'] = ItemExp

        ##Remove NaNs
        nans = np.isnan(Item['Mean'][:])
        nans = np.logical_not(nans)
        Item['Voltage'] = Item['Voltage'][nans]
        Item['StartPoint'] = Item['StartPoint'][nans]
        Item['EndPoint'] = Item['EndPoint'][nans]
        Item['Mean'] = Item['Mean'][nans]
        Item['STD'] = Item['STD'][nans]

        ## Restrict to set voltage Range
        if UseVoltageRange is not 0:
            print('Voltages Cut')
            ind = np.argwhere((Item['Voltage'] >= UseVoltageRange[0]) & (Item['Voltage'] <= UseVoltageRange[1]))
            Item['Voltage'] = Item['Voltage'][ind].flatten()
            Item['StartPoint'] = Item['StartPoint'][ind].flatten()
            Item['EndPoint'] = Item['EndPoint'][ind].flatten()
            Item['Mean'] = Item['Mean'][ind].flatten()
            Item['STD'] = Item['STD'][ind].flatten()

        ## FIT THE MEAN CALCULATED VALUES ###
        (a, b, sigma_a, sigma_b, b_save) = YorkFit(Item['Voltage'].flatten(), Item['Mean'].flatten(), 1e-9 * np.ones(len(Item['Voltage'])), Item['STD'].flatten())
        x_fit = np.linspace(min(Item['Voltage']), max(Item['Voltage']), 1000)
        y_fit = scipy.polyval([b, a], x_fit)
        Item['YorkFitValues'] = {'x_fit': x_fit, 'y_fit': y_fit, 'Yintercept': a, 'Slope': b, 'Sigma_Yintercept': sigma_a,
                         'Sigma_Slope': sigma_b, 'Parameter': b_save}

        All[current] = Item
        All['Currents'] = currents

    return All

def ExpFunc(x, a, b, c):
    return a * np.exp(-b * x) + c

def YorkFit(X, Y, sigma_X, sigma_Y, r=0):
    N_itermax=10 #maximum number of interations
    tol=1e-15 #relative tolerance to stop at
    N = len(X)
    temp = np.matrix([X, np.ones(N)])
    #make initial guess at b using linear squares

    tmp = np.matrix(Y)*lin.pinv(temp)
    b_lse = np.array(tmp)[0][0]
    #a_lse=tmp(2);
    b = b_lse #initial guess
    omega_X = np.true_divide(1,np.power(sigma_X,2))
    omega_Y = np.true_divide(1, np.power(sigma_Y,2))
    alpha=np.sqrt(omega_X*omega_Y)
    b_save = np.zeros(N_itermax+1) #vector to save b iterations in
    b_save[0]=b

    for i in np.arange(N_itermax):
        W=omega_X*omega_Y/(omega_X+b*b*omega_Y-2*b*r*alpha)

        X_bar=np.sum(W*X)/np.sum(W)
        Y_bar=np.sum(W*Y)/np.sum(W)

        U=X-X_bar
        V=Y-Y_bar

        beta=W*(U/omega_Y+b*V/omega_X-(b*U+V)*r/alpha)

        b=sum(W*beta*V)/sum(W*beta*U)
        b_save[i+1]=b
        if np.abs((b_save[i+1]-b_save[i])/b_save[i+1]) < tol:
            break

    a=Y_bar-b*X_bar
    x=X_bar+beta
    y=Y_bar+b*beta
    x_bar=sum(W*x)/sum(W)
    y_bar=sum(W*y)/sum(W)
    u=x-x_bar
    #%v=y-y_bar
    sigma_b=np.sqrt(1/sum(W*u*u))
    sigma_a=np.sqrt(1./sum(W)+x_bar*x_bar*sigma_b*sigma_b)
    return (a, b, sigma_a, sigma_b, b_save)

def DoubleExpFunc(x, Y0, percentFast, plateau, Kfast, Kslow):
    SpanFast = (Y0-plateau) * percentFast * 0.01
    SpanSlow = (Y0-plateau) *(100-percentFast) * 0.01
    return plateau + ((Y0-plateau) * percentFast * 0.01)*np.exp(-Kfast*x) + ((Y0-plateau) *(100-percentFast) * 0.01)*np.exp(-Kslow*x)

def MakeExponentialFit(xdata, ydata):
    try:
        popt, pcov = curve_fit(ExpFunc, xdata, ydata, method='lm', xtol = 1e-20, ftol = 1e-12, gtol = 1e-20)
        #popt, pcov = curve_fit(DoubleExpFunc, xdata, ydata, p0 = (ydata[0], 50, ydata[-1], 1 / 15, 1 / 60))
        return (popt, pcov)
    except RuntimeError:
        popt = (0,0,0)
        pcov = 0
        return (popt, pcov)

def CutDataIntoVoltageSegments(output):
    sweepedChannel = ''
    if output['type'] == 'ChimeraNotRaw' or (output['type'] == 'Axopatch' and not output['graphene']):
        ChangePoints = np.where(np.diff(output['v1']))[0]
        sweepedChannel = 'v1'
        if len(ChangePoints) is 0:
            print('No voltage sweeps in this file:\n{}'.format(output['filename']))
            return (0,0)
    elif (output['type'] == 'Axopatch' and output['graphene']):
        ChangePoints = np.where(np.diff(output['v1']))[0]
        sweepedChannel = 'v1'
        if len(ChangePoints) is 0:
            ChangePoints = np.where(np.diff(output['v2']))[0]
            if len(ChangePoints) is 0:
                print('No voltage sweeps in this file')
                return (0,0)
            else:
                sweepedChannel = 'v2'
    print('Cutting into Segments...\n{} change points detected in channel {}...'.format(len(ChangePoints), sweepedChannel))
    return (ChangePoints, sweepedChannel)

def CalculatePoreSize(G, L, s):
    return (G+np.sqrt(G*(G+16*L*s/np.pi)))/(2*s)

def CalculateConductance(L, s, d):
    return s*1/(4*L/(np.pi*d**2)+1/d)

def GetSurfaceChargeFromLeeEquation(s, c, D, G, L, A, B, version=1):
    lDu = 1/(8*np.pi*B*D*s)*(
        version * np.sqrt((4*np.pi*A*D**2*s+np.pi*B*D**2*s - 4*B*G*L - 8*np.pi*D*G)**2
                 - 16*np.pi*B*D*s*(np.pi*A*D**3*s - 4*A*D*G*L - 2*np.pi*D**2*G))
        - 4*np.pi*A*D**2*s - np.pi*B*D**2*s + 4*B*G*L + 8*np.pi*D*G
    )
    e = cst.elementary_charge
    Na = cst.Avogadro
    return lDu*(2 * c * Na * 1e3 * e)
def ConductivityFromConductanceModel(L, d, G):
    return G * (4 * L / (np.pi * d ** 2) + 1 / d)

def PlotIV(output, AllData, current = 'i1', unit=1e9, axis = '', WithFit = 1, PoreSize = [0,0], title = '', labeltxt = '', Polynomial=0, color='b', RectificationFactor=0, useEXP=0):
    if useEXP:
        current=current+'_Exp'

    if labeltxt is '':
        labeltxt = str(os.path.split(output['filename'])[1])[:-4]

    ind = np.argsort(AllData[current]['Voltage'])

    #Calculate Rectification Factor
    if RectificationFactor:
        parts = {'pos': AllData[current]['Voltage'] > 0, 'neg': AllData[current]['Voltage'] < 0}
        a={}; b={}; sigma_a={}; sigma_b={}; b_save={};x_values={}

        for part in parts:
            (a[part], b[part], sigma_a[part], sigma_b[part], b_save[part]) = YorkFit(AllData[current]['Voltage'][parts[part]], AllData[current]['Mean'][parts[part]], 1e-12 * np.ones(len(AllData[current]['Voltage'][parts[part]])), AllData[current]['STD'][parts[part]])
            x_values[part] = np.linspace(min(AllData[current]['Voltage'][parts[part]]), max(AllData[current]['Voltage'][parts[part]]), 1000)
            axis.plot(x_values[part], scipy.polyval([b[part], a[part]], x_values[part])*unit, color=color, ls='--')
        factor = b['neg'] / b['pos']
        labeltxt = labeltxt + ', R:{:4.2f}'.format(factor)

    #Polynomial is guide to the eye
    if not Polynomial:
        axis.errorbar(AllData[current]['Voltage'], AllData[current]['Mean']*unit, yerr=AllData[current]['STD']*unit, fmt='o', label=labeltxt, color=color)
    else:
        p = np.polyfit(AllData[current]['Voltage'][ind], AllData[current]['Mean'][ind]*unit, Polynomial)
        axis.errorbar(AllData[current]['Voltage'][ind], AllData[current]['Mean'][ind]*unit, yerr=AllData[current]['STD'][ind]*unit, fmt='o', label=labeltxt, color=color)
        axis.plot(AllData[current]['Voltage'][ind], np.polyval(p, AllData[current]['Voltage'][ind]), color=color)

    axis.set_ylabel('Current ' + current + ' [nA]')
    axis.set_xlabel('Voltage ' + AllData[current]['SweepedChannel'] +' [V]')

    if WithFit:
        axis.set_title(title + '\nG={}S, R={}Ohm'.format(pg.siFormat(AllData[current]['YorkFitValues']['Slope']), pg.siFormat(1/(AllData[current]['YorkFitValues']['Slope']))))
        axis.plot(AllData[current]['YorkFitValues']['x_fit'], AllData[current]['YorkFitValues']['y_fit']*unit, 'r--', label='Linear Fit of '+labeltxt)
    else:
        axis.set_title(title + 'IV Plot')

    if PoreSize[0]:
        textstr = 'Pore Size\nConductance: {}S/m\nLenght: {}m:\ndiameter: {}m'.format(pg.siFormat(PoreSize[0]), pg.siFormat(PoreSize[1]), pg.siFormat(CalculatePoreSize(AllData[current]['YorkFitValues']['Slope'], PoreSize[1], PoreSize[0])))
        axis.text(0.05, 0.95, textstr, transform=axis.transAxes, fontsize=12,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    return axis

def RefinedEventDetection(out, AnalysisResults, signals, limit):
    for sig in signals:
        #Choose correct reference
        if sig is 'i1_Up':
            sig1 = 'i1'
        elif sig is 'i2_Up':
            sig1 = 'i2'
        else:
            sig1 = sig
        if sig is 'i1' or 'i1_Up':
            volt='v1'
        elif sig is 'i2' or 'i2_Up':
            volt='v2'

        if len(AnalysisResults[sig]['RoughEventLocations']) > 1:
            startpoints = np.uint64(AnalysisResults[sig]['RoughEventLocations'][:, 0])
            endpoints = np.uint64(AnalysisResults[sig]['RoughEventLocations'][:, 1])
            localBaseline = AnalysisResults[sig]['RoughEventLocations'][:, 2]
            localVariance = AnalysisResults[sig]['RoughEventLocations'][:, 3]

            CusumBaseline=500
            numberofevents = len(startpoints)
            AnalysisResults[sig]['StartPoints'] = startpoints
            AnalysisResults[sig]['EndPoints'] = endpoints
            AnalysisResults[sig]['LocalBaseline'] = localBaseline
            AnalysisResults[sig]['LocalVariance'] = localVariance
            AnalysisResults[sig]['NumberOfEvents'] = len(startpoints)

            deli = np.zeros(numberofevents)
            dwell = np.zeros(numberofevents)
            fitvalue = np.zeros(numberofevents)
            AllEvents = []

            #####FIX THE UP AND DOWN STORY!! ALL THE PROBLEMS ARE IN HERE...
            for i in range(numberofevents):
                length = endpoints[i] - startpoints[i]
                if out[sig1][startpoints[i]] < 0 and (sig == 'i1_Up' or sig == 'i2_Up'):
                    if length <= limit and length > 3:
                        # Impulsion Fit to minimal value
                        fitvalue[i] = np.max(out[sig1][startpoints[i]+np.uint(1):endpoints[i]-np.uint(1)])
                    elif length > limit:
                        fitvalue[i] = np.mean(out[sig1][startpoints[i]+np.uint(5):endpoints[i]-np.uint(5)])
                    else:
                        fitvalue[i] = np.max(out[sig1][startpoints[i]:endpoints[i]])
                    deli[i] = -(localBaseline[i] - fitvalue[i])

                elif out[sig1][startpoints[i]] < 0 and (sig == 'i1' or sig == 'i2'):
                    if length <= limit and length > 3:
                        # Impulsion Fit to minimal value
                        fitvalue[i] = np.min(out[sig1][startpoints[i]+np.uint(1):endpoints[i]-np.uint(1)])
                    elif length > limit:
                        fitvalue[i] = np.mean(out[sig1][startpoints[i]+np.uint(5):endpoints[i]-np.uint(5)])
                    else:
                        fitvalue[i] = np.min(out[sig1][startpoints[i]:endpoints[i]])
                    deli[i] = -(localBaseline[i] - fitvalue[i])

                elif out[sig1][startpoints[i]] > 0 and (sig == 'i1_Up' or sig == 'i2_Up'):
                    if length <= limit and length > 3:
                        # Impulsion Fit to minimal value
                        fitvalue[i] = np.max(out[sig1][startpoints[i]+np.uint(1):endpoints[i]-np.uint(1)])
                    elif length > limit:
                        fitvalue[i] = np.mean(out[sig1][startpoints[i]+np.uint(5):endpoints[i]-np.uint(5)])
                    else:
                        fitvalue[i] = np.max(out[sig1][startpoints[i]:endpoints[i]])
                    deli[i] = (localBaseline[i] - fitvalue[i])

                elif out[sig1][startpoints[i]] > 0 and (sig == 'i1' or sig == 'i2'):
                    if length <= limit and length > 3:
                        # Impulsion Fit to minimal value
                        fitvalue[i] = np.min(out[sig1][startpoints[i]+np.uint(1):endpoints[i]-np.uint(1)])
                    elif length > limit:
                        fitvalue[i] = np.mean(out[sig1][startpoints[i]+np.uint(5):endpoints[i]-np.uint(5)])
                    else:
                        fitvalue[i] = np.min(out[sig1][startpoints[i]:endpoints[i]])
                    deli[i] = (localBaseline[i] - fitvalue[i])
                else:
                    print('Strange event that has to come from a neighboring universe...Computer will self-destruct in 3s!')
                dwell[i] = (endpoints[i] - startpoints[i]) / out['samplerate']
                AllEvents.append(out[sig1][startpoints[i]:endpoints[i]])
            frac = deli / localBaseline
            dt = np.array(0)
            dt = np.append(dt, np.diff(startpoints) / out['samplerate'])
            numberofevents = len(dt)
            #AnalysisResults[sig]['CusumFits'] = AllFits
            AnalysisResults[sig]['FractionalCurrentDrop'] = frac
            AnalysisResults[sig]['DeltaI'] = deli
            AnalysisResults[sig]['DwellTime'] = dwell
            AnalysisResults[sig]['Frequency'] = dt
            AnalysisResults[sig]['LocalVoltage'] = out[volt][startpoints]
            AnalysisResults[sig]['AllEvents'] = AllEvents
            AnalysisResults[sig]['FitLevel'] = fitvalue
        else:
            AnalysisResults[sig]['NumberOfEvents'] = 0

    return AnalysisResults

def CorrelateTheTwoChannels(AnalysisResults, DelayLimit, Ch1 = 'i1', Ch2 = 'i2'):
    if len(AnalysisResults[Ch1]['RoughEventLocations']) is not 0:
        i1StartP = np.int64(AnalysisResults[Ch1]['StartPoints'][:])
    else:
        i1StartP = []
    if len(AnalysisResults[Ch2]['RoughEventLocations']) is not 0:
        i2StartP = np.int64(AnalysisResults[Ch2]['StartPoints'][:])
    else:
        i2StartP = []

    # Common Events, # Take Longer
    CommonEventsi1Index = np.array([], dtype=np.uint64)
    CommonEventsi2Index = np.array([], dtype=np.uint64)

    for k in i1StartP:
        val = i2StartP[(i2StartP > k - DelayLimit) & (i2StartP < k + DelayLimit)]
        if len(val) == 1:
            CommonEventsi2Index = np.append(CommonEventsi2Index, np.uint64(np.where(i2StartP == val)[0][0]))
            CommonEventsi1Index = np.append(CommonEventsi1Index, np.uint64(np.where(i1StartP == k)[0][0]))
        if len(val) > 1:
            diff = np.absolute(val-k)
            CommonEventsi2Index = np.append(CommonEventsi2Index, np.uint64(np.where(i2StartP == val[np.argmin(diff)]))[0][0])
            CommonEventsi1Index = np.append(CommonEventsi1Index, np.uint64(np.where(i1StartP == k)[0][0]))
    # Only i1
    Onlyi1Indexes = np.delete(range(len(i1StartP)), CommonEventsi1Index)
    #Onlyi1Indexes=[]
    # Only i2
    Onlyi2Indexes = np.delete(range(len(i2StartP)), CommonEventsi2Index)
    #Onlyi2Indexes=[]

    CommonIndexes={}
    CommonIndexes[Ch1]=CommonEventsi1Index
    CommonIndexes[Ch2]=CommonEventsi2Index
    OnlyIndexes={}
    OnlyIndexes[Ch1] = Onlyi1Indexes
    OnlyIndexes[Ch2] = Onlyi2Indexes
    return (CommonIndexes, OnlyIndexes)

def PlotEvent(fig1, t1, i1, t2 = [], i2 = [], fit1 = np.array([]), fit2 = np.array([]), channel = 'i1'):
    if len(t2)==0:
        ax1 = fig1.add_subplot(111)
        ax1.plot(t1, i1*1e9, 'b')
        if len(fit1) is not 0:
            ax1.plot(t1, fit1*1e9, 'y')
        ax1.set_ylabel(channel + ' Current [nA]')
        ax1.set_xlabel(channel + ' Time [s]')
        ax1.ticklabel_format(useOffset=False)
        ax1.ticklabel_format(useOffset=False)
        return ax1
    else:
        ax1 = fig1.add_subplot(211)
        ax2 = fig1.add_subplot(212, sharex=ax1)
        ax1.plot(t1, i1*1e9, 'b')
        if len(fit1) is not 0:
            ax1.plot(t1, fit1*1e9, 'y')
        ax2.plot(t2, i2*1e9, 'r')
        if len(fit2) is not 0:
            ax2.plot(t2, fit2*1e9, 'y')
        ax1.set_ylabel('Ionic Current [nA]')
        #ax1.set_xticklabels([])
        ax2.set_ylabel('Transverse Current [nA]')
        ax2.set_xlabel('Time [s]')
        ax2.ticklabel_format(useOffset=False)
        ax2.ticklabel_format(useOffset=False)
        ax1.ticklabel_format(useOffset=False)
        ax1.ticklabel_format(useOffset=False)
        return ax1, ax2

def SaveAllPlots(CommonIndexes, OnlyIndexes, AnalysisResults, directory, out, buffer, withFit = 1):
    if len(CommonIndexes['i1']) is not 0:
        # Plot All Common Events
        pp = PdfPages(directory + '_SavedEventsCommon.pdf')
        ind1 = np.uint64(CommonIndexes['i1'])
        ind2 = np.uint64(CommonIndexes['i2'])

        t = np.arange(0, len(out['i1']))
        t = t / out['samplerate'] * 1e3
        count=1
        for eventnumber in range(len(ind1)):
            parttoplot = np.arange(AnalysisResults['i1']['StartPoints'][ind1[eventnumber]] - buffer,
                                   AnalysisResults['i1']['EndPoints'][ind1[eventnumber]] + buffer, 1, dtype=np.uint64)
            parttoplot2 = np.arange(AnalysisResults['i2']['StartPoints'][ind2[eventnumber]] - buffer,
                                    AnalysisResults['i2']['EndPoints'][ind2[eventnumber]] + buffer, 1, dtype=np.uint64)

            fit1 = np.concatenate([np.ones(buffer) * AnalysisResults['i1']['LocalBaseline'][ind1[eventnumber]],
                                   np.ones(AnalysisResults['i1']['EndPoints'][ind1[eventnumber]] - AnalysisResults['i1']['StartPoints'][
                                       ind1[eventnumber]]) * (
                                       AnalysisResults['i1']['LocalBaseline'][ind1[eventnumber]] - AnalysisResults['i1']['DeltaI'][ind1[eventnumber]]),
                                   np.ones(buffer) * AnalysisResults['i1']['LocalBaseline'][ind1[eventnumber]]])

            fit2 = np.concatenate([np.ones(buffer) * AnalysisResults['i2']['LocalBaseline'][ind2[eventnumber]],
                                   np.ones(AnalysisResults['i2']['EndPoints'][ind2[eventnumber]] - AnalysisResults['i2']['StartPoints'][
                                       ind2[eventnumber]]) * (
                                       AnalysisResults['i2']['LocalBaseline'][ind2[eventnumber]] - AnalysisResults['i2']['DeltaI'][ind2[eventnumber]]),
                                   np.ones(buffer) * AnalysisResults['i2']['LocalBaseline'][ind2[eventnumber]]])
            if withFit:
                fig = PlotEvent(t[parttoplot], out['i1'][parttoplot], t[parttoplot2], out['i2'][parttoplot2],
                                   fit1=fit1, fit2=fit2)
            else:
                fig = PlotEvent(t[parttoplot], out['i1'][parttoplot], t[parttoplot2], out['i2'][parttoplot2])

            if not np.mod(eventnumber+1,200):
                pp.close()
                pp = PdfPages(directory + '_SavedEventsCommon_' + str(count) + '.pdf')
                count+=1
            pp.savefig(fig)
            print('{} out of {} saved!'.format(str(eventnumber), str(len(ind1))))
            print('Length i1: {}, Fit i1: {}'.format(len(out['i1'][parttoplot]), len(fit1)))
            print('Length i2: {}, Fit i2: {}'.format(len(out['i2'][parttoplot2]), len(fit2)))
            fig.clear()
            plt.close(fig)
        pp.close()

    if len(OnlyIndexes['i1']) is not 0:
        # Plot All i1
        pp = PdfPages(directory + '_SavedEventsOnlyi1.pdf')
        ind1 = np.uint64(OnlyIndexes['i1'])

        t = np.arange(0, len(out['i1']))
        t = t / out['samplerate'] * 1e3
        count=1
        for eventnumber in range(len(ind1)):
            parttoplot = np.arange(AnalysisResults['i1']['StartPoints'][ind1[eventnumber]] - buffer,
                                   AnalysisResults['i1']['EndPoints'][ind1[eventnumber]] + buffer, 1, dtype=np.uint64)

            fit1 = np.concatenate([np.ones(buffer) * AnalysisResults['i1']['LocalBaseline'][ind1[eventnumber]],
                                   np.ones(AnalysisResults['i1']['EndPoints'][ind1[eventnumber]] - AnalysisResults['i1']['StartPoints'][
                                       ind1[eventnumber]]) * (
                                       AnalysisResults['i1']['LocalBaseline'][ind1[eventnumber]] - AnalysisResults['i1']['DeltaI'][ind1[eventnumber]]),
                                   np.ones(buffer) * AnalysisResults['i1']['LocalBaseline'][ind1[eventnumber]]])

            fig = PlotEvent(t[parttoplot], out['i1'][parttoplot], t[parttoplot], out['i2'][parttoplot], fit1=fit1)
            if not np.mod(eventnumber+1,200):
                pp.close()
                pp = PdfPages(directory + '_SavedEventsCommon_' + str(count) + '.pdf')
                count+=1
            pp.savefig(fig)
            print('{} out of {} saved!'.format(str(eventnumber), str(len(ind1))))
            fig.clear()
            plt.close(fig)
        pp.close()

    if len(OnlyIndexes['i2']) is not 0:
        # Plot All i2
        pp = PdfPages(directory + '_SavedEventsOnlyi2.pdf')
        ind1 = np.uint64(OnlyIndexes['i2'])

        t = np.arange(0, len(out['i2']))
        t = t / out['samplerate'] * 1e3
        count=1
        for eventnumber in range(len(ind1)):
            parttoplot = np.arange(AnalysisResults['i2']['StartPoints'][ind1[eventnumber]] - buffer,
                                   AnalysisResults['i2']['EndPoints'][ind1[eventnumber]] + buffer, 1, dtype=np.uint64)

            fit1 = np.concatenate([np.ones(buffer) * AnalysisResults['i2']['LocalBaseline'][ind1[eventnumber]],
                                   np.ones(AnalysisResults['i2']['EndPoints'][ind1[eventnumber]] - AnalysisResults['i2']['StartPoints'][
                                       ind1[eventnumber]]) * (
                                       AnalysisResults['i2']['LocalBaseline'][ind1[eventnumber]] - AnalysisResults['i2']['DeltaI'][ind1[eventnumber]]),
                                   np.ones(buffer) * AnalysisResults['i2']['LocalBaseline'][ind1[eventnumber]]])

            fig = PlotEvent(t[parttoplot], out['i1'][parttoplot], t[parttoplot], out['i2'][parttoplot], fit2=fit1)
            if not np.mod(eventnumber+1,200):
                pp.close()
                pp = PdfPages(directory + '_SavedEventsCommon_' + str(count) + '.pdf')
                count+=1
            pp.savefig(fig)
            print('{} out of {} saved!'.format(str(eventnumber), str(len(ind1))))
            fig.clear()
            plt.close(fig)
        pp.close()

    # Derivative
    if len(CommonIndexes['i1']) is not 0:
        # Plot All i1
        pp = PdfPages(directory + '_i1vsderivi2.pdf')
        ind1 = np.uint64(CommonIndexes['i1'])
        ind2 = np.uint64(CommonIndexes['i2'])

        t = np.arange(0, len(out['i1']))
        t = t / out['samplerate'] * 1e3
        count=1
        for eventnumber in range(len(ind1)):
            parttoplot = np.arange(AnalysisResults['i1']['StartPoints'][ind1[eventnumber]] - buffer,
                                   AnalysisResults['i1']['EndPoints'][ind1[eventnumber]] + buffer, 1, dtype=np.uint64)
            parttoplot2 = np.arange(AnalysisResults['i2']['StartPoints'][ind2[eventnumber]] - buffer,
                                    AnalysisResults['i2']['EndPoints'][ind2[eventnumber]] + buffer, 1, dtype=np.uint64)

            fig = PlotEvent(t[parttoplot], out['i1'][parttoplot], t[parttoplot2][:-1],
                               np.diff(out['i2'][parttoplot2]))

            if not np.mod(eventnumber+1,200):
                pp.close()
                pp = PdfPages(directory + '_SavedEventsCommon_' + str(count) + '.pdf')
                count+=1
            pp.savefig(fig)
            print('{} out of {} saved!'.format(str(eventnumber), str(len(ind1))))
            fig.clear()
            plt.close(fig)
        pp.close()

def PlotRecursiveLPResults(RoughEventLocations, inp, directory, buffer, channel='i2'):
    pp = PdfPages(directory + '_' + channel + '_DetectedEventsFromLPFilter.pdf')
    a=1
    for i in RoughEventLocations['RoughEventLocations']:
        startp = np.uint64(i[0]-buffer*inp['samplerate'])
        endp = np.uint64(i[1]+buffer*inp['samplerate'])
        t = np.arange(startp, endp)
        t = t / inp['samplerate'] * 1e3
        fig = PlotEvent(t, inp[channel][startp:endp], channel=channel)
        pp.savefig(fig)
        print('{} out of {} saved!'.format(str(a), str(len(RoughEventLocations['RoughEventLocations']))))
        a+=1
        fig.clear()
        plt.close(fig)
    pp.close()

def SaveAllAxopatchEvents(AnalysisResults, directory, out, buffer, withFit = 1):
    # Plot All Common Events
    pp = PdfPages(directory + '_SavedEventsAxopatch.pdf')
    t = np.arange(0, len(out['i1']))
    t = t / out['samplerate'] * 1e3

    for eventnumber in range(AnalysisResults['i1']['NumberOfEvents']):
        parttoplot = np.arange(AnalysisResults['i1']['StartPoints'][eventnumber] - buffer,
                               AnalysisResults['i1']['EndPoints'][eventnumber] + buffer, 1, dtype=np.uint64)

        fit1 = np.concatenate([np.ones(buffer) * AnalysisResults['i1']['LocalBaseline'][eventnumber],
                               np.ones(AnalysisResults['i1']['EndPoints'][eventnumber] -
                                       AnalysisResults['i1']['StartPoints'][
                                           eventnumber]) * (
                                   AnalysisResults['i1']['LocalBaseline'][eventnumber] -
                                   AnalysisResults['i1']['DeltaI'][eventnumber]),
                               np.ones(buffer) * AnalysisResults['i1']['LocalBaseline'][eventnumber]])

        if withFit:
            fig = PlotEvent(t[parttoplot], out['i1'][parttoplot], fit1=fit1)
        else:
            fig = PlotEvent(t[parttoplot], out['i1'][parttoplot])
        try:
            pp.savefig(fig)
        except:
            print('Problem at {} !'.format(str(eventnumber)))

        print('{} out of {} saved!'.format(str(eventnumber), str(AnalysisResults['i1']['NumberOfEvents'])))
        #print('Length i1: {}, Fit i1: {}'.format(len(out['i1'][parttoplot]), len(fit1)))
        fig.clear()
        plt.close(fig)
    pp.close()

def MakeIV(filenames, directory, conductance=10, title=''):
    if len(conductance) is 1:
        conductance=np.ones(len(filenames))*conductance
    for idx,filename in enumerate(filenames):
        print(filename)

        # Make Dir to save images
        output = OpenFile(filename)
        if not os.path.exists(directory):
            os.makedirs(directory)

        AllData = MakeIVData(output, delay=2)
        if AllData == 0:
            print('!!!! No Sweep in: ' + filename)
            continue

        # Plot Considered Part
        # figExtracteParts = plt.figure(1)
        # ax1 = figExtracteParts.add_subplot(211)
        # ax2 = figExtracteParts.add_subplot(212, sharex=ax1)
        # (ax1, ax2) = uf.PlotExtractedPart(output, AllData, current = 'i1', unit=1e9, axis = ax1, axis2=ax2)
        # plt.show()
        # figExtracteParts.savefig(directory + os.sep + 'PlotExtracted_' + str(os.path.split(filename)[1])[:-4] + '.eps')
        # figExtracteParts.savefig(directory + os.sep + 'PlotExtracted_' + str(os.path.split(filename)[1])[:-4] + '.png', dpi=150)


        # Plot IV
        if output['graphene']:
            figIV2 = plt.figure(3)
            ax2IV = figIV2.add_subplot(111)
            ax2IV = PlotIV(output, AllData, current='i2', unit=1e9, axis=ax2IV, WithFit=1, title=title[idx])
            figIV2.tight_layout()
            figIV2.savefig(directory + os.sep + str(os.path.split(filename)[1]) + '_IV_i2.png', dpi=300)
            figIV2.savefig(directory + os.sep + str(os.path.split(filename)[1]) + '_IV_i2.eps')
            figIV2.clear()

        figIV = plt.figure(2)
        ax1IV = figIV.add_subplot(111)
        ax1IV = PlotIV(output, AllData, current='i1', unit=1e9, axis=ax1IV, WithFit=1, PoreSize=[conductance[idx], 20e-9], title=title[idx])
        figIV.tight_layout()

        # Save Figures
        figIV.savefig(directory + os.sep + str(os.path.split(filename)[1]) + 'IV_i1.png', dpi=300)
        figIV.savefig(directory + os.sep + str(os.path.split(filename)[1]) + 'IV_i1.eps')
        figIV.clear()

def TwoChannelAnalysis(filenames, outdir, UpwardsOn=0):

    for filename in filenames:
        print(filename)
        inp = OpenFile(filename)
        file = os.sep + str(os.path.split(filename)[1][:-4])

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # Low Pass Event Detection
        AnalysisResults = {}

        if inp['graphene']:
            chan = ['i1', 'i2']  # For looping over the channels
        else:
            chan = ['i1']

        for sig in chan:
            AnalysisResults[sig] = {}
            AnalysisResults[sig]['RoughEventLocations'] = RecursiveLowPassFast(inp[sig], pm.coefficients[sig],
                                                                                 inp['samplerate'])
            if UpwardsOn:  # Upwards detection can be turned on or off
                AnalysisResults[sig + '_Up'] = {}
                AnalysisResults[sig + '_Up']['RoughEventLocations'] = RecursiveLowPassFastUp(inp[sig],
                                                                                               pm.coefficients[sig],
                                                                                               inp['samplerate'])

        # f.PlotRecursiveLPResults(AnalysisResults['i2'], inp, directory, 1e-3, channel='i2')
        # f.PlotRecursiveLPResults(AnalysisResults['i2_Up'], inp, directory+'UP_', 1e-3, channel='i2')

        # Refine the Rough Event Detection done by the LP filter and Add event infos
        AnalysisResults = RefinedEventDetection(inp, AnalysisResults, signals=chan,
                                                  limit=pm.MinimalFittingLimit * inp['samplerate'])

        # Correlate the two channels
        if inp['graphene']:
            (CommonIndexes, OnlyIndexes) = f.CorrelateTheTwoChannels(AnalysisResults, 10e-3 * inp['samplerate'])
            print('\n\nAnalysis Done...\nThere are {} common events\n{} Events on i1 only\n{} Events on i2 only'.format(
                len(CommonIndexes['i1']), len(OnlyIndexes['i1']), len(OnlyIndexes['i2'])))
            # Plot The Events
            f.SaveAllPlots(CommonIndexes, OnlyIndexes, AnalysisResults, directory, inp, pm.PlotBuffer)
        else:
            # Plot The Events
            f.SaveAllAxopatchEvents(AnalysisResults, directory, inp, pm.PlotBuffer)

def SaveToHDF5(inp_file, AnalysisResults, coefficients, outdir):
    file = str(os.path.split(inp_file['filename'])[1][:-4])
    f = h5py.File(outdir + file + '_OriginalDB.hdf5', "w")
    general = f.create_group("General")
    general.create_dataset('FileName', data=inp_file['filename'])
    general.create_dataset('Samplerate', data=inp_file['samplerate'])
    general.create_dataset('Machine', data=inp_file['type'])
    general.create_dataset('TransverseRecorded', data=inp_file['graphene'])
    general.create_dataset('TimeFileWritten', data=inp_file['TimeFileWritten'])
    general.create_dataset('TimeFileLastModified', data=inp_file['TimeFileLastModified'])
    general.create_dataset('ExperimentDuration', data=inp_file['ExperimentDuration'])

    segmentation_LP = f.create_group("LowPassSegmentation")
    for k,l in AnalysisResults.items():
        set1 = segmentation_LP.create_group(k)
        lpset1 = set1.create_group('LowPassSettings')
        for o, p in coefficients[k].items():
             lpset1.create_dataset(o, data=p)
        for m, l in AnalysisResults[k].items():
            if m is 'AllEvents':
                eventgroup = set1.create_group(m)
                for i, val in enumerate(l):
                    eventgroup.create_dataset('{:09d}'.format(i), data=val)
            elif m is 'Cusum':
                eventgroup = set1.create_group(m)
                for i1, val1 in enumerate(AnalysisResults[k]['Cusum']):
                    cusevent = eventgroup.create_group('{:09d}'.format(i1))
                    cusevent.create_dataset('NumberLevels', data=np.uint64(len(AnalysisResults[k]['Cusum'][i1]['levels'])))
                    if len(AnalysisResults[k]['Cusum'][i1]['levels']):
                        cusevent.create_dataset('up', data=AnalysisResults[k]['Cusum'][i1]['up'])
                        cusevent.create_dataset('down', data=AnalysisResults[k]['Cusum'][i1]['down'])
                        cusevent.create_dataset('both', data=AnalysisResults[k]['Cusum'][i1]['both'])
                        cusevent.create_dataset('fit', data=AnalysisResults[k]['Cusum'][i1]['fit'])
                        # 0: level number, 1: current, 2: length, 3: std
                        cusevent.create_dataset('levels_current', data=AnalysisResults[k]['Cusum'][i1]['levels'][1])
                        cusevent.create_dataset('levels_length', data=AnalysisResults[k]['Cusum'][i1]['levels'][2])
                        cusevent.create_dataset('levels_std', data=AnalysisResults[k]['Cusum'][i1]['levels'][3])
            else:
                set1.create_dataset(m, data=l)

def PlotExtractedPart(output, AllData, current = 'i1', unit=1e9, axis = '', axis2 = ''):
    time = np.arange(0, len(output[current])) / output['samplerate']
    axis.plot(time, output[current] * unit, 'b', label=str(os.path.split(output['filename'])[1])[:-4])
    axis.set_ylabel('Current ' + current + ' [nA]')
    axis.set_title('Time Trace')
    for i in range(0, len(AllData[current]['StartPoint'])):
        axis.plot(time[AllData[current]['StartPoint'][i]:AllData[current]['EndPoint'][i]],
                 output[current][AllData[current]['StartPoint'][i]:AllData[current]['EndPoint'][i]] * unit, 'r')
    axis2.plot(time, output[AllData[current]['SweepedChannel']], 'b', label=str(os.path.split(output['filename'])[1])[:-4])
    axis2.set_ylabel('Voltage ' + AllData[current]['SweepedChannel'] + ' [V]')
    axis2.set_xlabel('Time')
    return (axis, axis2)

def xcorr(x, y, k, normalize=True):
    n = x.shape[0]

    # initialize the output array
    out = np.empty((2 * k) + 1, dtype=np.double)
    lags = np.arange(-k, k + 1)

    # pre-compute E(x), E(y)
    mu_x = x.mean()
    mu_y = y.mean()

    # loop over lags
    for ii, lag in enumerate(lags):

        # use slice indexing to get 'shifted' views of the two input signals
        if lag < 0:
            xi = x[:lag]
            yi = y[-lag:]
        elif lag > 0:
            xi = x[:-lag]
            yi = y[lag:]
        else:
            xi = x
            yi = y

        # x - mu_x; y - mu_y
        xdiff = xi - mu_x
        ydiff = yi - mu_y

        # E[(x - mu_x) * (y - mu_y)]
        out[ii] = xdiff.dot(ydiff) / n

        # NB: xdiff.dot(ydiff) == (xdiff * ydiff).sum()

    if normalize:
        # E[(x - mu_x) * (y - mu_y)] / (sigma_x * sigma_y)
        out /=  np.std(x) * np.std(y)

    return lags, out

def cusum_detection(data, threshhold, minlength, maxstates):
    print(threshhold)
    s = timer()
    logp = 0  # instantaneous log-likelihood for positive jumps
    logn = 0  # instantaneous log-likelihood for negative jumps
    cpos = np.zeros(len(data), dtype='float64')  # cumulative log-likelihood function for positive jumps
    cneg = np.zeros(len(data), dtype='float64')  # cumulative log-likelihood function for negative jumps
    gpos = np.zeros(2, dtype='float64')  # decision function for positive jumps
    gneg = np.zeros(2, dtype='float64')  # decision function for negative jumps
    edges = np.array([0],
                     dtype='int64')  # initialize an array with the position of the first subevent - the start of the event
    real_start = np.array([],
                          dtype='int64')  # initialize an array with the position of the first subevent - the start of the event
    real_end = np.array([],
                        dtype='int64')  # initialize an array with the position of the first subevent - the start of the event
    real_Depth = np.array([],
                          dtype='int64')  # initialize an array with the position of the first subevent - the start of the event

    anchor = 0  # the last detected change
    length = len(data)
    data_std = np.std(data)
    h = threshhold / data_std
    k = 0
    nStates = np.uint64(0)
    varM = data[0]
    varS = 0
    mean = data[0]
    print('length data =' + str(length))
    v = np.zeros(length, dtype='float64')
    while k < length - 100:
        k += 1
        if nStates == 0:
            variance = np.var(data[anchor:k])  # initial params for pattern region
        mean = np.mean(data[anchor:k])
        if variance == 0: break
        logp = threshhold / variance * (data[
                                            k] - mean - threshhold / 2.)  # instantaneous log-likelihood for current sample assuming local baseline has jumped in the positive direction
        logn = -threshhold / variance * (data[
                                             k] - mean + threshhold / 2.)  # instantaneous log-likelihood for current sample assuming local baseline has jumped in the negative direction
        cpos[k] = cpos[k - 1] + logp  # accumulate positive log-likelihoods
        cneg[k] = cneg[k - 1] + logn  # accumulate negative log-likelihoods
        gpos[1] = max(gpos[0] + logp, 0)  # accumulate or reset positive decision function
        gneg[1] = max(gneg[0] + logn, 0)  # accumulate or reset negative decision function
        if (gpos[1] > h or gneg[1] > h):

            if (gpos[1] > h):  # significant positive jump detected
                jump = np.uint64(anchor + np.argmin(cpos[anchor:k + 1]))  # find the location of the start of the jump
                if jump - edges[np.uint64(nStates)] > minlength and np.abs(data[np.uint64(jump + minlength)] - data[jump]) > threshhold / 4:
                    edges = np.append(edges, jump)
                    nStates += 1
                    # print('EVENT!!!!! at ='+str(self.t[jump]))
                    anchor = k  # no data meaning at bad points!
                    # away from bad point more!
                    cpos[0:len(cpos)] = 0  # reset all decision arrays
                    cneg[0:len(cneg)] = 0
                    gpos[0:len(gpos)] = 0
                    gneg[0:len(gneg)] = 0
            if (gneg[1] > h):  # significant negative jump detected
                jump = anchor + np.argmin(cneg[anchor:k + 1])
                if jump - edges[np.uint64(nStates)] > minlength and np.abs(data[np.uint64(jump + minlength)] - data[jump]) > threshhold / 4:
                    edges = np.append(edges, jump)
                    nStates += 1
                    # print('EVENT!!!!! at ='+str(self.t[jump] ))
                    anchor = k  # no data meaning at bad points!
                    # away from bad point more!
                    cpos[0:len(cpos)] = 0  # reset all decision arrays
                    cneg[0:len(cneg)] = 0
                    gpos[0:len(gpos)] = 0
                    gneg[0:len(gneg)] = 0

        gpos[0] = gpos[1]
        gneg[0] = gneg[1]
        if maxstates > 0:
            if nStates > maxstates:
                print('too sensitive')
                nStates = 0
                k = 0
                threshhold = threshhold * 1.1
                h = h * 1.1
                logp = 0  # instantaneous log-likelihood for positive jumps
                logn = 0  # instantaneous log-likelihood for negative jumps
                cpos = np.zeros(len(data), dtype='float64')  # cumulative log-likelihood function for positive jumps
                cneg = np.zeros(len(data), dtype='float64')  # cumulative log-likelihood function for negative jumps
                gpos = np.zeros(2, dtype='float64')  # decision function for positive jumps
                gneg = np.zeros(2, dtype='float64')  # decision function for negative jumps
                edges = np.array([0],
                                 dtype='int64')  # initialize an array with the position of the first subevent - the start of the event
                anchor = 0  # the last detected change
                length = len(data)
                mean = data[0]
                nStates = 0
                mean = data[0]
    edges = np.append(edges, len(data) - 1)  # mark the end of the event as an edge
    nStates += 1

    cusum = dict()
    print('Events = ' + str([edges]))
    for i in range(len(edges) - 1):
        real_start = np.append(real_start, edges[i])
        real_end = np.append(real_end, edges[i + 1])
        real_Depth = np.append(real_Depth, np.mean(data[np.uint64(edges[i]):np.uint64(edges[i + 1])]))

    cusum['Real_Start'] = real_start
    cusum['Real_End'] = real_end
    cusum['Real_Depth'] = real_Depth
    print('Real Start =' + str(cusum['Real_Start']))
    print('Real End =' + str(cusum['Real_End']))
    cusum['CurrentLevels'] = [np.average(data[np.uint64(edges[i]) + minlength:np.uint64(edges[i + 1])]) for i in
                              range(np.uint64(nStates))]  # detect current levels during detected sub-event
    print('Edges[-1] = ' + str(edges[-1]))
    cusum['EventDelay'] = edges  # locations of sub-events in the data
    cusum['Threshold'] = threshhold  # record the threshold used
    print('Event = ' + str(cusum['EventDelay']))
    cusum['jumps'] = np.diff(cusum['CurrentLevels'])

    e = timer()
    print('cusum took = ' + str(e - s) + 's')
    return cusum

def detect_cusum(data, var, threshhold, minlength, maxstates):
    logp = 0  # instantaneous log-likelihood for positive jumps
    logn = 0  # instantaneous log-likelihood for negative jumps
    cpos = np.zeros(len(data), dtype='float64')  # cumulative log-likelihood function for positive jumps
    cneg = np.zeros(len(data), dtype='float64')  # cumulative log-likelihood function for negative jumps
    gpos = np.zeros(len(data), dtype='float64')  # decision function for positive jumps
    gneg = np.zeros(len(data), dtype='float64')  # decision function for negative jumps
    edges = np.array([0],
                     dtype='int64')  # initialize an array with the position of the first subevent - the start of the event
    anchor = 0  # the last detected change
    length = len(data)
    mean = data[0]
    h = threshhold / var
    k = self.safety_reg - self.control1
    nStates = 0
    varM = data[0]
    varS = 0
    mean = data[0]
    print('length data =' + str(length))
    v = np.zeros(length, dtype='float64')
    while k < length - 100:
        k += 1
        if nStates == 0: variance = np.var(data[anchor:k])
        mean = np.mean(data[anchor:k])
        logp = threshhold / variance * (data[
                                            k] - mean - threshhold / 2.)  # instantaneous log-likelihood for current sample assuming local baseline has jumped in the positive direction
        logn = -threshhold / variance * (data[
                                             k] - mean + threshhold / 2.)  # instantaneous log-likelihood for current sample assuming local baseline has jumped in the negative direction
        cpos[k] = cpos[k - 1] + logp  # accumulate positive log-likelihoods
        cneg[k] = cneg[k - 1] + logn  # accumulate negative log-likelihoods
        gpos[k] = max(gpos[k - 1] + logp, 0)  # accumulate or reset positive decision function
        gneg[k] = max(gneg[k - 1] + logn, 0)  # accumulate or reset negative decision function
        if (gpos[k] > h or gneg[k] > h):

            if (gpos[k] > h):  # significant positive jump detected
                jump = anchor + np.argmin(cpos[anchor:k + 1])  # find the location of the start of the jump
                if jump - edges[nStates] > minlength and np.abs(data[jump + minlength] - data[jump]) > threshhold / 4:
                    edges = np.append(edges, jump)
                    nStates += 1
                    print('EVENT!!!!! at =' + str((self.control1 + jump) / self.out['samplerate']))
                    anchor = k  # no data meaning at bad points!
                    # away from bad point more!
                    cpos[0:len(cpos)] = 0  # reset all decision arrays
                    cneg[0:len(cneg)] = 0
                    gpos[0:len(gpos)] = 0
                    gneg[0:len(gneg)] = 0
            if (gneg[k] > h):  # significant negative jump detected
                jump = anchor + np.argmin(cneg[anchor:k + 1])
                if jump - edges[nStates] > minlength and np.abs(data[jump + minlength] - data[jump]) > threshhold / 4:
                    edges = np.append(edges, jump)
                    nStates += 1
                    print('EVENT!!!!! at =' + str((self.control1 + jump) / self.out['samplerate']))
                    anchor = k  # no data meaning at bad points!
                    # away from bad point more!
                    cpos[0:len(cpos)] = 0  # reset all decision arrays
                    cneg[0:len(cneg)] = 0
                    gpos[0:len(gpos)] = 0
                    gneg[0:len(gneg)] = 0

        if maxstates > 0:
            if nStates > 10:
                print('too sensitive')
                nStates = 0
                k = 0
                threshhold = threshhold * 1.1
                h = h * 1.1
                logp = 0  # instantaneous log-likelihood for positive jumps
                logn = 0  # instantaneous log-likelihood for negative jumps
                cpos = np.zeros(len(data), dtype='float64')  # cumulative log-likelihood function for positive jumps
                cneg = np.zeros(len(data), dtype='float64')  # cumulative log-likelihood function for negative jumps
                gpos = np.zeros(len(data), dtype='float64')  # decision function for positive jumps
                gneg = np.zeros(len(data), dtype='float64')  # decision function for negative jumps
                edges = np.array([0],
                                 dtype='int64')  # initialize an array with the position of the first subevent - the start of the event
                anchor = 0  # the last detected change
                length = len(data)
                mean = data[0]
                nStates = 0
                mean = data[0]
    edges = np.append(edges, len(data))  # mark the end of the event as an edge
    nStates += 1

    cusum = dict()
    cusum['CurrentLevels'] = [np.average(data[edges[i] + minlength:edges[i + 1]]) for i in
                              range(nStates)]  # detect current levels during detected sub-events
    cusum['EventDelay'] = edges  # locations of sub-events in the data
    cusum['Threshold'] = threshhold  # record the threshold used
    cusum['jumps'] = np.diff(cusum['CurrentLevels'])
    # self.__recordevent(cusum)

    return cusum

def CUSUM(input, delta, h, startp):
    Nd = 0
    kd = []
    krmv = []
    krmvdown = []
    both = []
    k0 = 0
    k = 0
    l = len(input)
    m = np.zeros(l)
    m[k0] = input[k0]
    v = np.zeros(l)
    sp = np.zeros(l)
    Sp = np.zeros(l)
    gp = np.zeros(l)
    sn = np.zeros(l)
    Sn = np.zeros(l)
    gn = np.zeros(l)
    mc = []
    while k < l-1:
        k += 1
        m[k] = np.mean(input[k0:k+1])
        v[k] = np.var(input[k0:k+1])

        sp[k] = delta / v[k] * (input[k] - m[k] - delta / 2)
        sn[k] = -delta / v[k] * (input[k] - m[k] + delta / 2)

        Sp[k] = Sp[k - 1] + sp[k]
        Sn[k] = Sn[k - 1] + sn[k]

        gp[k] = np.max([gp[k - 1] + sp[k], 0])
        gn[k] = np.max([gn[k - 1] + sn[k], 0])

        if gp[k] > h or gn[k] > h:
            if gp[k] > h:
                kmin = np.argmin(Sp[k0:k])
                krmv.append(kmin + k0 - 1)
                both.append(kmin + k0 - 1)
            else:
                kd.append(k)
                kmin = np.argmin(Sn[k0:k])
                krmvdown.append(kmin + k0 - 1)
                both.append(kmin + k0 - 1)
            k0 = k
            m[k0] = input[k0]
            v[k0] = 0
            sp[k0] = 0
            Sp[k0] = 0
            gp[k0] = 0
            sn[k0] = 0
            Sn[k0] = 0
            gn[k0] = 0
            Nd = Nd + 1
    '''
    if Nd == 0:
        mc = np.mean(input) * np.ones(k)
    elif Nd == 1:
        mc = [m[krmv[0]] * np.ones(krmv[0]), m[k] * np.ones(k - krmv[0])]
    else:
        mc.append(m[krmv[0]] * np.ones(krmv[0]))
        for ii in range(1, Nd-1):
            mc.append(m[krmv[ii]] * np.ones(krmv[ii] - krmv[ii - 1]))
        mc.append(m[k] * np.ones(k - krmv[Nd-1]))
    '''
    ## Calculate Inner Levels
    if len(both) >= 2:
        levels = [np.zeros(len(both)-1), np.zeros(len(both)-1),
                  np.zeros(len(both)-1), np.zeros(len(both)-1)]
        # 0: level number, 1: current, 2: length, 3: std
        fit = np.array([])
        for ind, g in enumerate(both[:-1]):
            levels[0][ind] = ind
            levels[1][ind] = np.mean(input[g:both[ind + 1]])
            levels[2][ind] = both[ind + 1] - g
            levels[3][ind] = np.std(input[g:both[ind + 1]])
            np.append(fit, levels[1][ind] * np.ones(np.uint64(levels[2][ind])))
    else:
        levels = []
        fit = []
    out = {'up': krmv+startp, 'down': krmvdown+startp, 'both': both+startp, 'levels':levels, 'fit': fit}
    return out

def CondPoreSizeTick(x, pos):
    formCond = EngFormatter(unit='S')
    formPoreSize = EngFormatter(unit='m')
    outstring = '{}, {}'.format(formCond.format_eng(x), formPoreSize.format_eng(CalculatePoreSize(x, 1e-9, 10)))
    return outstring

def EredoxBefore14062018():
    Eredox = np.array([31.7, 82.9, 135, 185], dtype=float)
    Eredox = Eredox * 1e-3
    return (Eredox-Eredox[0])

def DebyeHuckelParam(c, T = 20):
    ## Valid Only for same valency and c1 = c2
    epsilon_r = 87.740-0.40008*T+9.398e-4*T**2-1.410e-6*T**3
    print('Dielectric of Water is {}'.format(epsilon_r))
    n = c*1e3*cst.Avogadro
    formCond = EngFormatter(unit='m')
    k = np.sqrt((cst.e**2*2*n)/(cst.k*(T+273.15)*cst.epsilon_0*epsilon_r))
    return [k, 1/k]#, '{}'.format(formCond.format_eng(1/k))]

def ElectrostaticPotential(T, sigma, Tau=5e18, pK=7.9, pH=11):
    return cst.k * T / cst.e * (np.log(-sigma / (cst.e * Tau + sigma)) + (pK - pH) * np.log(10))

def GetTempRedoxPotential(Cdiff=1e-3/1e-1, T=293.15):
    return cst.R * T / cst.physical_constants['Faraday constant'][0] * np.log(Cdiff)

def GetRedoxError(cmax, cmin, ErrorPercent = 0.1, T = 293.15):
    cmaxErr = cmax*ErrorPercent
    cminErr = cmin*ErrorPercent
    return cst.R * T / cst.physical_constants['Faraday constant'][0] * np.sqrt((1/cmax**2 * cmaxErr**2 + 1/cmin**2 * cminErr**2))

##Result: Negative ion divided by positive ion
#Goldman–Hodgkin–Katz voltage equation
def GetIonSelectivityWithoutPotential(c_trans, c_cis, Erev, T):
    n_trans = c_trans#*1e3#*cst.Avogadro
    n_cis = c_cis#*1e3#*cst.Avogadro
    A = Erev * cst.physical_constants['Faraday constant'][0] / (cst.R*T)
    return -(n_trans-n_cis*np.exp(-A))*(1-np.exp(A))/((n_trans-n_cis*np.exp(A))*(1-np.exp(-A)))

def GetIonSelectivityWithPotential(c_trans, c_cis, Erev, T, phi):
    sel1 = GetIonSelectivityWithoutPotential(c_trans, c_cis, Erev, T)
    return sel1*np.exp((-2*phi*cst.physical_constants['Faraday constant'][0])/(cst.R * T))

def Selectivity(ConcGradient, V):
    return V * cst.physical_constants['Faraday constant'][0]/((cst.R*293.15) * np.log(ConcGradient))
