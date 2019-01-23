import os
import platform
import shelve
from tkinter import filedialog
import h5py
import numpy as np
import pyabf
import scipy
import scipy.signal as sig
from PyQt5 import QtGui
from scipy import io
from scipy import signal
import Functions


def ImportABF(datafilename):
    abf = pyabf.ABF(datafilename)
    #abf.info()  # shows what is available
    #output={'type': 'Clampfit', 'graphene': 0, 'samplerate': abf.pointsPerSec, 'i1': -20000./65536 * abf.dataY, 'v1': abf.dataC, 'filename': datafilename}
    output = {'type': 'Clampfit', 'graphene': 0, 'samplerate': abf.dataRate, 'i1': abf.data[0] * 1e-12,
              'v1': abf.data[1], 'filename': datafilename}
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
    if 'blockLength' in matfile:
        blockLength=np.int32(matfile['blockLength'])
    else:
        blockLength=1048576

    closedloop_gain = TIAgain * preADCgain
    bitmask = (2 ** 16 - 1) - (2 ** (16 - ADCbits) - 1)
    data = -ADCvref + (2 * ADCvref) * (data & bitmask) / 2 ** 16
    data = (data / closedloop_gain + currentoffset)
    data.shape = [data.shape[1], ]
    output = {'matfilename': str(os.path.splitext(datafilename)[0]),'i1raw': data, 'v1': np.float64(matfile['SETUP_biasvoltage']), 'samplerateRaw': np.int64(samplerate), 'type': 'ChimeraRaw', 'filename': datafilename, 'graphene': 0, 'blockLength':blockLength}
    return output


def ImportChimeraData(datafilename):
    matfile = io.loadmat(str(os.path.splitext(datafilename)[0]))
    samplerate = matfile['ADCSAMPLERATE']
    if samplerate<4e6:
        data = np.fromfile(datafilename, np.dtype('float64'))
        buffersize = int(matfile['DisplayBuffer'])
        out = Functions.Reshape1DTo2D(data, buffersize)
        output = {'i1': out['i1'], 'v1': out['v1'], 'samplerate':float(samplerate), 'type': 'ChimeraNotRaw', 'filename': datafilename,'graphene': 0}
    else:
        output = ImportChimeraRaw(datafilename)
    return output

def ImportCSV(datafilename):
    """ 
     Function used to import files in other classic formats than .abf, .dat, .log, have here new simple
    functions ImportCSV to deal with .txt and .csv data files if OpenFile is called.
    Here the functions read the txt and csv files from the 7th row. We surrpose that the 6 first ones
    contain the settings and as the first column in the file being the current measured and the second column
    being the voltage measured (Supposedly constant).
    
    Returns a dictionary output with in kews: 
        'type' : string with the type of file read, here a text file
        'graphene' : boolean indicating if the recording was made with a transverse current measurement (1 or True) or not (0 or False)
        'i1' a list of the currents
        'v1' a list of the voltages
        'samplerate' float of sampling frequency
        'filename' : string with the input filename
    """
    x=np.loadtxt(datafilename, np.dtype('>f4'), skiprows=1)
#    f=open(datafilename, encoding='utf-8')
#    for i in range(0, ): 
#        #The 7 first lines in the file contain informations on the parameters, change here 
#        a=str(f.readline())
#        print(a)
#        if 'Acquisition' in a or 'Sample Rate' in a:
#            samplerate=int(''.join(i for i in a if i.isdigit()))/1000
#        if 'FEMTO preamp Bandwidth' in a:
#            femtoLP=int(''.join(i for i in a if i.isdigit())) 
    i1 = np.array(x[:,1])#reads the current values in the first column of the txt file
    v1 = np.array(x[:,3]) #reads the voltage values in the second column of the txt file
    output={'type': 'TextFile', 'graphene': 0, 'samplerate': 100e3, 'i1': i1, 'v1': v1, 'filename': datafilename}
    return output 


def OpenFile(filename = '', ChimeraLowPass = 10e3,approxImpulseResponse=False,Split=False,verbose=False):
    """ 
    Function used to read data. It extracts the currents and voltage signals from the file in input
    by calling the import function corresponding to the file format.
    
    Returns a dictionary output with in keys:
        'i1' a list of the currents
        'v1' a list of the voltages
        'samplerate' float of sampling frequency
        
    """
    
    if ChimeraLowPass==None:
        ChimeraLowPass=10e3
    if filename == '':
        datafilename = QtGui.QFileDialog.getOpenFileName()
        datafilename=datafilename[0]
        print(datafilename)
    else:
        datafilename=filename

    if verbose:
        print('Loading file... ' +filename)

    if datafilename[-3::] == 'dat':
        isdat = 1
        output = ImportAxopatchData(datafilename)
    elif datafilename[-3::] == 'log':
        isdat = 0
        output = ImportChimeraData(datafilename)
        if output['type'] is 'ChimeraRaw':  # Lowpass and downsample
            if verbose:
                print('length: ' + str(len(output['i1raw'])))
            Wn = round(2 * ChimeraLowPass / output['samplerateRaw'], 4)  # [0,1] nyquist frequency
            b, a = signal.bessel(4, Wn, btype='low', analog=False)  # 4-th order digital filter


            if approxImpulseResponse:
                z, p, k = signal.tf2zpk(b, a)
                eps = 1e-9
                r = np.max(np.abs(p))
                approx_impulse_len = int(np.ceil(np.log(eps) / np.log(r)))
                Filt_sig=(signal.filtfilt(b, a, output['i1raw'], method='gust', irlen=approx_impulse_len))
            else:
                Filt_sig=(signal.filtfilt(b, a, output['i1raw'], method='gust'))

            ds_factor = np.ceil(output['samplerateRaw'] / (5 * ChimeraLowPass))
            output['i1'] = scipy.signal.resample(Filt_sig, int(len(output['i1raw']) / ds_factor))
            output['samplerate'] = output['samplerateRaw'] / ds_factor
            output['v1'] = output['v1']*np.ones(len(output['i1']))
            if verbose:
                print('Samplerate after filtering:' + str(output['samplerate']))
                print('new length: ' + str(len(output['i1'])))

            if Split:
                splitNr=len(output['i1raw'])/output['blockLength']
                assert(splitNr.is_integer())
                rawSplit=np.array_split(output['i1raw'], splitNr)
                Filt_sigSplit=[]
                for raw in rawSplit:
                    if approxImpulseResponse:
                        z, p, k = signal.tf2zpk(b, a)
                        eps = 1e-9
                        r = np.max(np.abs(p))
                        approx_impulse_len = int(np.ceil(np.log(eps) / np.log(r)))
                        signalSplit=signal.filtfilt(b, a, raw, method='gust', irlen=approx_impulse_len)
                    else:
                        signalSplit=signal.filtfilt(b, a, raw, method = 'gust')
                    Filt_sigSplit.append(scipy.signal.resample(signalSplit, int(len(raw) / ds_factor)))
                output['i1Cut']=Filt_sigSplit

        if max(abs(output['i1']))>1e-6:
            if verbose:
                print('converting to SI units')
            output['i1']=1e-9*output['i1']
            output['v1']=1e-3*output['v1']
            
    elif datafilename[-3::] == 'abf':
        output = ImportABF(datafilename)
        if verbose:
            print('length: ' + str(len(output['i1'])))
            
    elif datafilename[-3::] == 'csv' or datafilename[-3::] == 'txt':
        output = ImportCSV(datafilename)
        if verbose:
            print('lenght: ' + str(len(output['i1'])))

    st = os.stat(datafilename)
    if platform.system() == 'Darwin':
        if verbose:
            print('Platform is ' + platform.system())
        output['TimeFileWritten'] = st.st_birthtime
        output['TimeFileLastModified'] = st.st_mtime
        output['ExperimentDuration'] = st.st_mtime - st.st_birthtime
    elif platform.system() == 'Windows':
        if verbose:
            print('Platform is Windows')
        output['TimeFileWritten'] = st.st_ctime
        output['TimeFileLastModified'] = st.st_mtime
        output['ExperimentDuration'] = st.st_mtime - st.st_ctime
    else:
        if verbose:
            print('Platform is ' + platform.system() +
                ', might not get accurate results.')
        try:
            output['TimeFileWritten'] = st.st_ctime
            output['TimeFileLastModified'] = st.st_mtime
            output['ExperimentDuration'] = st.st_mtime - st.st_ctime
        except:
            raise Exception('Platform not detected')

    return output

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

def SaveVariables(savename, **kwargs):
    if os.path.isdir(savename):
        savefile=os.path.join(savename,os.path.basename(savename)+'_Events')
    else:
        #cut of .dat extension
        if savename.lower().endswith('.dat'):
            savename = savename[0:-4]

        #Check if file already exists, otherwise popup dialog
        if os.path.isfile(savename + '.dat'):
            #root = tkinter.Tk()
            #root.withdraw()
            savename = filedialog.asksaveasfile(mode='w', defaultextension=".dat")
            if savename is None:  # asksaveasfile return `None` if dialog closed with "cancel".
                return
            savefile=base=os.path.splitext(savename.name)[0]
            # raise IOError('File ' + savename + '.dat already exists.')
        else:
            savefile = savename

    #Check if directory exists
    directory = os.path.dirname(savefile)
    if not os.path.exists(directory):
        os.makedirs(directory)

    shelfFile=shelve.open(savefile)
    for arg_name in kwargs:
        shelfFile[arg_name]=kwargs[arg_name]
    shelfFile.close()
    print('saved as: ' + savefile + '.dat')

def LoadVariables(loadname, variableName):
    if not isinstance(loadname, str):
        raise Exception('The second argument must be a string')
    #cut of .dat extension
    if loadname.lower().endswith('.dat'):
        loadname = loadname[0:-4]

    if not os.path.isfile(loadname + '.dat'):
        raise Exception('File does not exist')

    shelfFile = shelve.open(loadname)

    try:
        Variable  = shelfFile[variableName]
    except KeyError as e:
        message = 'Key ' + variableName + ' does not exist, available Keys: \n' + "\n".join(list(shelfFile.keys()))
        print(message)
        raise
    else:
        shelfFile.close()
        print('Loaded  ' + variableName + 'from ' + loadname + '.dat')
        return Variable
