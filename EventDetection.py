import NanoporeClasses as NC
import shelve
import os
import glob
from time import sleep
import Functions
from pprint import pprint
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
import numpy as np
import sys
import argparse
import datetime
from tkinter.filedialog import askopenfilenames,askdirectory
import timeit
import LoadData


timeInSec = EngFormatter(unit='s', places=2)
Amp = EngFormatter(unit='A', places=2)

#Default parameters
extension = '*.log'
coefficients = {'a': 0.999,
                'E': 0,
                'S': 5,
                'maxEventLength': 200e-3,  # maximal event length to be considered an event
                'minEventLength': 600e-6,  # maximal event length to be considered an impulse
                'fitLength': 3e-3,  # minimal event length to be fitted for impulses
                'dt': 25,  # go back dt points for fitting of impulses
                'hbook': 1,
                'delta': 0.2e-9,
                'deltaRel': None,  # If set overrides delta, calculates delta relative to baseline
                'ChimeraLowPass': 10e3}

def GetParameters():
    print("Usage:")
    print("run(inputData, newExtension=None, newCoefficients={}, outputFile=None, "
          "force=False, cut=False, verboseLevel=0)")
    print()
    print("Default extension:")
    print(extension)
    print()
    print("Default coefficients:")
    pprint(coefficients)

def eventdetectionwithfilegeneration(file, coefficients, verboseLevel=1, forceRun=False, CutTraces=False):
    # Create new class that contains all events
    AllEvents = NC.AllEvents()
    # Loop over all files in folder
    # Extract filename and generate filename to save located events
    if verboseLevel >= 1:
        print('analysing ' + file)
    directory = os.path.dirname(file) + os.sep + 'analysisfiles'
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)
    filename, file_extension = os.path.splitext(os.path.basename(file))
    savefilename = directory + os.sep + filename + 'data'
    shelfFile = shelve.open(savefilename)
    if ~forceRun:
        try: #Check if file can be loaded
            coefficientsloaded=shelfFile['coefficients']
            tEL = shelfFile['TranslocationEvents']
            # If coefficients before are not identical, analysis needs to run again
            if (coefficientsloaded == coefficients):
                if verboseLevel >= 1:
                    print('loaded from file')
            else:
                forceRun = True
        except:
            # Except if cannot be loaded, analysis needs to run
            forceRun = True
            pass

    if forceRun:
        # Extract list of events for this file
        tEL = eventdetection(file, coefficients, verboseLevel, CutTraces)
        if verboseLevel >= 1:
            print('Saved {} events'.format(len(tEL.events)))

        # Open savefile and save events for this file
        shelfFile['TranslocationEvents'] = tEL
        shelfFile['coefficients'] = coefficients
        if verboseLevel >= 1:
            print('saved to file')
    shelfFile.close()
    # Add events to the initially created class that contains all events
    AllEvents.AddEvent(tEL)
    return AllEvents

def batcheventdetection(folder, extension, coefficients, verboseLevel=1, forceRun=False, CutTraces=False):
    # Create new class that contains all events
    AllEvents = NC.AllEvents()

    print('found ' + str(len(glob.glob(os.path.join(folder, extension)))) + ' files.')

    # Loop over all files in folder
    for fullfilename in glob.glob(os.path.join(folder, extension)):

        # Extract filename and generate filename to save located events
        if verboseLevel >= 1:
            print('analysing '+fullfilename)
        filename, file_extension = os.path.splitext(os.path.basename(fullfilename))
        directory = os.path.dirname(fullfilename) + os.sep + 'analysisfiles'
        if not os.path.exists(directory):
            os.makedirs(directory,exist_ok=True)

        savefilename = directory + os.sep + filename + 'data'
        shelfFile = shelve.open(savefilename)
        if ~forceRun:
            try: #Check if file can be loaded
                coefficientsloaded=shelfFile['coefficients']
                tEL = shelfFile['TranslocationEvents']

                # If coefficients before are not identical, analysis needs to run again
                if (coefficientsloaded == coefficients):
                    if verboseLevel >= 1:
                        print('loaded from file')
                else:
                    forceRun = True
            except:
                # Except if cannot be loaded, analysis needs to run
                forceRun = True
                pass

        if forceRun:
            # Extract list of events for this file
            tEL = eventdetection(fullfilename, coefficients, verboseLevel, CutTraces)
            if verboseLevel >= 1:
                print('Saved {} events'.format(len(tEL.events)))

            # Open savefile and save events for this file
            shelfFile['TranslocationEvents'] = tEL
            shelfFile['coefficients'] = coefficients
            if verboseLevel >= 1:
                print('saved to file')

        shelfFile.close()

        # Add events to the initially created class that contains all events
        AllEvents.AddEvent(tEL)

    return AllEvents


def eventdetection(fullfilename, coefficients, verboseLevel=1, CutTraces=False, showFigures=False):
    """ 
    Function used to find the events of TranslocationEvents class in the raw data in file 'filename'. 
    It calls the function RecursiveLowPassFast to approximatively locate rough events in the data. 
    If a short TranslocationEvent object is detected its type attribute will be changed to 'Impulse' and the 
    meanTrace attribute will take the value of the minimal current value within the event. 
    
    Then the CUSUM function will be called  to build the CUSUM-fit and assign values to the different 
    attributes of the TranslocationEvent objects. 
    
    Depending on how the CUSUM was able to fit the trace inside and around the event, the type attribute of the TransocationEvent will be set to 'Real'
    (if the CUSUM fit went well) or 'Rough' (if the CUSUM was not able to fit the trace).
    
    Parameters
    ----------
        fullfilename : str
            Full path to data file.
        coefficients : dict
            Contains the parameters for the analysis.
        verboseLevel : bool, optional
            1 by default. It will print strings indicating the progress of the function in the console with different levels of depth.
        CutTraces : bool, optional
            False by default. If True, will cut the signal traces around the events to avoid having appended chunks detected as events.
        showFigures : bool , optional
            False by default. If True, it will display a simple figure with the shape of the signal.
    Returns
    -------
    list of TranslocationEvent
        All the events in the signal. 
    
    """

    if 'ChimeraLowPass' in coefficients:
        ChimeraLowPass = coefficients['ChimeraLowPass']
    else:
        ChimeraLowPass = None

    loadedData = LoadData.OpenFile(fullfilename, ChimeraLowPass, True, CutTraces)
    maxTime = coefficients['maxEventLength']
    minTime = coefficients['minEventLength']
    fitTime = coefficients['fitLength']
    IncludedBaseline = int(1e-2 * loadedData['samplerate'])
    delta = coefficients['delta']
    hbook = coefficients['hbook']
    dt = coefficients['dt']
    deltaRel = coefficients['deltaRel']
    samplerate = loadedData['samplerate']


    if verboseLevel >= 1:
        print('\n Recursive lowpass...', end='')

    # Call RecursiveLowPassFast to detect events in current trace
    start_time = timeit.default_timer()
    if 'i1Cut' in loadedData:
        events = []
        for cutTrace in loadedData['i1Cut']:
            events.extend(Functions.RecursiveLowPassFast(cutTrace, coefficients, loadedData['samplerate']))
    else:
        events = Functions.RecursiveLowPassFast(loadedData['i1'], coefficients, loadedData['samplerate'])
    if verboseLevel >= 1:
        print('done. Calculation took {}'.format(timeInSec.format_data(timeit.default_timer() - start_time)))
        print('Roughly detected events: {}'.format(len(events)))

    # Make a new class translocationEventList that contains all the found events in this file
    translocationEventList = NC.AllEvents()

    # Plot current file
    if showFigures:
        plt.figure(1)
        plt.cla()
        timeVals=np.linspace(0, len(loadedData['i1'])/loadedData['samplerate'], num=len(loadedData['i1']))
        plt.plot(timeVals,loadedData['i1'])

        plt.draw()
        plt.pause(0.001)

    if verboseLevel >= 1:
        print('Fine tuning...', end='')

    # Call RecursiveLowPassFast to detect events in current trace
    start_time = timeit.default_timer()
    cusumEvents = 0

    # Loop over all roughly detected events
    for event in events:
        beginEvent = event[0]
        endEvent = event[1]
        localBaseline = event[2]
        stdEvent = event[3]

        if 'blockLength' in loadedData:
            voltI = int(beginEvent / loadedData['blockLength'])
        else:
            voltI = int(0)

        # Check if the start of the event is later than the first few ms in IncludedBaseline
        if beginEvent > IncludedBaseline:

            # Extract trace and extract baseline
            traceforfitting = loadedData['i1'][int(beginEvent - dt):int(endEvent + dt)]

            # If event is shorter than fitTime
            if 3 < (endEvent-beginEvent) <= (fitTime * samplerate):
                # Try to fit and estimate real
                output = Functions.ImpulseFitting(traceforfitting, localBaseline, samplerate,
                                                  verbose=(verboseLevel >= 3))

            elif (endEvent-beginEvent) < (maxTime * samplerate):
                output = 2
            else:
                if verboseLevel >= 2:
                    print('Too long event {t:4.0f} ms\n'.format(t=1000*(endEvent-beginEvent)/samplerate))
                output = -1  # Don't include in rest of analysis

            # if output -1 fit failed, if output -2 differential failed
            if not isinstance(output, int) and len(output) > 1:
                (ps1, pe1, pe2, rsquared_event, Idrop) = output
                if verboseLevel >= 3:
                    print('ps1: {}, pe1: {}, pe2: {}, rsquared: {}'.format(ps1, pe1, pe2, rsquared_event))

                # length small enough and r squared big enough
                if (pe2-ps1)/samplerate < minTime and rsquared_event > 0.7:
                    newEvent = NC.TranslocationEvent(fullfilename, 'Impulse')
                    trace = traceforfitting[ps1:pe2]

                    ps1c = beginEvent - 25 + ps1
                    pe2c = beginEvent - 25 + pe2
                    traceBefore = loadedData['i1'][int(ps1c) - IncludedBaseline:int(ps1c)]
                    traceAfter = loadedData['i1'][int(pe2c):int(pe2c) + IncludedBaseline]
                    newEvent.SetEvent(trace, ps1c, localBaseline, samplerate, currentDrop=Idrop)
                    newEvent.SetCoefficients(coefficients, loadedData['v1'][voltI])
                    newEvent.SetBaselineTrace(traceBefore, traceAfter)
                    newEvent.eventLength = (pe1-ps1)/samplerate

                    if verboseLevel >= 2:
                        print('Fitted impulse of {t:1.3f} ms and {i:2.2f} nA and {r:1.2f} R-squared\n'.format(
                            t=newEvent.eventLength*1e3, i=newEvent.currentDrop*1e9, r=rsquared_event))
                elif (pe2-ps1)/samplerate < 3 * minTime and rsquared_event < 0.5:
                    if verboseLevel >= 2:
                        print('Bad fit {r:1.2f} R-squared, event ignored.\n'.format(r=rsquared_event))
                    output = -1  # Don't include in rest of analysis
                else:
                    if verboseLevel >= 3:
                        print('Too long event for impulse')
                    output = 2

            # Good enough for CUSUM fitting
            if isinstance(output, int) and abs(output) == 2:
                # CUSUM fit
                sigma = np.sqrt(stdEvent)

                # If deltaRel is set calculate delta based on relative value to baseline
                if deltaRel is not None:
                    delta = deltaRel * localBaseline
                    if verboseLevel >= 2:
                        print('Using relative delta for CUSUM fitting: {}'.format(Amp.format_data(delta)))

                h = hbook * delta / sigma
                (mc, kd, krmv) = Functions.CUSUM(traceforfitting, delta, h, verbose=(verboseLevel >= 3))
                krmv = [krmvVal + int(beginEvent) - dt + 1 for krmvVal in krmv]

                if len(krmv) > 1:
                    trace = loadedData['i1'][int(krmv[0]):int(krmv[-1])]
                    traceBefore = loadedData['i1'][int(krmv[0]) - IncludedBaseline:int(krmv[0])]
                    traceAfter = loadedData['i1'][int(krmv[-1]):int(krmv[-1]) + IncludedBaseline]
                    beginEvent = krmv[0]

                    # Adding CUSUM fitted event
                    newEvent = NC.TranslocationEvent(fullfilename, 'CUSUM')
                    newEvent.SetEvent(trace, beginEvent, localBaseline, samplerate)
                    newEvent.SetBaselineTrace(traceBefore, traceAfter)
                    newEvent.SetCoefficients(coefficients, loadedData['v1'][voltI])
                    newEvent.SetCUSUMVariables(mc, kd, krmv)

                    cusumEvents += 1
                    if verboseLevel >= 2:
                        print('Fitted CUSUM of {t:1.3f} ms and {i:2.2f} nA\n'.format(
                            t=newEvent.eventLengthCUSUM * 1e3, i=newEvent.currentDropCUSUM * 1e9))
                else:
                    trace = loadedData['i1'][int(beginEvent):int(endEvent)]
                    traceBefore = loadedData['i1'][int(beginEvent) - IncludedBaseline:int(beginEvent)]
                    traceAfter = loadedData['i1'][int(endEvent):int(endEvent) + IncludedBaseline]

                    # Adding roughly detected event
                    newEvent = NC.TranslocationEvent(fullfilename, 'Rough')
                    newEvent.SetEvent(trace, beginEvent, localBaseline, samplerate)
                    newEvent.SetBaselineTrace(traceBefore, traceAfter)
                    newEvent.SetCoefficients(coefficients, loadedData['v1'][voltI])

                    if verboseLevel >= 2:
                        print('CUSUM failed. Adding only roughly located event of {t:1.3f} ms and {i:2.2f} nA\n'.format(
                            t=len(trace) * 1e3/samplerate, i=(np.mean(np.append(traceBefore,traceAfter)) - np.mean(trace) )* 1e9))


            if output != -1:
                # Add event to TranslocationList
                translocationEventList.AddEvent(newEvent)

    if verboseLevel >= 1:
        print('done. Total fitting took {}'.format(timeInSec.format_data(timeit.default_timer() - start_time)))
        print('{} events fitted with CUSUM\n'.format(cusumEvents))

    #Plot events if True
    if showFigures:
        minVal = np.max(loadedData['i1'])#-1.05e-9 #np.min(loadedData['i1'])
        maxVal = np.max(loadedData['i1'])+0.05e-9#-1.05e-9 #np.min(loadedData['i1'])

        for i in range(len(events)):
            beginEvent=events[i][0]
            endEvent = events[i][1]
            if beginEvent>100 and (endEvent - beginEvent) >= (minTime * loadedData['samplerate']) and (endEvent - beginEvent) < (
                    coefficients['maxEventLength'] * loadedData['samplerate']):
                plt.plot([beginEvent/loadedData['samplerate'], beginEvent/loadedData['samplerate']], [minVal, maxVal], 'y-', lw=1)

        plt.show()
    return translocationEventList


def LoadEvents(loadname):
    # Check if loadname is a directory or not
    if os.path.isdir(loadname):
        # create standard filename for loading
        loadFile = os.path.join(loadname,os.path.basename(loadname)+'_Events')
    else:
        loadFile, file_extension = os.path.splitext(loadname)

    # Open file and extract TranslocationEvents
    shelfFile = shelve.open(loadFile)
    TranslocationEvents = shelfFile['TranslocationEvents']
    shelfFile.close()
    AllEvents = NC.AllEvents()
    AllEvents.AddEvent(TranslocationEvents)
    AllEvents.SetFolder(loadname)
    return AllEvents


def run(inputData, newExtension=None, newCoefficients={}, outputFile=None, force=False, cut=False, verboseLevel=0):
    """ 
    Function used to call all the other functions in the module 
    needed to find the events in raw nanopore experiment data.  
    
    Parameters
    -----------
    inputData : str
        Full path to data file.
    newExtension : str, optional
        None by default. NewExtension for input directory.
    newCoefficients : dict
        Contains the default parameters for the analysis.
    outputFile : str, optional
        None by default. Full path to output file.
    force : bool, optional
        False by default.
    cut : bool, optional
        False by default. False by default. If True, will cut the signal traces around the events to avoid having appended chunks detected as events.
    verboseLevel : int, optional
        0 by default. 0 by default. If higher, it will print various outputs during running

    Returns
    -------
    AllEvents object
        All the events.
    
    """
    
    if newExtension is None:
        newExtension = extension

    for key in newCoefficients:
        coefficients[key]=newCoefficients[key]

    if os.path.isdir(inputData):
        TranslocationEvents = batcheventdetection(inputData, newExtension, coefficients, verboseLevel, force, cut)
    else:
        TranslocationEvents = eventdetection(inputData,coefficients, verboseLevel, cut)

    #Check if list is empty
    if outputFile is not None and TranslocationEvents.events:
        if os.path.isdir(inputData):
            outputData = inputData + os.sep + 'Data' + os.sep + 'Data' + datetime.date.today().strftime("%Y%m%d")
        LoadData.SaveVariables(outputFile, TranslocationEvents=TranslocationEvents)
        return True
    else:
        return TranslocationEvents


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input directory or file')
    parser.add_argument('-e', '--ext', help='Extension for input directory')
    parser.add_argument('-o', '--output', help='Output file for saving')
    parser.add_argument('-c', '--coeff',
                        help='Coefficients for selecting events [-C filter E standarddev maxlength minlength',
                        nargs='+')
    parser.add_argument('-u', '--cut',
                        help='Cut Traces before detecting event, prevent detecting appended chunks as event',
                        action='store_true')
    parser.add_argument('-f', '--force', help='Force analysis to run (don''t load from file', action='store_true')

    args = parser.parse_args()
    inputData = args.input
    if inputData is None:
        inputData = askdirectory()

    outputData = args.output
    if outputData is None:
        if os.path.isdir(inputData):
            outputData = inputData + os.sep + 'Data' + os.sep + 'Data' + datetime.date.today().strftime("%Y%m%d")
        else:
            outputData = os.path.dirname(inputData) + os.sep + 'Data' + os.sep + 'Data' + datetime.date.today().strftime("%Y%m%d")

    if args.coeff is not None:
        if len(args.coeff) % 2 == 0:
            for i in range(0, len(args.coeff), 2):
                if i <= len(args.coeff):
                    coefficients[str(args.coeff[i])]=float(args.coeff[i+1])

    if args.ext is not None:
        extension = args.ext

    print('Loading from: ' + inputData)
    print('Saving to: ' + outputData)
    print('Coefficients are: ', end='')
    pprint(coefficients)



    if os.path.isdir(inputData):
        print('extension is: ' + extension +'\nStarting.... \n')
        TranslocationEvents = batcheventdetection(inputData, extension, newcoefficients=coefficients, force=args.force, cut=args.cut)
    else:
        print('Starting.... \n')
        TranslocationEvents = eventdetection(inputData,coefficients, args.cut)

    #Check if list is empty
    if TranslocationEvents.events:
        LoadData.SaveVariables(outputData, TranslocationEvents=TranslocationEvents)
