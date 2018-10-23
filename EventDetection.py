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


timeInSec= EngFormatter(unit='s', places=2)


def batcheventdetection(folder,extension,coefficients, forceRun=False, CutTraces=False):
    #Create new class that contains all events
    AllEvents=NC.AllEvents()

    #Loop over all files in folder
    for fullfilename in glob.glob(os.path.join(folder,extension)):

        #Extract filename and generate filename to save located events
        print('analysing '+fullfilename)
        filename, file_extension = os.path.splitext(os.path.basename(fullfilename))
        directory=os.path.dirname(fullfilename) + os.sep + 'analysisfiles'
        if not os.path.exists(directory):
            os.makedirs(directory,exist_ok=True)

        savefilename=directory + os.sep + filename + 'data'
        shelfFile = shelve.open(savefilename)
        if ~forceRun:
            try: #Check if file can be loaded
                coefficientsloaded=shelfFile['coefficients']
                tEL=shelfFile['translocationEvents']

                # If coefficients before are not identical, analysis needs to run again
                if (coefficientsloaded==coefficients):
                    print('loaded from file')
                else:
                    forceRun = True
            except:
                #Except if cannot be loaded, analysis needs to run
                forceRun=True
                pass

        if forceRun:
            #Extract list of events for this file
            tEL=eventdetection(fullfilename,coefficients,CutTraces)
            print('Saved {} events'.format(len(tEL.events)))

            #Open savefile and save events for this file
            shelfFile['translocationEvents']=tEL
            shelfFile['coefficients']=coefficients
            print('saved to file')

        shelfFile.close()

        #Add events to the initially created class that contains all events
        AllEvents.AddEvent(tEL)

    return AllEvents





def eventdetection(fullfilename, coefficients, CutTraces=False, showFigures = False):
    if 'ChimeraLowPass' in coefficients:
        ChimeraLowPass=coefficients['ChimeraLowPass']
    else:
        ChimeraLowPass=None
    loadedData=Functions.OpenFile(fullfilename, ChimeraLowPass, True, CutTraces)
    minTime=coefficients['minEventLength']
    IncludedBaseline = int(1e-2 * loadedData['samplerate'])
    delta=coefficients['delta']
    hbook=coefficients['hbook']


    print('Recursive lowpass...', end='')
    #Call RecursiveLowPassFast to detect events in current trace
    start_time = timeit.default_timer()
    if 'i1Cut' in loadedData:
        events=[]
        for cutTrace in loadedData['i1Cut']:
            events.extend(Functions.RecursiveLowPassFast(cutTrace,coefficients,loadedData['samplerate']))
    else:
        events=Functions.RecursiveLowPassFast(loadedData['i1'],coefficients,loadedData['samplerate'])
    print('done. Calculation took {}'.format(timeInSec.format_data(timeit.default_timer() - start_time)))
    print('nr events: {}'.format(len(events)))

    #CusumParameters = [delta sigma Normalize ImpulsionLimit IncludedBaseline TooMuchLevels IncludedBaselinePlots 1
    #                   SamplingFrequency];
    #[EventDatabase] = event_detection(RawSignal, CusumParameters, RoughEventLocations, lowpasscutoff);

    #Make a new class translocationEventList that contains all the found events in this file
    translocationEventList=NC.AllEvents()

    #Plot current file
    if showFigures:
        plt.figure(1)
        plt.cla()
        timeVals=np.linspace(0, len(loadedData['i1'])/loadedData['samplerate'], num=len(loadedData['i1']))
        plt.plot(timeVals,loadedData['i1'])

        plt.draw()
        plt.pause(0.001)

    print('CUSUM fitting...', end='')
    #Call RecursiveLowPassFast to detect events in current trace
    start_time = timeit.default_timer()
    #Loop over all detected events
    for event in events:
        beginEvent = event[0]
        endEvent = event[1]
        localBaseline = event[2]
        stdEvent = event[3]

        lengthevent=endEvent-beginEvent

        # Check if the start of the event is > IncludedBaseline then compare eventtime to minimal and maxima;
        if beginEvent > IncludedBaseline:
            # Extract trace and extract baseline
            Trace = loadedData['i1'][int(beginEvent):int(endEvent)]
            traceBefore = loadedData['i1'][int(beginEvent) - IncludedBaseline:int(beginEvent)]
            traceAfter = loadedData['i1'][int(endEvent):int(endEvent) + IncludedBaseline]

            if lengthevent<=(minTime*loadedData['samplerate']) and lengthevent>3:

                newEvent = NC.TranslocationEvent(fullfilename,'Impulse')

                # Add Trace, mean of the event,the samplerate, coefficients and baseline to the New Event class
                newEvent.SetEvent(Trace, beginEvent,localBaseline, loadedData['samplerate'])
                if 'blockLength' in loadedData:
                    voltI = int(beginEvent/loadedData['blockLength'])
                else:
                    voltI = int(0)
                newEvent.SetCoefficients(coefficients, loadedData['v1'][voltI])
                newEvent.SetBaselineTrace(traceBefore, traceAfter)

                # Add event to TranslocationList
                translocationEventList.AddEvent(newEvent)
            elif lengthevent>(minTime*loadedData['samplerate']) and lengthevent<(coefficients['eventlengthLimit']*loadedData['samplerate']):
                #Make new event class

                #CUSUM fit
                sigma = np.sqrt(stdEvent)
                h = hbook * delta / sigma
                (mc, kd, krmv)= Functions.CUSUM(loadedData['i1'][int(beginEvent) - IncludedBaseline: int(endEvent) + IncludedBaseline], delta, h)
                krmv = [krmvVal + int(beginEvent) - IncludedBaseline + 1 for krmvVal in krmv]
                #Add Trace, mean of the event,the samplerate, coefficients and baseline to the New Event class
                if len(krmv)>2:
                    newEvent = NC.TranslocationEvent(fullfilename, 'Real')
                    Trace=loadedData['i1'][int(krmv[0]):int(krmv[-1])]
                    traceBefore = loadedData['i1'][int(krmv[0]) - IncludedBaseline:int(krmv[0])]
                    traceAfter = loadedData['i1'][int(krmv[-1]):int(krmv[-1]) + IncludedBaseline]
                    beginEvent=krmv[0]
                else:
                    newEvent = NC.TranslocationEvent(fullfilename, 'Rough')

                newEvent.SetEvent(Trace,beginEvent,localBaseline,loadedData['samplerate'])
                newEvent.SetBaselineTrace(traceBefore, traceAfter)
                newEvent.SetCUSUMVariables(mc, kd, krmv)

                if 'blockLength' in loadedData:
                    voltI = int(beginEvent/loadedData['blockLength'])
                else:
                    voltI = int(0)
                newEvent.SetCoefficients(coefficients,loadedData['v1'][voltI])

                #Add event to TranslocationList
                translocationEventList.AddEvent(newEvent)

    print('done. Calculation took {}'.format(timeInSec.format_data(timeit.default_timer() - start_time)))

    #Plot events if True
    if showFigures:
        minVal=np.max(loadedData['i1'])#-1.05e-9 #np.min(loadedData['i1'])
        maxVal=np.max(loadedData['i1'])+0.05e-9#-1.05e-9 #np.min(loadedData['i1'])

        for i in range(len(events)):
            beginEvent=events[i][0]
            endEvent = events[i][1]
            if beginEvent>100 and (endEvent - beginEvent) >= (minTime * loadedData['samplerate']) and (endEvent - beginEvent) < (
                    coefficients['eventlengthLimit'] * loadedData['samplerate']):
                plt.plot([beginEvent/loadedData['samplerate'], beginEvent/loadedData['samplerate']], [minVal, maxVal], 'y-', lw=1)

        #plt.draw()
        plt.show()
        #plt.pause(1)
    return translocationEventList




def LoadEvents(loadname):
    #Check if loadname is a directory or not
    if os.path.isdir(loadname):
        #create standard filename for loading
        loadFile = os.path.join(loadname,os.path.basename(loadname)+'_Events')
    else:
        loadFile, file_extension = os.path.splitext(loadname)

    #Open file and extract TranslocationEvents
    shelfFile=shelve.open(loadFile)
    TranslocationEvents=shelfFile['TranslocationEvents']
    shelfFile.close()
    AllEvents=NC.AllEvents()
    AllEvents.AddEvent(TranslocationEvents)
    AllEvents.SetFolder(loadname)
    return AllEvents

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input directory or file')
    parser.add_argument('-e', '--ext', help='Extension for input directory')
    parser.add_argument('-o', '--output', help='Output file for saving')
    parser.add_argument('-c', '--coeff', help='Coefficients for selecting events [-C filter E standarddev maxlength minlength', nargs = '+')
    parser.add_argument('-u', '--cut', help='Cut Traces before detecting event, prevent detecting appended chunks as event', action='store_true')
    parser.add_argument('-f', '--force', help='Force analysis to run (don''t load from file', action='store_true')

    args = parser.parse_args()
    inputData=args.input
    if inputData==None:
        inputData=askdirectory()

    outputData=args.output
    if outputData==None:
        if os.path.isdir(inputData):
            outputData = inputData + os.sep + 'Data' + os.sep + 'Data' + datetime.date.today().strftime("%Y%m%d")
        else:
            outputData=os.path.dirname(inputData) + os.sep + 'Data' + os.sep + 'Data' + datetime.date.today().strftime("%Y%m%d")


    coefficients = {'a': 0.999, 'E': 0, 'S': 5, 'eventlengthLimit': 200e-3, 'minEventLength': 500e-6, 'hbook':1,'delta':0.2e-9,'ChimeraLowPass':10e3}
    if args.coeff is not None:
        if len(args.coeff) % 2 == 0:
            for i in range(0, len(args.coeff), 2):
                if i <= len(args.coeff):
                    coefficients[str(args.coeff[i])]=float(args.coeff[i+1])

    extension=args.ext
    if extension==None:
        extension='*.log'


    print('Loading from: ' + inputData)
    print('Saving to: ' + outputData)
    print('Coefficients are: ', end='')
    pprint(coefficients)



    if os.path.isdir(inputData):
        print('extension is: ' + extension +'\nStarting.... \n')
        TranslocationEvents=batcheventdetection(inputData, extension, coefficients, args.force, args.cut)
    else:
        print('Starting.... \n')
        TranslocationEvents=eventdetection(inputData,coefficients, args.cut)

    #Check if list is empty
    if TranslocationEvents.events:
        Functions.SaveVariables(outputData,TranslocationEvents=TranslocationEvents)