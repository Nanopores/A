import NanoporeClasses as NC
import shelve
import os
import glob
from time import sleep
import Functions
from pprint import pprint
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import datetime
from tkinter.filedialog import askopenfilenames,askdirectory
import timeit



def batcheventdetection(folder,extension,coefficients):
    #Create new class that contains all events
    AllEvents=NC.AllEvents()

    #Loop over all files in folder
    for fullfilename in glob.glob(os.path.join(folder,extension)):

        #Extract filename and generate filename to save located events
        print('analysing '+fullfilename)
        filename, file_extension = os.path.splitext(fullfilename)
        savefilename=filename+'data'
        shelfFile = shelve.open(savefilename)
        try:
            coefficientsloaded=shelfFile['coefficients']
            tEL=shelfFile['translocationEvents']
            assert(coefficientsloaded==coefficients)
            print('loaded from file')
        except KeyError or AssertionError:
            #Extract list of events for this file
            tEL=eventdetection(fullfilename,coefficients)
            pprint('found {} events'.format(len(tEL.events)))

            #Open savefile and save events for this file
            shelfFile['translocationEvents']=tEL
            shelfFile['coefficients']=coefficients
            print('saved to file')
            pass
        shelfFile.close()

        #Add events to the initially created class that contains all events
        AllEvents.AddEvent(tEL)

    return AllEvents





def eventdetection(fullfilename,coefficients,showFigures = False):
    loadedData=Functions.OpenFile(fullfilename,None,True)
    minTime=coefficients['minEventLength']

    print('eventdetection...', end='')
    #Call RecursiveLowPassFast to detect events in current trace
    start_time = timeit.default_timer()
    events=Functions.RecursiveLowPassFast(loadedData['i1'],coefficients,loadedData['samplerate'])
    print('done. Calculation took ' + str(timeit.default_timer() - start_time) + 's')
    print('nr events: {}'.format(len(events)))

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

    # for i in range(length(events)):

    #print('min ticks:{}'.format(minTime*loadedData['samplerate']))

    #Loop over all detected events
    for i in range((len(events))):
        beginEvent=events[i][0]
        endEvent=events[i][1]
        meanEvent=events[i][2]
        stdEvent=events[i][3]

        # Check if the start of the event is > 100 and eventtime is more then the minimal and less than the maxima;
        if beginEvent>100 and (endEvent-beginEvent)>=(minTime*loadedData['samplerate']) and (endEvent-beginEvent)<(coefficients['eventlengthLimit']*loadedData['samplerate']):
            #Make new event class
            newEvent=NC.TranslocationEvent(fullfilename)
            #Extract trace and extract baseline
            Trace=loadedData['i1'][int(beginEvent):int(endEvent)]
            traceBefore=loadedData['i1'][int(beginEvent)-100:int(beginEvent)+1]
            traceAfter=loadedData['i1'][int(endEvent)-1:int(endEvent)+100]

            #Add Trace, mean of the event,the samplerate, coefficients and baseline to the New Event class
            newEvent.SetEvent(Trace,meanEvent,loadedData['samplerate'])
            newEvent.SetCoefficients(coefficients)
            newEvent.SetBaselineTrace(traceBefore,traceAfter)

            #Add event to TranslocationList
            translocationEventList.AddEvent(newEvent)

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
    return AllEvents

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input directory or file')
    parser.add_argument('-e', '--ext', help='extension for input directory')
    parser.add_argument('-o', '--output', help='output file for saving')
    parser.add_argument('-c', '--coeff', help='Coefficients for selecting events [-C filter E standarddev maxlength minlength', type = float, nargs = 5)

    args = parser.parse_args()
    inputData=args.input
    if inputData==None:
        inputData=askdirectory()

    outputData=args.output
    if outputData==None:
        outputData=os.path.dirname(inputData)+'_data_'+datetime.date.today().strftime("%Y%m%d")

    if args.coeff==None:
        coefficients= {'a': 0.99, 'E': 0, 'S': 5, 'eventlengthLimit': 0.5,'minEventLength': 100e-6}
    else:
        coefficients = {'a': args.coeff[0] , 'E': args.coeff[1], 'S': args.coeff[2], 'eventlengthLimit': args.coeff[3], 'minEventLength': args.coeff[4]}

    extension=args.ext
    if extension==None:
        extension='*.log'

    print('Loading from: ' + inputData)
    print('Saving to: ' + outputData)
    print('Coefficients are: ', end='')
    pprint(coefficients)



    if os.path.isdir(inputData):
        print('extension is: ' + extension +'\nStarting.... \n')
        TranslocationEvents=batcheventdetection(inputData, extension, coefficients)
    else:
        print('Starting.... \n')
        TranslocationEvents=eventdetection(inputData,coefficients)

    #Check if list is empty
    if TranslocationEvents.events:
        TranslocationEvents.SaveEvents(outputData)