import NanoporeClasses as NC
import shelve
import os
import glob
from time import sleep
import Functions
from pprint import pprint
import matplotlib.pyplot as plt
import numpy as np



def batcheventdetection(folder,extension,coefficients):
    AllEvents=NC.AllEvents()
    for fullfilename in glob.glob(os.path.join(folder,extension)):
        print('analysing '+fullfilename)
        filename, file_extension = os.path.splitext(fullfilename)
        savefilename=filename+'data'
        tEL=eventdetection(fullfilename,coefficients)
        shelfFile=shelve.open(savefilename)
        shelfFile['translocationEvents']=tEL
        #tEL=shelfFile['translocationEvents']
        pprint('found {} events'.format(len(tEL.events)))
        shelfFile.close()
        AllEvents.AddEvent(tEL)

    return AllEvents





def eventdetection(fullfilename,coefficients,showFigures = False):
    loadedData=Functions.OpenFile(fullfilename)
    #pprint(loadedData)
    minTime=coefficients['minEventLength']

    print('eventdetection...', end='')
    events=Functions.RecursiveLowPassFast(loadedData['i1'],coefficients,loadedData['samplerate'])
    # print('nr events: {}'.format(len(events)))

    translocationEventList=NC.AllEvents()

    if showFigures:
        plt.figure(1)
        plt.cla()
        timeVals=np.linspace(0, len(loadedData['i1'])/loadedData['samplerate'], num=len(loadedData['i1']))
        plt.plot(timeVals,loadedData['i1'])

        plt.draw()
        plt.pause(0.001)

    # for i in range(length(events)):

    #print('min ticks:{}'.format(minTime*loadedData['samplerate']))
    for i in range((len(events))):
        beginEvent=events[i][0]
        endEvent=events[i][1]
        meanEvent=events[i][2]
        stdEvent=events[i][3]

        if beginEvent>100 and (endEvent-beginEvent)>=(minTime*loadedData['samplerate']) and (endEvent-beginEvent)<(coefficients['eventlengthLimit']*loadedData['samplerate']):
            newEvent=NC.TranslocationEvent(fullfilename)
            Trace=loadedData['i1'][int(beginEvent):int(endEvent)]
            traceBefore=loadedData['i1'][int(beginEvent)-100:int(beginEvent)+1]
            traceAfter=loadedData['i1'][int(endEvent)-1:int(endEvent)+100]
            newEvent.SetEvent(Trace,meanEvent,loadedData['samplerate'])
            newEvent.SetCoefficients(coefficients)
            newEvent.SetBaselineTrace(traceBefore,traceAfter)
            translocationEventList.AddEvent(newEvent)
            print(newEvent.currentDrop)


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
    if os.path.isdir(loadname):
        loadFile = os.path.join(loadname,os.path.basename(loadname)+'_Events')
    else:
        loadFile, file_extension = os.path.splitext(loadname)

    shelfFile=shelve.open(loadFile)
    TranslocationEvents=shelfFile['TranslocationEvents']
    shelfFile.close()
    AllEvents=NC.AllEvents()
    AllEvents.AddEvent(TranslocationEvents)
    return AllEvents
