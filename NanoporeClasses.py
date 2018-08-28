import numpy as np
from pprint import pprint
import shelve
import os
import math
from Plotting.NanoporePlots import PlotCurrentTrace,PlotCurrentTraceBaseline


class TranslocationEvent:
    def __init__(self, filename):
        self.filename = filename

    def SetEvent(self,eventTrace,baseline,samplerate):
        self.eventTrace=eventTrace
        self.baseline=baseline
        self.samplerate=samplerate

        self.meanTrace=np.mean(eventTrace)
        self.lengthEvents=len(eventTrace)/samplerate
        self.currentDrop=baseline-np.mean(eventTrace)

    def SetCoefficients(self,coefficients):
        self.coefficients=coefficients

    def SetBaselineTrace(self, before,after):
        self.before=before
        self.after=after

    def PlotEvent(self):
        part1=np.append(self.before,self.eventTrace)
        currentTrace=np.append(part1,self.after)
        #PlotCurrentTrace(currentTrace,self.samplerate)
        PlotCurrentTraceBaseline(self.before,self.eventTrace,self.after,self.samplerate)



class AllEvents:
    def __init__(self):
        self.events=[]

    def AddEvent(self, translocationEvent):
        if isinstance(translocationEvent,AllEvents):
            self.events.extend(translocationEvent.events)
        elif isinstance(translocationEvent,list):
            self.events.extend(translocationEvent)
        else:
            self.events.append(translocationEvent)

    def GetAllLengths(self):
        Lengths=[event.lengthEvents for event in self.events]
        return Lengths

    def GetAllIdrops(self):
        currentDrops=[event.currentDrop for event in self.events]
        return currentDrops

    def SaveEvents(self,folder):
        savefile=os.path.join(folder,os.path.basename(folder)+'_Events')
        shelfFile=shelve.open(savefile)
        shelfFile['TranslocationEvents']=self.events
        shelfFile.close()
        print('saved as: ' + savefile)

    def GetEventsMinCondition(self,minCurrent=-math.inf,maxCurrent=math.inf,minLength=0,maxLength=math.inf):
        minCurrent = -math.inf if not minCurrent else minCurrent
        maxCurrent = math.inf if not maxCurrent else maxCurrent
        minLength = 0 if not minLength else minLength
        maxLength = math.inf if not maxLength else maxLength

        newEvents=AllEvents()
        for event in self.events:
            if minCurrent<event.currentDrop<maxCurrent and minLength<event.lengthEvents<maxLength:
                newEvents.AddEvent(event)
        print('selected {} events from {}'.format(len(newEvents.events), len(self.events)))
        return newEvents

    def PlotAllEvents(self):
        for event in self.events:
            event.PlotEvent()
