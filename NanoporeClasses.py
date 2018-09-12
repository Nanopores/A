import numpy as np
from pprint import pprint
import shelve
import os
import math
from Plotting.NanoporePlots import PlotCurrentTrace,PlotCurrentTraceBaseline
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter

Amp = EngFormatter(unit='A', places=2)
Time = EngFormatter(unit='s', places=2)
Volt = EngFormatter(unit='V', places=2)


class TranslocationEvent:
    def __init__(self, filename,type='roughEvent'):
        self.filename = filename
        self.type=type

    def SetEvent(self,eventTrace,baseline,samplerate):
        self.eventTrace=eventTrace
        self.baseline=baseline
        self.samplerate=samplerate

        self.meanTrace=np.mean(eventTrace)
        self.minTrace = np.min(eventTrace)
        self.lengthEvents=len(eventTrace)/samplerate

        if self.type=='Real':
            self.currentDrop=baseline-self.meanTrace
        else:
            self.currentDrop = baseline - self.minTrace

    def SetCoefficients(self,coefficients,voltage):
        self.coefficients=coefficients
        self.voltage=voltage

    def SetBaselineTrace(self, before,after):
        self.before=before
        self.after=after

    def SetCUSUMVariables(self, segmentedSignal, kd, changeTimes):
        self.segmentedSignal=segmentedSignal
        self.changeTimes=changeTimes

        self.kd=kd
        # self.krmv=krmv
        # self.mc=mc

    def PlotEvent(self):
        part1=np.append(self.before,self.eventTrace)

        fn=filename_w_ext = os.path.basename(self.filename)

        plotTitle = fn + '\n' + 'Event length: {}\nCurrent drop: {}'.format(Time.format_data(self.lengthEvents), Amp.format_data(self.currentDrop))
        #PlotCurrentTraceBaseline(self.before, self.eventTrace, self.after, self.samplerate, titleplot)

        timeVals1 = np.linspace(0, len(self.before) / self.samplerate, num=len(self.before))
        timeVals2 = np.linspace(0 + max(timeVals1), len(self.eventTrace) / self.samplerate + max(timeVals1),
                                num=len(self.eventTrace))
        timeVals3 = np.linspace(0 + max(timeVals2), len(self.after) / self.samplerate + max(timeVals2), num=len(self.after))

        # plt.figure(figsize=(10, 6))
        fig, ax = plt.subplots(figsize=(10, 6))

        ax.plot(timeVals1, self.before, color='tomato')
        ax.plot(timeVals2, self.eventTrace, color='mediumslateblue')
        ax.plot(timeVals3, self.after, color='tomato')

        beforeBaseline=np.full(len(self.before), self.baseline)
        ax.plot(timeVals1,beforeBaseline, '--', color='tomato')
        afterBaseline = np.full(len(self.after), self.baseline)
        ax.plot(timeVals3,afterBaseline, '--', color='tomato')

        meanTrace = np.full(len(self.eventTrace), self.meanTrace)
        ax.plot(timeVals2,meanTrace, '--', color='mediumslateblue')

        ax.set_xlabel('time (s)')
        ax.set_ylabel('current (A)')
        ax.xaxis.set_major_formatter(Time)
        ax.yaxis.set_major_formatter(Amp)

        if plotTitle:
            plt.title(plotTitle)

        plt.show()

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

    def GetAllIdropsNorm(self):
        currentDrops = [event.currentDrop/ event.baseline for event in self.events]
        return currentDrops

    def SaveEvents(self,savename):
        if os.path.isdir(savename):
            savefile=os.path.join(savename,os.path.basename(savename)+'_Events')
        else:
            if os.path.isfile(savename + '.dat'):
                raise IOError('File ' + savename + '.dat already exists.')
            else:
                savefile = savename

        #Check if directory exists
        directory = os.path.dirname(savefile)
        if not os.path.exists(directory):
            os.makedirs(directory)

        shelfFile=shelve.open(savefile)
        shelfFile['TranslocationEvents']=self.events
        shelfFile.close()
        print('saved as: ' + savefile + '.dat')
        # if os.path.exists(savefile + '.bak'):
        #     os.remove(savefile + '.bak')
        # if os.path.exists(savefile + '.dir'):
        #     os.remove(savefile + '.dir')

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

    def PlotIEvent(self,i):
        event=self.events[i]
        event.Plotevent()
