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

    def SetEvent(self,eventTrace,beginEvent,baseline,samplerate):
        self.eventTrace=eventTrace
        self.baseline=baseline
        self.samplerate=samplerate
        self.beginEvent=beginEvent
        self.endEvent=beginEvent+len(eventTrace)

        self.meanTrace=np.mean(eventTrace)
        self.minTrace = np.min(eventTrace)
        self.eventLength=len(eventTrace)/samplerate

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
        self.segmentedSignal=segmentedSignal
        self.changeTimes=changeTimes
        self.beginEventCUSUM=changeTimes[0]
        self.endEventCUSUM = changeTimes[-1]
        self.eventLengthCUSUM=(changeTimes[-1]-changeTimes[0])/self.samplerate

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

        ax.plot(np.append(timeVals1,timeVals2[0]), np.append(self.before,self.eventTrace[0]), color='tomato')
        ax.plot(timeVals2, self.eventTrace, color='mediumslateblue')
        ax.plot(np.append(timeVals2[-1],timeVals3), np.append(self.eventTrace[-1],self.after), color='tomato')

        beforeBaseline=np.full(len(self.before), self.baseline)
        ax.plot(timeVals1,beforeBaseline, '--', color='tomato')
        afterBaseline = np.full(len(self.after), self.baseline)
        ax.plot(timeVals3,afterBaseline, '--', color='tomato')

        meanTrace = np.full(len(self.eventTrace), self.baseline-self.currentDrop)
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

        self.savefile=savefile

        #Check if directory exists
        directory = os.path.dirname(savefile)
        if not os.path.exists(directory):
            os.makedirs(directory)

        shelfFile=shelve.open(savefile)
        shelfFile['TranslocationEvents']=self.events
        shelfFile.close()
        print('saved as: ' + savefile + '.dat')


    def SetFolder(self,loadname):
        self.savefile=loadname

    def GetEventsMinCondition(self,minCurrent=-math.inf,maxCurrent=math.inf,minLength=0,maxLength=math.inf):
        minCurrent = -math.inf if not minCurrent else minCurrent
        maxCurrent = math.inf if not maxCurrent else maxCurrent
        minLength = 0 if not minLength else minLength
        maxLength = math.inf if not maxLength else maxLength

        newEvents=AllEvents()
        for event in self.events:
            if minCurrent<event.currentDrop<maxCurrent and minLength<event.lengthEvents<maxLength:
                newEvents.AddEvent(event)
        newEvents.SetFolder(self.savefile)
        print('selected {} events from {}'.format(len(newEvents.events), len(self.events)))
        return newEvents

    def PlotAllEvents(self):
        for event in self.events:
            event.PlotEvent()

    def PlotIEvent(self,i):
        event=self.events[i]
        event.Plotevent()

    def PlotHistogram(self):
        appendTrace=[]
        for event in self.events:
            appendTrace.extend(event.eventTrace)

        fig, ax = plt.subplots(figsize=(6, 10))

        weights = np.ones_like(appendTrace) / float(len(appendTrace))
        ax.hist(appendTrace,weights=weights,bins=40,orientation='horizontal')
        ax.yaxis.set_major_formatter(Amp)

        savefile=self.savefile
        # Check if directory exists
        directory = os.path.dirname(savefile)
        fig.savefig(directory + os.sep + 'Histogram.pdf', transparent=True)
        plt.show()


    def PlotI_tau(self,  normalized=False):

        current = [event.currentDrop for event in self.events]
        tau=[event.eventLength for event in self.events]

        # definitions for the axes
        left, width = 0.15, 0.55
        bottom, height = 0.1, 0.6
        left_h = left + width + 0.015
        bottom_h = bottom + height + 0.015

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.2]
        rect_histy = [left_h, bottom, 0.2, height]

        # start with a rectangular Figure
        fig = plt.figure(1, figsize=(10, 8))

        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx, sharex=axScatter)
        axHisty = plt.axes(rect_histy, sharey=axScatter)

        # the scatter plot:
        axScatter.scatter(tau, current, color='tomato', marker='o', s=30, linewidths=0.1, edgecolors='red') # added some stuff here to improve aesthetics

        # now determine nice limits by hand:
        extra = 0.1  # 0.1 = 10%
        tauRange = np.max(tau) - np.min(tau)
        currentClean = [x for x in current if str(x) != 'nan']
        currentRange = np.max(currentClean) - np.min(currentClean)
        #
        # tauLim = (int(tauMax/binwidth) + 1) * binwidth
        # currentLim = (int(currentMax/binwidth) + 1) * binwidth

        axScatter.set_xlim((np.min(tau) - extra * tauRange, np.max(tau) + extra * tauRange))
        axScatter.set_ylim((np.min(currentClean) - extra * currentRange, np.max(currentClean) + extra * currentRange))

        # tauBins = np.arange(-tauLim, tauLim + binwidth, binwidth)
        # currentBins = np.arange(-currentLim, currentLim + binwidth, binwidth)
        axHistx.hist(tau, bins=50, color='tomato')
        plt.setp(axHistx.get_xticklabels(), visible=False)
        axHisty.hist(currentClean, bins=50, orientation='horizontal', color='tomato')   #increased binning here to avoid overlapping bars in the histogram, added color option
        plt.setp(axHisty.get_yticklabels(), visible=False)

        axScatter.set_xlabel('Event length (s)')
        axScatter.xaxis.set_major_formatter(Time)
        if normalized:
            axScatter.set_ylabel('current drop (normalized)')
        else:
            axScatter.set_ylabel('current drop (A)')
            axScatter.yaxis.set_major_formatter(Amp)


        savefile=self.savefile
        # Check if directory exists
        directory = os.path.dirname(savefile)
        fig.savefig(directory + os.sep + 'PlotITau.pdf', transparent=True)

        plt.show()

