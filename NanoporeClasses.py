import numpy as np
from pprint import pprint
import shelve
import os
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter

Amp = EngFormatter(unit='A', places=2)
Time = EngFormatter(unit='s', places=2)
Volt = EngFormatter(unit='V', places=2)
Cond = EngFormatter(unit='S', places=2)


class TranslocationEvent:
    """
    Class used to represent an event.

    Attributes
    ----------
    filename : str
        directory of the source file from which the event was extracted
    type: str
        classification of the event 'CUSUM' (CUSUM-fitted), 'Rough' (non-CUSUM fitted), 'Impulse' (short events)
    evenTrace: lst(int)
        list of data points within the event
    baseline:float
        mean current value in the basline around the event
    samplerate:float
        sampling frequency
    beginEvent: int
        start coordinate of the event in the signal
    endEvent: int
        end coordinate of the event in the signal
    meanTrace: float
        mean current in the event
    minTrace: float
        minimum current in the event
    eventLength: float
        lenght of the event in seconds
    currentDrop: float
        current drop defined as the difference of the baseline and meanTrace or minTrace dependeng on the event type
    before: lst(float)
        list of data points before the event, used for plotting
    after: lst(float)
        list of data points after the event, used for plotting
    changeTimes: lst(int)
        list of coordinates where the current level changes in the event, only for multilevel events
    kd: ****
    segmentedSignal: lst(float)
        CUSUM fit
    beginEventCUSUM:
        more precise start coordinate of the event detected with the CUSUM algorithm
    currentDropCUSUM:
        more precise start coordinate of the event detected with the CUSUM algorithm
    coefficients:
        coefficients used in the CUSUM algorithm
    voltage:
        voltage applied across the nanopore during the data recording

    """

    def __init__(self, filename, type='roughEvent'):
        self.filename = filename
        self.type = type

    def SetEvent(self, eventTrace, beginEvent, baseline, samplerate, currentDrop=None):
        self.eventTrace = eventTrace
        self.baseline = baseline
        self.samplerate = samplerate
        self.beginEvent = beginEvent
        self.endEvent = beginEvent+len(eventTrace)

        self.meanTrace = np.mean(eventTrace)
        self.minTrace = np.min(eventTrace)
        self.eventLength = len(eventTrace)/samplerate

        if currentDrop is None:
            self.currentDrop = baseline - self.meanTrace
        else:
            self.currentDrop = currentDrop

    def SetCoefficients(self,coefficients,voltage):
        self.coefficients = coefficients
        self.voltage = voltage

    def SetBaselineTrace(self, before,after):
        self.before = before
        self.after = after

        self.baseline=np.mean(np.append(before,after))


    def SetCUSUMVariables(self, segmentedSignal, kd, changeTimes):
        self.changeTimes = changeTimes
        self.kd = kd
        self.segmentedSignal = segmentedSignal
        self.changeTimes = changeTimes
        if len(changeTimes) >= 1:
            self.beginEventCUSUM = changeTimes[0]
            self.currentDropCUSUM = max(segmentedSignal)-min(segmentedSignal)

        if len(changeTimes) >= 2:
            self.endEventCUSUM = changeTimes[-1]
            self.eventLengthCUSUM = (changeTimes[-1]-changeTimes[0])/self.samplerate
            if hasattr(self,'before') and hasattr(self,'after') and hasattr(self,'eventTrace'):
                self.mcbefore=np.mean(self.before)*np.ones(len(self.before))
                self.mcafter = np.mean(self.after) * np.ones(len(self.after))
                self.mctrace=np.array([])
                for ii in range(1,len(changeTimes)):
                    self.mctrace=np.append(self.mctrace,np.mean(self.eventTrace[changeTimes[ii-1]-changeTimes[0]:changeTimes[ii]-changeTimes[0]])*np.ones(changeTimes[ii]-changeTimes[ii-1]))



class AllEvents:
    """"
    Class used to represent all the events in a nanopore experiment output as a list of events.

    Attributes
    ----------

    events : lst(events)
        list of events containing all the events detected by the eventdetection function

    """

    def __init__(self):
        self.events=[]

    def AddEvent(self, translocationEvent):
        if isinstance(translocationEvent,AllEvents):
            if len(self.events) == 0:
                self.events = translocationEvent.events
            else:
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

    def GetEventTypes(self,eventType):
        events = []
        events = [event for event in self.events if str.upper(event.type) == str.upper(eventType)]
        return events

    def GetEventsforVoltages(self,voltage):
        events = [event for event in self.events if event.voltage == voltage]
        return events

    def GetAllVoltages(self):
        voltages = [event.voltage for event in self.events]
        voltages = set(voltages)
        return voltages

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
        event = self.events[i]
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
