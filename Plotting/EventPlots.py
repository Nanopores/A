import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import EngFormatter
from matplotlib.widgets import CheckButtons, Button
import argparse
from tkinter.filedialog import askopenfilenames,askdirectory
import shelve
import os
from pprint import pprint
import Functions
import LoadData


Amp = EngFormatter(unit='A', places=2)
Time = EngFormatter(unit='s', places=2)
Volt = EngFormatter(unit='V', places=2)
Cond = EngFormatter(unit='S', places=2)


def PlotG_tau(events, savefile = None, showCurrentInstead=False, normalized=False,showCUSUM=True):
    #Clean up events
    events = [event for event in events if event.voltage is not 0]

    #categorize events in three types
    #CUSUM fitted events
    if any(hasattr(event, 'changeTimes') and len(event.changeTimes)>2 for event in events):
        CUSUMIndices, CUSUMEdevents = zip(*[(ind, event) for ind, event in enumerate(events) if
                                        hasattr(event, 'changeTimes') and len(event.changeTimes) > 2])
    else:
        CUSUMIndices = CUSUMEdevents = []

    #Non-fitted events
    if any(event.type =='Rough' for event in events):
        nonFittedIndices, nonFittedEvents = zip(*[(ind, event) for ind, event in enumerate(events) if
                                            event.type == 'Rough'])
    else:
        nonFittedIndices = nonFittedEvents = []

    #Impulse events
    if any(event.type == 'Impulse' for event in events):
        impulseIndices, impulseEvents = zip(*[(ind, event) for ind, event in enumerate(events) if
                                              event.type == 'Impulse'])
    else:
        impulseIndices = impulseEvents = []

    catEvents=(CUSUMEdevents, nonFittedEvents, impulseEvents)
    catIndices=(CUSUMIndices,nonFittedIndices,impulseIndices)

    #Extract y and tau out of events
    def extractytau(filteredevents):
        if showCurrentInstead:
            yVals = [event.currentDrop for event in filteredevents]
            tau = [event.eventLength for event in filteredevents]
        else:
            yVals = [event.currentDrop / event.voltage for event in filteredevents if event.voltage > 0]
            tau = [event.eventLength for event in filteredevents if event.voltage > 0]
        return tau,yVals

    #Save figure
    def SavePlot(event):
        # Check if directory exists
        directory = os.path.dirname(savefile)
        if showCurrentInstead:
            fig.savefig(directory + os.sep + 'PlotITau.pdf', transparent=True)
        else:
            fig.savefig(directory + os.sep + 'PlotGTau.pdf', transparent=True)

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

    #define axes and link histogram to scatterplot
    axScatter = plt.axes(rect_scatter)

    axHistx = plt.axes(rect_histx, sharex=axScatter)
    axHisty = plt.axes(rect_histy, sharey=axScatter)

    #Checkboxes to turn on or off events
    rax = plt.axes([0.75, 0.73, 0.14, 0.15])
    visBool=[True,True,True]
    labelsCheckBox=('CUSUM-fitted', 'Not fitted', 'impulse')
    check = CheckButtons(rax,labelsCheckBox , visBool)

    bax = plt.axes([0.77, 0.9, 0.1, 0.03])
    bnext = Button(bax, 'Save figure')
    bnext.on_clicked(SavePlot)

    #Show labels
    def setlabels():
        plt.setp(axHistx.get_xticklabels(), visible=False)
        plt.setp(axHisty.get_yticklabels(), visible=False)

        axScatter.set_xlabel('Event length (s)')
        axScatter.xaxis.set_major_formatter(Time)
        if normalized:
            axScatter.set_ylabel('current drop (normalized)')
        else:
            if showCurrentInstead:
                axScatter.set_ylabel('current drop (A)')
                axScatter.yaxis.set_major_formatter(Amp)
            else:
                axScatter.set_ylabel('Conductance drop (G)')
                axScatter.yaxis.set_major_formatter(Cond)



    # Determine nice limits by hand:
    def limits(tau, yVals, axScatter):
        extra = 0.1  # 0.1 = 10%
        tauRange = np.max(tau) - np.min(tau)
        yClean = [y for y in yVals if str(y) != 'nan']
        yRange = np.max(yClean) - np.min(yClean)

        axScatter.set_xlim((np.min(tau) - extra * tauRange, np.max(tau) + extra * tauRange))
        axScatter.set_ylim((np.min(yClean) - extra * yRange, np.max(yClean) + extra * yRange))


    #define colors of the 3 classes
    colors=['tomato', 'lightgreen','skyblue']
    linecolors=['red','green', 'blue']

    # the scatter plot:
    scatters=[None]*3
    def PlotEvents(visBool):
        #clear axes
        axScatter.clear()
        axHistx.clear()
        axHisty.clear()

        #lists for setting the limits
        alltau=[]
        allyVals=[]
        for i in range(len(catEvents)):
            #If checkbox is True, plot events
            if visBool[i]:

                #Extract Tau and Y
                events=catEvents[i]
                tau,yVals = extractytau(events)
                scatters[i] = axScatter.scatter(tau, yVals, color=colors[i], marker='o', s=30,
                                             linewidths=0.1,edgecolors=linecolors[i], picker=5,visible=visBool[i])  # added some stuff here to improve aesthetics
                #axScatter.set_xscale('log')
                axHistx.hist(tau, bins=50, color=colors[i], visible=visBool[i])
                axHisty.hist(yVals, bins=50, orientation='horizontal', color=colors[i], visible=visBool[i])

                alltau.extend(tau)
                allyVals.extend(yVals)

        #set limits
        if len(alltau)>0:
            limits(alltau,allyVals,axScatter)
        setlabels()
        plt.title('{} number of events'.format(len(alltau)))

    PlotEvents(visBool)

    #If click on checkbox, switch Boolean and replot events
    def func(label):
        for i in range(len(labelsCheckBox)):
            if label == labelsCheckBox[i]:
                visBool[i] = not visBool[i]
        PlotEvents(visBool)
        plt.draw()


    #WHen clicking on event
    def onpick(event):
        for i in range(len(catEvents)):
            if event.artist == scatters[i]:

                N = len(event.ind)
                if not N: return True
                figi = plt.figure(figsize=(10, 6))
                for subplotnum, dataind in enumerate(event.ind):
                    ax = figi.add_subplot(N, 1, subplotnum + 1)
                    PlotEvent(catEvents[i][dataind], ax, savefile, showCUSUM)
                figi.show()
        return True

    fig.canvas.mpl_connect('pick_event', onpick)
    check.on_clicked(func)

    plt.show()

def PlotEvent(event,ax=None, savefile=os.getcwd(), showCUSUM=False):

    #Link event to axes to keep it around

    if ax is None:
        #plt.figure(figsize=(10, 6))
        fig, ax = plt.subplots(figsize=(10, 6))
    ax._event=event

    def SavePlot(eventMouse):
        # Check if directory exists
        directory = os.path.dirname(savefile)
        savename=directory + os.sep + 'event.pdf'
        i=1

        while os.path.exists(directory + os.sep + 'event_{}.pdf'.format(i)):
            i += 1
        savename = directory + os.sep + 'event_{}.pdf'.format(i)

        eventMouse.inaxes.figure.savefig(savename, transparent=True)

    def ShowFullTrace(eventMouse):
        event=eventMouse.inaxes.figure.axes[0]._event
        ShowEventInTrace(event)




    if showCUSUM and hasattr(event, 'changeTimes') and len(event.changeTimes)>2:
        eventLength = event.eventLengthCUSUM
        currentDrop = event.currentDropCUSUM
    else:
        showCUSUM=False
        eventLength = event.eventLength
        currentDrop = event.currentDrop


    fn=filename_w_ext = os.path.basename(event.filename)

    plotTitle = fn + '\n' + 'Event length: {}\nConductance drop: {} with voltage {}'.format(Time.format_data(eventLength), Cond.format_data(currentDrop/event.voltage),Volt.format_data(event.voltage))

    ax.set_xlabel('time (s)')
    ax.set_ylabel('current (A)')
    ax.xaxis.set_major_formatter(Time)
    ax.yaxis.set_major_formatter(Amp)

    if plotTitle:
        plt.title(plotTitle)


    #Add buttons

    #Save button
    bax = plt.axes([0.77, 0.95, 0.15, 0.03])
    bsave = Button(bax, 'Save figure')
    bsave.on_clicked(SavePlot)
    #Link button to axes to preserve function
    ax._bsave = bsave

    #Show original trace button
    bax2 = plt.axes([0.77, 0.9, 0.15, 0.03])
    bfull = Button(bax2, 'Show original Trace')
    # Link button to axes to preserve function
    ax._bfull = bfull
    bfull.on_clicked(ShowFullTrace)

    #Plotting
    timeVals1 = np.linspace(0, len(event.before) / event.samplerate, num=len(event.before))
    timeVals2 = np.linspace(0 + max(timeVals1), len(event.eventTrace) / event.samplerate + max(timeVals1),
                            num=len(event.eventTrace))
    timeVals3 = np.linspace(0 + max(timeVals2), len(event.after) / event.samplerate + max(timeVals2), num=len(event.after))




    ax.plot(np.append(timeVals1,timeVals2[0]), np.append(event.before,event.eventTrace[0]), color='tomato')
    ax.plot(timeVals2, event.eventTrace, color='mediumslateblue')
    ax.plot(np.append(timeVals2[-1],timeVals3), np.append(event.eventTrace[-1],event.after), color='tomato')

    if showCUSUM:
        timeVals = np.linspace(0, len(event.segmentedSignal) / event.samplerate, num=len(event.segmentedSignal))

        if hasattr(event,'mcbefore') and hasattr(event,'mcafter') and hasattr(event,'mctrace'):
            ax.plot(timeVals1, event.mcbefore,'--', color='tomato')
            x=np.append(np.append(timeVals1[-1],timeVals2),timeVals3[0])
            y=np.append(np.append(event.mcbefore[-1],event.mctrace),event.mcafter[0])
            ax.plot(x, y, color='yellow')
            ax.plot(timeVals3, event.mcafter,'--', color='tomato')
        else:
            ax.plot(timeVals,event.segmentedSignal, color='yellow') #,timeVals3[0],event.mcafter[0]
    else:
        beforeBaseline=np.full(len(event.before), event.baseline)
        ax.plot(timeVals1,beforeBaseline, '--', color='tomato')
        afterBaseline = np.full(len(event.after), event.baseline)
        ax.plot(timeVals3,afterBaseline, '--', color='tomato')

        meanTrace = np.full(len(event.eventTrace), event.baseline-event.currentDrop)
        ax.plot(timeVals2,meanTrace, '--', color='mediumslateblue')

    if 'fig' in locals():
         plt.show()

def ShowEventInTrace(event):
    filename=event.filename
    loadedData = LoadData.OpenFile(filename,1e3,True) #, ChimeraLowPass, True, CutTraces)

    fig, ax = plt.subplots(figsize=(10, 6))

    FullTrace = loadedData['i1']

    times=np.linspace(0, len(FullTrace) / event.samplerate, num=len(FullTrace))
    ax.plot(times,FullTrace,zorder=1)

    ax.set_xlabel('time (s)')
    ax.set_ylabel('current (A)')
    ax.xaxis.set_major_formatter(Time)
    ax.yaxis.set_major_formatter(Amp)


    # Create a Rectangle patch
    if hasattr(event,'changeTimes') and len(event.changeTimes)>2:
        start_i = (event.beginEventCUSUM - len(event.before))/event.samplerate
        end_i = (event.endEventCUSUM + len(event.after))/event.samplerate
    else:
        start_i = (event.beginEvent - len(event.before))/event.samplerate
        end_i = (event.endEvent + len(event.after))/event.samplerate
    minE=np.min(np.append(np.append(event.eventTrace,event.before),event.after))
    maxE=np.max(np.append(np.append(event.eventTrace,event.before),event.after))
    rect = patches.Rectangle((start_i, minE-0.1*(maxE-minE)), end_i-start_i, maxE+0.2*(maxE-minE)-minE, linestyle='--', linewidth=1, edgecolor='r', facecolor='none',zorder=10)

    # Add the patch to the Axes
    ax.add_patch(rect)

    plt.title(os.path.basename(filename))

    plt.show()


def PlotCurrentTrace(currentTrace, samplerate):
    timeVals = np.linspace(0, len(currentTrace) / samplerate, num=len(currentTrace))
    fig,ax=plt.subplots(figsize=(10, 6))

    ax.plot(timeVals, currentTrace)
    ax.set_xlabel('time (s)')
    ax.set_ylabel('current (A)')
    ax.xaxis.set_major_formatter(Time)
    ax.yaxis.set_major_formatter(Amp)

    plt.show()


def PlotCurrentTraceBaseline(before, currentTrace, after, samplerate, plotTitle=''):
    timeVals1 = np.linspace(0, len(before) / samplerate, num=len(before))
    timeVals2 = np.linspace(0 + max(timeVals1), len(currentTrace) / samplerate + max(timeVals1), num=len(currentTrace))
    timeVals3 = np.linspace(0 + max(timeVals2), len(after) / samplerate + max(timeVals2), num=len(after))

    #plt.figure(figsize=(10, 6))
    fig,ax=plt.subplots(figsize=(10, 6))

    ax.plot(timeVals1, before, color='red')
    ax.plot(timeVals2, currentTrace)
    ax.plot(timeVals3, after, color='red')
    ax.set_xlabel('time (s)')
    ax.set_ylabel('current (A)')
    ax.xaxis.set_major_formatter(Time)
    ax.yaxis.set_major_formatter(Amp)

    if plotTitle:
        plt.title(plotTitle)

    plt.show()

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input file')

    args = parser.parse_args()
    inputData=args.input
    if inputData==None:
        inputData=askopenfilenames(filetypes=[('data files', 'Data*.dat')])
        if inputData:
            inputData=os.path.splitext(inputData[0])[0]

    if inputData:
        shelfFile=shelve.open(inputData)
        translocationEvents=shelfFile['translocationEvents']
        shelfFile.close()

        PlotG_tau(translocationEvents.events,inputData)
