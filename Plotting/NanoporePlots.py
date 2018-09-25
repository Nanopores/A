import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
import argparse
from tkinter.filedialog import askopenfilenames,askdirectory
import shelve
from matplotlib.ticker import EngFormatter
import os
from pprint import pprint
from matplotlib.widgets import CheckButtons


Amp = EngFormatter(unit='A', places=2)
Time = EngFormatter(unit='s', places=2)
Volt = EngFormatter(unit='V', places=2)
Cond = EngFormatter(unit='S', places=2)


def PlotG_tau(events, savefile, showCurrentInstead=False, normalized=False,showCUSUM=True):


    #Filter events
    if showCUSUM:
        indices, filteredevents = zip(*[(ind, event) for ind, event in enumerate(events) if hasattr(event, 'changeTimes') and len(event.changeTimes) > 2 ])
    else:
        filteredevents=events

    if showCurrentInstead:
        yVals = [event.currentDrop for event in filteredevents]
        tau = [event.eventLength for event in filteredevents]
    else:
        shownevents = [ event for event in filteredevents if event.voltage > 0 ]
        yVals = [event.currentDrop / event.voltage for event in filteredevents if event.voltage > 0]
        tau = [event.eventLength for event in filteredevents if event.voltage > 0]

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

    # rax = plt.axes([0.75, 0.8, 0.14, 0.15])
    # check = CheckButtons(rax, ('CUSUM-fitted', 'Not fitted', 'impulse'), (True, False, False))

    # the scatter plot:
    scatters = axScatter.scatter(tau, yVals, color='tomato', marker='o', s=30,
                                 linewidths=0.1,edgecolors='red', picker=5)  # added some stuff here to improve aesthetics

    # now determine nice limits by hand:
    extra = 0.1  # 0.1 = 10%
    tauRange = np.max(tau) - np.min(tau)
    yClean = [y for y in yVals if str(y) != 'nan']
    yRange = np.max(yClean) - np.min(yClean)
    #
    # tauLim = (int(tauMax/binwidth) + 1) * binwidth
    # currentLim = (int(currentMax/binwidth) + 1) * binwidth

    axScatter.set_xlim((np.min(tau) - extra * tauRange, np.max(tau) + extra * tauRange))
    axScatter.set_ylim((np.min(yClean) - extra * yRange, np.max(yClean) + extra * yRange))

    # tauBins = np.arange(-tauLim, tauLim + binwidth, binwidth)
    # currentBins = np.arange(-currentLim, currentLim + binwidth, binwidth)
    axHistx.hist(tau, bins=50, color='tomato')
    plt.setp(axHistx.get_xticklabels(), visible=False)
    axHisty.hist(yClean, bins=50, orientation='horizontal',
                 color='tomato')  # increased binning here to avoid overlapping bars in the histogram, added color option
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

    # Check if directory exists
    directory = os.path.dirname(savefile)
    if showCurrentInstead:
        fig.savefig(directory + os.sep + 'PlotITau.pdf', transparent=True)
    else:
        fig.savefig(directory + os.sep + 'PlotGTau.pdf', transparent=True)

    def onpick(event):
        if event.artist != scatters: return True

        N = len(event.ind)
        if not N: return True
        for subplotnum, dataind in enumerate(event.ind):
            figi = plt.figure()
            ax = figi.add_subplot(N, 1, subplotnum + 1)
            PlotEvent(shownevents[dataind], ax, showCUSUM)
            figi.show()

        # figi = plt.figure()

        # for subplotnum, dataind in enumerate(event.ind):
        #     ax = figi.add_subplot(N, 1, subplotnum + 1)
        #     ax.plot(X[dataind])
        #     ax.text(0.05, 0.9, 'mu=%1.3f\nsigma=%1.3f' % (xs[dataind], ys[dataind]),
        #             transform=ax.transAxes, va='top')
        #     ax.set_ylim(-0.5, 1.5)
        # figi.show()


        return True

    fig.canvas.mpl_connect('pick_event', onpick)
    plt.show()

def PlotEvent(event,ax=None,showCUSUM=False):

    if ax is None:
        plt.figure(figsize=(10, 6))
        fig, ax = plt.subplots(figsize=(10, 6))


    if showCUSUM and hasattr(event, 'changeTimes') and len(event.changeTimes)>2:
        eventLength = event.eventLengthCUSUM
        currentDrop = event.currentDropCUSUM
    else:
        showCUSUM=False
        eventLength = event.eventLength
        currentDrop = event.currentDrop


    part1=np.append(event.before,event.eventTrace)

    fn=filename_w_ext = os.path.basename(event.filename)

    plotTitle = fn + '\n' + 'Event length: {}\nCurrent drop: {}'.format(Time.format_data(eventLength), Amp.format_data(currentDrop))
    #PlotCurrentTraceBaseline(self.before, self.eventTrace, self.after, self.samplerate, titleplot)

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

    ax.set_xlabel('time (s)')
    ax.set_ylabel('current (A)')
    ax.xaxis.set_major_formatter(Time)
    ax.yaxis.set_major_formatter(Amp)

    if plotTitle:
        plt.title(plotTitle)

    if 'fig' in locals():
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
        inputData=askopenfilenames()

    shelfFile=shelve.open(inputData)
    translocationEvents=shelfFile['TranslocationEvents']
    shelfFile.close()

    PlotG_tau(translocationEvents.events,inputData)