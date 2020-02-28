# -*- coding: utf-8 -*-
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import EngFormatter
from matplotlib.widgets import CheckButtons, Button
import argparse
import platform
import shelve
import os
import Functions
import LoadData
import NanoporeClasses as NC

from holoviews.operation import histogram
import holoviews as hv
from holoviews import opts

import seaborn as sns
sns.set()

hv.extension('bokeh')

# needed for plotting in jupyter
from bokeh.resources import INLINE
import bokeh.io
bokeh.io.output_notebook(INLINE)

# temporary, ignore bokeh warnings
import warnings
warnings.simplefilter(action='ignore', category=DeprecationWarning)


Amp = EngFormatter(unit='A', places=2)
Time = EngFormatter(unit='s', places=2)
Volt = EngFormatter(unit='V', places=2)
Cond = EngFormatter(unit='S', places=2)


# Extract y and tau out of events
def extractytau(filteredevents, showCurrent, normalized = False):
    if showCurrent or normalized:
        if normalized:
            yVals = [event.currentDrop/event.baseline for event in filteredevents]
        else:
            yVals = [event.currentDrop for event in filteredevents]
        tau = [event.eventLength for event in filteredevents]
    else:
        yVals = [event.currentDrop / event.voltage for event in filteredevents if event.voltage > 0]
        tau = [event.eventLength for event in filteredevents if event.voltage > 0]
    return tau, yVals

def PlotGTau(eventClass, xLim = None, yLim = None, showCurrent = False, normalized = False, showLog=False, showNonFitted=False):
    """
    Function used to produce scatter plots and histograms of the events.
    The figure produced has 3 subplots:
        Up: Histogram with the number of events per event length
        Right: Histogram with the number of events per conductance drop [nS]
        Center: Scatter-plot of all the events wiht in on x-coordinates the event length [s]
            and in y-coordinates the conductance drop [nS].

    In the plots, the events were distributed into the 3 types of events:
        CUSUM-fitted in red ('Real' type)
        Impulse in blue ('Impulse' type)
        Non-fitted in green ('Rough' type)

    Parameters
    ----------
    eventClass : AllEvents object
        All the events to be plotted.
    xLim : 2D list, optional
        limits on x axis
    yLim : 2D list, optional
        limits on y axis
    showCurrent : bool, optional
        False by default. If True, it will change the SI unit in the y-axis from siemens [S] to ampers [A].
        So instead of units in conductance drop, it will have current drop.
    normalized : bool, optional
        False by default. If True, it will change in the y-axis the unit from siemens [S] to normalized current drop without unit.
    showLog :   bool, optional
        False by default. If True, it will change the axis to logarithmic
    showNonFitted:  bool, optional
        False by default. If True, it will show non fitted events

    """
    # categorize events in three types
    # CUSUM fitted events
    CUSUMEvents = eventClass.GetEventTypes('CUSUM')

    if len(CUSUMEvents) == 0:
        CUSUMEvents = eventClass.GetEventTypes('Real')

    # Non-fitted events
    if showNonFitted:
        nonFittedEvents = eventClass.GetEventTypes('Rough')
    else:
        nonFittedEvents = []

    # Impulse events
    impulseEvents = eventClass.GetEventTypes('Impulse')

    # select min max
    events = [event for event in eventClass.events if showCurrent or event.voltage is not 0]
    allTau, allYVals = extractytau(events, showCurrent, normalized)

    maxy = 1 if normalized else 1.1 * np.percentile(allYVals, 99)
    minx = 0.9 * np.percentile(allTau, 1) if showLog else 0

    brx = (minx, 1.1 * np.percentile(allTau, 99)) if xLim is None else tuple(xLim)
    bry = (0, maxy) if yLim is None else tuple(yLim)

    br = dict()
    br['x'] = brx
    br['y'] = bry

    showlogHist = dict()
    showlogHist['x'] = showLog
    showlogHist['y'] = False

    xlabel = 'Dwell time (s)'
    ylabel = 'Current drop (A)' if showCurrent else 'Conductance drop (S)'

    tau1, yVals1 = extractytau(impulseEvents, showCurrent, normalized)
    points = hv.Points((tau1, yVals1), label='impulse')

    tau2, yVals2 = extractytau(CUSUMEvents, showCurrent, normalized)
    points2 = hv.Points((tau2, yVals2), label='CUSUM')

    tau3, yVals3 = extractytau(nonFittedEvents, showCurrent, normalized)
    points3 = hv.Points((tau3, yVals3), label='non-fitted')

    xhist, yhist = (histogram(points3, bin_range=br[dim], dimension=dim, log=(showlogHist[dim] and len(tau3) > 0),normed=False) *
                    histogram(points2, bin_range=br[dim], dimension=dim, log=(showlogHist[dim] and len(tau2) > 0),normed=False) *
                    histogram(points, bin_range=br[dim], dimension=dim, log=(showlogHist[dim] and len(tau1) > 0),normed=False)
                    for dim in 'xy')

    return ((points3 * points2 * points).opts(logx=showLog, xlabel=xlabel, ylabel=ylabel, width=500, height=500, xlim=brx,
                                       ylim=bry) << yhist.opts(
        width=200) << xhist.opts(height=150, logx=showLog)).opts(
        opts.Histogram(xlabel='', ylabel='', alpha=0.3, show_legend=False))


def PlotG_tau(translocationEvents, fig=None, savefile=None, showCurrent=False, normalized=False, showCUSUM=True, showLog=False):
    """
    Function used to produce scatter plots and histograms of the events.
    The figure produced has 3 subplots:
        Up: Histogram with the number of events per event length
        Right: Histogram with the number of events per conductance drop [nS]
        Center: Scatter-plot of all the events wiht in on x-coordinates the event length [s]
            and in y-coordinates the conductance drop [nS].

    In the plots, the events were distributed into the 3 types of events:
        CUSUM-fitted in red ('Real' type)
        Non-fitted in blue ('Rough' type)
        Impulse in green ('Impulse' type)

    Parameters
    ----------
    events : AllEvents object
        All the events to be plotted.
    savefile : str, optional
        Full path to file where the plots will be saved.
    showCurrent : bool, optional
        False by default. If True, it will change the SI unit in the y-axis from siemens [S] to ampers [A].
        So instead of units in conductance drop, it will have current drop.
    normalized : bool, optional
        False by default. If True, it will change in the y-axis the unit from siemens [S] to normalized current drop without unit.
    showCUSUM : bool, optional
        True by default.

    """

    # categorize events in three types
    # CUSUM fitted events
    CUSUMEvents = translocationEvents.GetEventTypes('CUSUM')

    if len(CUSUMEvents) == 0:
        CUSUMEvents = translocationEvents.GetEventTypes('Real')

    # Non-fitted events
    nonFittedEvents = translocationEvents.GetEventTypes('Rough')

    # Impulse events
    impulseEvents = translocationEvents.GetEventTypes('Impulse')

    catEvents = (CUSUMEvents, nonFittedEvents, impulseEvents)

    # Save figure
    def SavePlot(event):
        # Check if directory exists
        directory = os.path.dirname(savefile)
        if showCurrent:
            fig.savefig(directory + os.sep + 'PlotITau.pdf', transparent=True)
        else:
            fig.savefig(directory + os.sep + 'PlotGTau.pdf', transparent=True)

    def ShowLog(event):
        name = axScatter.xaxis._scale.name
        bshowlog.label.set_text('Show ' + name)
        axScatter.set_xscale('linear') if name == 'log' else axScatter.set_xscale('log')

    # definitions for the axes
    left, width = 0.15, 0.55
    bottom, height = 0.1, 0.6
    left_h = left + width + 0.015
    bottom_h = bottom + height + 0.015

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    if fig is None:
        fig = plt.figure(1, figsize=(10, 8))

    # define axes and link histogram to scatterplot
    axScatter = fig.add_axes(rect_scatter)

    axHistx = fig.add_axes(rect_histx, sharex=axScatter)
    axHisty = fig.add_axes(rect_histy, sharey=axScatter)

    # Checkboxes to turn on or off events
    rax = plt.axes([0.75, 0.73, 0.14, 0.15])
    visBool = [True, False, True]
    labelsCheckBox = ('CUSUM-fitted', 'Not fitted', 'impulse')
    check = CheckButtons(rax, labelsCheckBox, visBool)
    # Link button to axes to preserve function
    rax._check = check

    bax = plt.axes([0.77, 0.9, 0.1, 0.03])
    bnext = Button(bax, 'Save figure')
    bnext.on_clicked(SavePlot)
    # Link button to axes to preserve function
    bax._bnext = bnext

    baxl = plt.axes([0.77, 0.95, 0.1, 0.03])
    txt = 'log' if showLog else 'linear'
    bshowlog = Button(baxl, 'Show ' + txt)
    bshowlog.on_clicked(ShowLog)
    # Link button to axes to preserve function
    bax._bshowlog = bshowlog

    # Show labels
    def setlabels():
        plt.setp(axHistx.get_xticklabels(), visible=False)
        plt.setp(axHisty.get_yticklabels(), visible=False)

        axScatter.set_xlabel('Event length (s)')
        axScatter.xaxis.set_major_formatter(Time)
        if normalized:
            axScatter.set_ylabel('current drop (normalized)')
        else:
            if showCurrent:
                axScatter.set_ylabel('current drop (A)')
                axScatter.yaxis.set_major_formatter(Amp)
            else:
                axScatter.set_ylabel('Conductance drop (G)')
                axScatter.yaxis.set_major_formatter(Cond)

    # Set limits
    allEvents = [event for event in translocationEvents.events if showCurrent or event.voltage is not 0]

    bins = 100
    extra = 0.1  # 0.1 = 10%
    allTau, allYVals = extractytau(allEvents, showCurrent)
    yClean = [y for y in allYVals if str(y) != 'nan']
    taurangeHist = np.linspace(min(allTau), max(allTau), num=bins)
    yValsrangeHist = np.linspace(min(yClean), max(yClean), num=bins)

    # define colors of the 3 classes
    colors = ['tomato', 'lightgreen', 'skyblue']
    linecolors = ['red', 'green', 'blue']

    # the scatter plot:
    scatters = [None]*3

    def PlotEvents(visBool):
        # clear axes
        axScatter.clear()
        axHistx.clear()
        axHisty.clear()
        nrEvents = 0

        for i in range(len(catEvents)):
            # If checkbox is True, plot events
            if visBool[i]:

                # Extract Tau and Y
                events = catEvents[i]
                tau, yVals = extractytau(events, showCurrent)
                scatters[i] = axScatter.scatter(tau, yVals, color=colors[i], marker='o', s=30,
                                                linewidths=0.1, edgecolors=linecolors[i], picker=5, visible=visBool[i])

                axHistx.hist(tau, bins=taurangeHist, color=colors[i], visible=visBool[i])
                axHisty.hist(yVals, bins=yValsrangeHist, orientation='horizontal', color=colors[i], visible=visBool[i])
                nrEvents += len(events)

        if showLog:
            axScatter.set_xscale('log')

        # set limits
        if len(allTau) > 0:
            tauRange = np.max(allTau) - np.min(allTau)
            yRange = np.max(yClean) - np.min(yClean)
            axScatter.set_xlim((np.max([np.min(allTau) - extra * tauRange, 1e-6]), np.max(allTau) + extra * tauRange))
            axScatter.set_ylim((np.min(yClean) - extra * yRange, np.max(yClean) + extra * yRange))

        setlabels()
        fig.suptitle('{} number of events'.format(nrEvents))

    PlotEvents(visBool)

    # If click on checkbox, switch Boolean and replot events
    def func(label):
        for i in range(len(labelsCheckBox)):
            if label == labelsCheckBox[i]:
                visBool[i] = not visBool[i]
        PlotEvents(visBool)
        plt.draw()

    # When clicking on event
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

    if not hasattr(fig, '_AnalysisUI__attached'):
        plt.show()


def PlotGTauVoltage(eventClass, bins=20, xLim=None, yLim=None, showCurrent=False,
                    voltageLimits=None, showLog=False, heatmap=False):
    voltagesList = eventClass.GetAllVoltages()
    if voltageLimits is not None:
        voltagesList = [x for x in voltagesList if x >= voltageLimits[0] and x <= voltageLimits[1]]

    # select min max
    events = [event for event in eventClass.events if (showCurrent or event.voltage is not 0) and event.voltage in voltagesList]
    allTau, allYVals = extractytau(events, showCurrent)

    maxy = 1.1 * np.percentile(allYVals, 99)
    minx = 0.9 * np.percentile(allTau, 1) if showLog else 0

    brx = (minx, 1.1 * np.percentile(allTau, 99)) if xLim is None else tuple(xLim)
    bry = (0, maxy) if yLim is None else tuple(yLim)
    brxbins = np.linspace(brx[0], brx[-1], bins).tolist()
    brybins = np.linspace(bry[0], bry[-1], bins).tolist()

    # define of variables
    colors = sns.color_palette(n_colors=len(voltagesList))
    colorlist = ["#%02x%02x%02x" % (int(color[0]*255), int(color[1]*255), int(color[2]*255)) for color in colors]

    xlabel = 'Dwell time (s)'
    ylabel = 'Current drop (A)' if showCurrent else 'Conductance drop (S)'

    # Filter points specifically for the heatmap
    allTauf = [xval for (xval, yval) in zip(allTau, allYVals) if
               xval > brx[0] and xval < brx[-1] and yval > bry[0] and yval < bry[-1]]
    allYValsf = [yval for (xval, yval) in zip(allTau, allYVals) if
               xval > brx[0] and xval < brx[-1] and yval > bry[0] and yval < bry[-1]]

    hexT = hv.HexTiles((allTauf, allYValsf))

    for i, voltage in enumerate(voltagesList):
        selectEvents = eventClass.GetEventsforVoltages(voltage)
        tau, yVals = extractytau(selectEvents, showCurrent)

        points_ = hv.Points((tau, yVals), label=str(Volt.format_data(voltage)))
        xhist_ = histogram(points_, bins=brxbins, dimension='x', log=(showLog and len(tau) > 0), normed=False)
        yhist_ = histogram(points_, bins=brybins, dimension='y', normed=False)
        if i == 0:
            points = points_
            xhist = xhist_
            yhist = yhist_
        else:
            points = points * points_
            xhist = xhist * xhist_
            yhist = yhist * yhist_

    if heatmap and not showLog:
        return ((hexT*points).opts(logx=showLog, xlabel=xlabel, ylabel=ylabel, width=500, height=500, xlim=brx, ylim=bry) \
               << yhist.opts(width=200) << xhist.opts(height=150, logx=showLog)).opts(
            opts.Histogram(fill_color=hv.Cycle(colorlist), xlabel='', ylabel='', alpha=0.3, show_legend=False),
            opts.HexTiles(min_count=0, tools=['hover']),
            opts.Points(color=hv.Cycle(colorlist), alpha=1, size=1, muted_alpha=.1))
    else:
        return (points.opts(logx=showLog, xlabel=xlabel, ylabel=ylabel, width=500, height=500, xlim=brx, ylim=bry)\
            << yhist.opts(width=200) << xhist.opts(height=150, logx=showLog)).opts(
            opts.Histogram(fill_color=hv.Cycle(colorlist), xlabel='', ylabel='', alpha=0.3, show_legend=False),
            opts.Points(color=hv.Cycle(colorlist), alpha=0.4, size=6, muted_alpha=.05))

def PlotEvent(event, ax=None, savefile=os.getcwd(), showCUSUM=True, showCurrent=False, showButtons = True, axisFormatter = True, plotTitleBool = True):
    """
    Function used to plot a single event passed in argument. The event will be represented
    in a blue trace and the baseline in a red trace.

    Parameters
    ----------
    event : TranslocationEvent object
        Event to be plotted.
    ax :axes.Axes object or array of Axes object
        Axes of the plot.
    savefile : str, optional
        By default, the current working directory. Full path to the directory to save the plot.
    showCUSUM : bool, optional
        False by default. If True, it will plot the CUSUM-fit overlayed in yellow.

    """

    #Link event to axes to keep it around

    if ax is None:
        #plt.figure(figsize=(10, 6))
        fig, ax = plt.subplots(figsize=(10, 6))
    ax._event = event

    def SavePlot(eventMouse):
        # Check if directory exists
        directory = os.path.dirname(savefile)
        savename=directory + os.sep + 'event.pdf'
        i = 1

        while os.path.exists(directory + os.sep + 'event_{}.pdf'.format(i)):
            i += 1
        savename = directory + os.sep + 'event_{}.pdf'.format(i)

        eventMouse.inaxes.figure.savefig(savename, transparent=True)

    def ShowFullTrace(eventMouse):
        event=eventMouse.inaxes.figure.axes[0]._event
        ShowEventInTrace(event)


    if showCUSUM and hasattr(event, 'changeTimes') and len(event.changeTimes)>=2:
        eventLength = event.eventLengthCUSUM
        currentDrop = event.currentDropCUSUM
    else:
        showCUSUM=False
        eventLength = event.eventLength
        currentDrop = event.currentDrop


    fn=filename_w_ext = os.path.basename(event.filename)
    if showCurrent:
        plotTitle = fn + '\n' + 'Event length: {}\nCurrent drop: {} with voltage {}'.format(
            Time.format_data(eventLength), Cond.format_data(currentDrop),
            Volt.format_data(event.voltage))
    else:
        plotTitle = fn + '\n' + 'Event length: {}\nConductance drop: {} with voltage {}'.\
            format(Time.format_data(eventLength), Cond.format_data(currentDrop/event.voltage),Volt.format_data(event.voltage))

    ax.set_xlabel('time (s)')
    ax.set_ylabel('current (A)')
    if axisFormatter:
        ax.xaxis.set_major_formatter(Time)
        ax.yaxis.set_major_formatter(Amp)

    if plotTitleBool:
        plt.title(plotTitle)

    if showButtons:
        # Add buttons

        # Save button
        bax = plt.axes([0.77, 0.95, 0.15, 0.03])
        bsave = Button(bax, 'Save figure')
        bsave.on_clicked(SavePlot)
        # Link button to axes to preserve function
        ax._bsave = bsave

        # Show original trace button
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

    ax.plot(np.append(timeVals1, timeVals2[0]), np.append(event.before, event.eventTrace[0]), color='tomato')
    ax.plot(timeVals2, event.eventTrace, color='mediumslateblue')
    ax.plot(np.append(timeVals2[-1], timeVals3), np.append(event.eventTrace[-1], event.after), color='tomato')

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

def ShowEventInTrace_SignalPreloaded(FullTrace, AllData, eventnumber, ax, line = None, firstCall = True, dscorrection = None):
    times = np.linspace(0, len(FullTrace) / AllData.events[eventnumber].samplerate, num=len(FullTrace))
    if line:
        line.set_ydata(FullTrace)
        line.set_xdata(times)
        if firstCall:
            ax.set_xlim(np.min(times), np.max(times))
            ax.set_ylim(np.min(FullTrace), np.max(FullTrace))
        print('updated lines')
    else:
        ax.plot(times, FullTrace, zorder=1)
        print('Updated plot')

    ax.set_xlabel('time (s)')
    ax.set_ylabel('current (A)')

    # Create a Rectangle patch
    if dscorrection:
        if hasattr(AllData.events[eventnumber], 'changeTimes') and len(AllData.events[eventnumber].changeTimes) > 2:
            start_i = (AllData.events[eventnumber].beginEventCUSUM - len(AllData.events[eventnumber].before)) / AllData.events[eventnumber].samplerate * dscorrection
            end_i = (AllData.events[eventnumber].endEventCUSUM + len(AllData.events[eventnumber].after)) / AllData.events[eventnumber].samplerate * dscorrection
        else:
            start_i = (AllData.events[eventnumber].beginEvent - len(AllData.events[eventnumber].before)) / AllData.events[eventnumber].samplerate * dscorrection
            end_i = (AllData.events[eventnumber].endEvent + len(AllData.events[eventnumber].after)) / AllData.events[eventnumber].samplerate * dscorrection
    else:
        if hasattr(AllData.events[eventnumber], 'changeTimes') and len(AllData.events[eventnumber].changeTimes) > 2:
            start_i = (AllData.events[eventnumber].beginEventCUSUM - len(AllData.events[eventnumber].before)) / \
                      AllData.events[eventnumber].samplerate
            end_i = (AllData.events[eventnumber].endEventCUSUM + len(AllData.events[eventnumber].after)) / \
                    AllData.events[eventnumber].samplerate
        else:
            start_i = (AllData.events[eventnumber].beginEvent - len(AllData.events[eventnumber].before)) / \
                      AllData.events[eventnumber].samplerate
            end_i = (AllData.events[eventnumber].endEvent + len(AllData.events[eventnumber].after)) / AllData.events[
                eventnumber].samplerate

    minE = np.min(np.append(np.append(AllData.events[eventnumber].eventTrace, AllData.events[eventnumber].before), AllData.events[eventnumber].after))
    maxE = np.max(np.append(np.append(AllData.events[eventnumber].eventTrace, AllData.events[eventnumber].before), AllData.events[eventnumber].after))
    rect = patches.Rectangle((start_i, minE - 0.1 * (maxE - minE)), end_i - start_i, maxE + 0.2 * (maxE - minE) - minE,
                             linestyle='--', linewidth=1, edgecolor='r', facecolor='none', zorder=10)
    # Add the patch to the Axes
    print(ax.patches.clear())
    ax.add_patch(rect)

def ShowEventInTrace(event, lowPass = 1e3):
    """
    Function used to show the event with it's location framed in red in the original full signal trace in blue.

    Parameters
    ----------
    event : TranslocationEvent object
        Event to be plotted.
        LowPass low pass filter, default set to 1 kHz

    """

    filename = event.filename
    loadedData = LoadData.OpenFile(filename, lowPass, True)

    fig, ax = plt.subplots(figsize=(10, 6))

    if loadedData['samplerate'] > lowPass:
        output = Functions.LowPass(loadedData['i1'], loadedData['samplerate'], lowPass)
        FullTrace = output['data']
        samplerate = output['samplerate']
    else:
        FullTrace = loadedData['i1']
        samplerate = loadedData['samplerate']

    times = np.linspace(0, len(FullTrace) / samplerate, num=len(FullTrace))
    ax.plot(times, FullTrace, zorder=1)

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
    """
    Function used in TranslocationEvent class methods in the NanoporeClasses module to plot events.
    It will plot the TranslocationEvent's currrentTrace passed in argument.

    Parameters
    ----------
    currentTrace : list of float
        Data points in current to be plotted.
    samplerate : float
        Sampling frequency of the data acquisition.

    """

    timeVals = np.linspace(0, len(currentTrace) / samplerate, num=len(currentTrace))
    fig,ax=plt.subplots(figsize=(10, 6))

    ax.plot(timeVals, currentTrace)
    ax.set_xlabel('time (s)')
    ax.set_ylabel('current (A)')
    ax.xaxis.set_major_formatter(Time)
    ax.yaxis.set_major_formatter(Amp)

    plt.show()


def PlotCurrentTraceBaseline(before, currentTrace, after, samplerate, plotTitle=''):
    """
    Function used in TranslocationEvent class methods in the NanoporeClasses module to plot events.
    It will plot TranslocationEvent's currentTrace and the surrounding baseline (after and before)
    passed in argument.

    Parameters
    ----------
    before : list of float
        Data points in current of the baseline trace before the event.
    currentTrace : list of float
        Data points in current of the event trace.
    after : list of float
        Data points in current of the baseline trace after the event.
    samplerate : float
        Sampling frequency of the data aquisition.
    plotTitle : str, optional
        Plot title.

    """

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
    from tkinter import Tk
    from tkinter.filedialog import askopenfilenames, askdirectory

    matplotlib.use('TkAgg')

    if (platform.system() == 'Darwin'):
        root = Tk()
        root.withdraw()

    if (platform.system() == 'Darwin'):
        os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input file')

    args = parser.parse_args()
    inputData = args.input
    if inputData == None:
        inputData = askopenfilenames(filetypes=[('data files', '*.dat')])  #for Mac systems, replace 'Data*.dat' with >> '*.dat'
        #if inputData:
        #    inputData=os.path.splitext(inputData[0])[0]
        if (platform.system() == 'Darwin'):
            root.update()

    translocationEvents = NC.AllEvents()
    if inputData:
        for filename in inputData:
            shelfFile = shelve.open(os.path.splitext(filename)[0])
            translocationEventstemp = shelfFile['TranslocationEvents']
            shelfFile.close()
            translocationEvents.AddEvent(translocationEventstemp)
        PlotG_tau(translocationEvents, savefile=inputData)


