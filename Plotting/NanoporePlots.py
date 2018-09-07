import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter


Amp = EngFormatter(unit='A', places=2)
Time = EngFormatter(unit='s', places=2)
Volt = EngFormatter(unit='V', places=2)


def PlotI_tau(current,tau, normalized=False):
    # definitions for the axes
    left, width = 0.15, 0.55
    bottom, height = 0.1, 0.6
    left_h = left + width + 0.015
    bottom_h = bottom + height + 0.015

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(10, 8))

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
