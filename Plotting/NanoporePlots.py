﻿import numpy as np
import matplotlib.pyplot as plt





def PlotI_tau(current,tau):
    # definitions for the axes
    left, width = 0.1, 0.6
    bottom, height = 0.1, 0.6
    bottom_h = left_h = left + width + 0.05

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]


    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # the scatter plot:
    axScatter.scatter(tau, current)

    # now determine nice limits by hand:
    extra = 0.1 #0.1 = 10%
    tauRange=np.max(tau)-np.min(tau)
    currentRange=np.max(current)-np.min(current)
    #
    # tauLim = (int(tauMax/binwidth) + 1) * binwidth
    # currentLim = (int(currentMax/binwidth) + 1) * binwidth


    axScatter.set_xlim((np.min(tau)-extra*tauRange, np.max(tau)+extra*tauRange))
    axScatter.set_ylim((np.min(current)-extra*currentRange, np.max(current)+extra*currentRange))

    #tauBins = np.arange(-tauLim, tauLim + binwidth, binwidth)
    #currentBins = np.arange(-currentLim, currentLim + binwidth, binwidth)
    axHistx.hist(tau, bins=10)
    axHisty.hist(current, bins=10, orientation='horizontal')

    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())

    axScatter.set_xlabel('Event length (s)')
    axScatter.set_ylabel('current drop (A)')


    plt.show()

def PlotCurrentTrace(currentTrace,samplerate):
    timeVals=np.linspace(0, len(currentTrace)/samplerate, num=len(currentTrace))
    plt.figure(figsize=(10, 6))

    plt.plot(timeVals,currentTrace)
    plt.xlabel('time (s)')
    plt.ylabel('current (A)')
    plt.show()

def PlotCurrentTraceBaseline(before,currentTrace,after,samplerate):
    timeVals1=np.linspace(0, len(before)/samplerate, num=len(before))
    timeVals2=np.linspace(0+max(timeVals1), len(currentTrace)/samplerate+max(timeVals1), num=len(currentTrace))
    timeVals3=np.linspace(0+max(timeVals2), len(after)/samplerate+max(timeVals2), num=len(after))

    plt.figure(figsize=(10, 6))

    plt.plot(timeVals1,before,color='red')
    plt.plot(timeVals2,currentTrace)
    plt.plot(timeVals3,after,color='red')
    plt.xlabel('time (s)')
    plt.ylabel('current (A)')
    plt.show()