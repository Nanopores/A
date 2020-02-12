from sys import platform as sys_pf

if sys_pf == 'darwin':
    import matplotlib

    matplotlib.use("TkAgg")

import numpy as np
import pandas as pd
import scipy.signal as sig
import Functions as f
import shelve
from sklearn.mixture import GaussianMixture as GMM
import os
import matplotlib.pyplot as plt
import matplotlib
import h5py
import seaborn as sns

from Plotting import PlotData
# -*- coding: utf-8 -*-



import NanoporeClasses

from Plotting import EventPlots
import importlib
importlib.reload(EventPlots)
from bokeh.io import export_svgs
import EventDetection

TLfolder = r'/Volumes/lben-archives/2019 - CURRENT/Mukesh_Archives/Axopatch200B/30102019/TLfiles'


extension='*.dat'
coeff = {'a': 0.999,
         'E': 0,
         'S': 5,
         'maxEventLength': 50e-3,
         'minEventLength': 100e-6,
         'dt': 25,
         'hbook': 1,
         'delta': 1.2e-9,
         'ChimeraLowPass': 10e3}

Data = EventDetection.run(TLfolder, newExtension=extension, newCoefficients=coeff, verboseLevel=0)

CUSUMevents = Data.GetEventTypes('CUSUM')

CUSUMdata = NanoporeClasses.AllEvents()
CUSUMdata.AddEvent(CUSUMevents)


(p, ph, pv) = EventPlots.PlotGTauVoltageOld(CUSUMdata, voltageLimits = [0.49, 0.51],
                           xLim = [0, 1e-3], yLim = [0, 4e-9],
                           showCurrent = True)





name = 'Test'
plots = [p, ph, pv]

for i, plot in enumerate(plots):
    print(i)
    plot.output_backend = "svg"
    export_svgs(plot, filename=name + str(i) + "_plot_500mV_DevA.svg")