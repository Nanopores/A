import Functions
import os
import LoadData
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy import signal

from matplotlib.ticker import EngFormatter

from ipywidgets import interact
import numpy as np

from bokeh.io import push_notebook
from bokeh.plotting import figure, show, output_notebook
from bokeh.models import FuncTickFormatter, NumeralTickFormatter

# Needed for plotting in jupyter
from bokeh.resources import INLINE
import bokeh.io
bokeh.io.output_notebook(INLINE)

import math

Amp = EngFormatter(unit='A', places=2)
Time = EngFormatter(unit='s', places=2)
Volt = EngFormatter(unit='V', places=2)
Cond = EngFormatter(unit='S', places=2)


def getdata(filename,lowPass):
    if isinstance(filename, str):
        data = LoadData.OpenFile(filename,lowPass,True) #, ChimeraLowPass, True, CutTraces)
    elif isinstance(filename,dict):
        data = filename
        filename = data['filename']
    else:
        raise Exception('Incorrect input')
    return data

def custom_formatterA():
    units = [('m', 1e-3),
             ('µ', 1e-6),
             ('n', 1e-9),
             ('p', 1e-12)]
    if tick == 0:
        return '0'
    for u in units:
        if Math.abs(tick) >= u[1]:
            return '{0:.1f} {}A'.format(tick / u[1], u[0]) #Note, Javascript Math
    return '{:.2E}'.format(tick)


def custom_formattersec():
    units = [
        ('µs', 1e6),
        ('ms', 1e3),
        ('s', 1.0),
        ('min', 1/60),
        ('hour', 1/3600),
    ]
    for u in units:
        if tick <= 100/u[1]:
            return '{0:.1f}{}'.format(tick * u[1], u[0])


def SimpleTracePlot(filename, lowPass=10e3):
    loadedData = getdata(filename, lowPass)

    if loadedData['samplerate'] > lowPass:
        output = Functions.LowPass(loadedData['i1'], loadedData['samplerate'], lowPass)
        FullTrace = output['data']
        samplerate = output['samplerate']
    else:
        FullTrace = loadedData['i1']
        samplerate = loadedData['samplerate']

    output_notebook()

    p = figure(plot_height=300, plot_width=900,tools='pan,box_zoom,xwheel_zoom,reset,save', active_scroll='xwheel_zoom')

    times = np.linspace(0, len(FullTrace) / samplerate, num=len(FullTrace))
    p.line(times,FullTrace)

    p.xaxis.axis_label = 'time (s)'
    p.yaxis.axis_label = 'current (A)'

    p.yaxis[0].formatter = FuncTickFormatter.from_py_func(custom_formatterA)

    p.title.text = os.path.basename(filename)

    show(p)

def PlotPSD(inputdata):
    if isinstance(inputdata, str):
        filename = inputdata
        loadedData = LoadData.OpenFile(filename, approxImpulseResponse=True) #, ChimeraLowPass, True, CutTraces)
    else:
        filename = ''
        loadedData = inputdata
    frequencies, P_den = Functions.GetPSD(loadedData)

    output_notebook()

    p = figure(plot_height=300, plot_width=900, x_axis_type="log", y_axis_type="log", tools='pan,box_zoom,xwheel_zoom,reset,save')

    p.line(frequencies,P_den*1e24)

    p.xaxis.axis_label = 'Frequencies (Hz)'
    p.yaxis.axis_label = 'PSD pA^2/Hz)'

    p.title.text = 'PSD :' + os.path.basename(filename)

    show(p)
