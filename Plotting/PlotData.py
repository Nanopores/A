import Functions
import os
import LoadData
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib.ticker import EngFormatter

from ipywidgets import interact
import numpy as np

from bokeh.io import push_notebook
from bokeh.plotting import figure, show, output_notebook
from bokeh.models import FuncTickFormatter


Amp = EngFormatter(unit='A', places=2)
Time = EngFormatter(unit='s', places=2)
Volt = EngFormatter(unit='V', places=2)
Cond = EngFormatter(unit='S', places=2)

def custom_formatterA():
    units = [
        ('pA', 1e12),
        ('nA', 1e9),
        ('µA', 1e6),
        ('mA', 1e3),
        ('A', 1.0),
    ]
    for u in units:
        if tick <= 1000/u[1]:
            return '{0:.2f}{}'.format(tick * u[1], u[0])

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

def SimpleTracePlot(filename, lowPass = 10e3):
    loadedData = LoadData.OpenFile(filename,lowPass,True) #, ChimeraLowPass, True, CutTraces)
    if loadedData['samplerate'] > lowPass:
        output = Functions.LowPass(loadedData['i1'], loadedData['samplerate'], lowPass)
        FullTrace = output['data']
        samplerate = output['samplerate']
    else:
        FullTrace = loadedData['i1']
        samplerate = loadedData['samplerate']

    output_notebook()

    p = figure(plot_height=300, plot_width=900,tools='pan,box_zoom,xwheel_zoom,reset,save', active_scroll='xwheel_zoom')

    times=np.linspace(0, len(FullTrace) / samplerate, num=len(FullTrace))
    p.line(times,FullTrace)
    #ax.plot(times,FullTrace,zorder=1)

    #ax.set_xlabel('time (s)')
    p.xaxis.axis_label = 'time (s)'
    #ax.set_ylabel('current (A)')
    p.yaxis.axis_label = 'current (A)'
    #ax.xaxis.set_major_formatter(Time)
    #ax.yaxis.set_major_formatter(Amp)

    p.yaxis[0].formatter = FuncTickFormatter.from_py_func(custom_formatterA)
    #p.xaxis[0].formatter = FuncTickFormatter.from_py_func(custom_formattersec)

    p.title.text = os.path.basename(filename)

    show(p)

