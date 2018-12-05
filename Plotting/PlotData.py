import Functions
import os
import LoadData
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib.ticker import EngFormatter


Amp = EngFormatter(unit='A', places=2)
Time = EngFormatter(unit='s', places=2)
Volt = EngFormatter(unit='V', places=2)
Cond = EngFormatter(unit='S', places=2)


def SimpleTracePlot(filename, lowPass = 10e3):
    loadedData = LoadData.OpenFile(filename,lowPass,True) #, ChimeraLowPass, True, CutTraces)
    if loadedData['samplerate'] > lowPass:
        output = Functions.LowPass(loadedData['i1'], loadedData['samplerate'], lowPass)
        FullTrace = output['data']
        samplerate = output['samplerate']
    else:
        FullTrace = loadedData['i1']
        samplerate = loadedData['samplerate']

    fig, ax = plt.subplots(figsize=(10, 6))

    times=np.linspace(0, len(FullTrace) / samplerate, num=len(FullTrace))
    ax.plot(times,FullTrace,zorder=1)

    ax.set_xlabel('time (s)')
    ax.set_ylabel('current (A)')
    ax.xaxis.set_major_formatter(Time)
    ax.yaxis.set_major_formatter(Amp)

    plt.title(os.path.basename(filename))

    plt.show()
