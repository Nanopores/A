import UsefulFunctions as uf
import pyqtgraph as pg
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import scipy
import numpy as np
import os
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
from matplotlib.ticker import EngFormatter
matplotlib.rcParams['agg.path.chunksize'] = 100000

fontP = FontProperties()
fontP.set_size('small')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
Tk().withdraw()
os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')

#filenames = ['/Volumes/backup/2017/Michael/Axopatch/20170512/10mMKClInFlowCellORingPore1mm.dat']
filenames = ['/Volumes/lben-archives/2016/Michael/Axopatch/21112016/17B_100mMCis1MtransKCl_80mer_9.dat']
expname = 'All'

#filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file
#filenames={'/Users/migraf/Desktop/22B_UpperRibbonOverNight_.dat'}

#With Voltages
volt = 1

for filename in filenames:
    print(filename)
    out = uf.ImportAxopatchData(filename)
    fig = plt.figure(1, figsize=(8.7,5.96) )

    if volt:
        ax = fig.add_subplot(311)
        ax2 = fig.add_subplot(312, sharex = ax)
        ax3 = fig.add_subplot(313, sharex = ax)
        ax3.yaxis.set_major_formatter(EngFormatter(unit='V'))
        ax3.set_ylabel('Voltage')
    else:
        ax = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex = ax)

    ax.plot(np.arange(0, len(out['i1']))/out['samplerate'], out['i1'], 'b')
    ax2.plot(np.arange(0, len(out['i2']))/out['samplerate'], out['i2'], 'r')

    if volt:
        ax3.plot(np.arange(0, len(out['i1']))/out['samplerate'], out['v1'], 'b')
        ax3.plot(np.arange(0, len(out['i1'])) / out['samplerate'], out['v2'], 'r')
    #ax.set_xlim(41000, 42600)
    #ax2.set_xlim(41000, 42600)
    #ax.set_ylim(-5,15)

    ax.set_xlabel('Time')
    ax2.set_xlabel('Time')

    ax.set_ylabel('Ionic Current')
    ax2.set_ylabel('Transverse Current')

    ax.xaxis.set_major_formatter(EngFormatter(unit='s'))
    ax.yaxis.set_major_formatter(EngFormatter(unit='A'))
    ax2.yaxis.set_major_formatter(EngFormatter(unit='A'))


    #fig1.savefig(filename[:-4] + '.png', dpi=300)
    plt.show()
    #fig1.clear()
