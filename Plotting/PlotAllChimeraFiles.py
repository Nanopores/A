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

fontP = FontProperties()
fontP.set_size('small')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

Tk().withdraw()
os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')

#filenames = ['/Volumes/backup/2017/Michael/Axopatch/20170512/10mMKClInFlowCellORingPore1mm.dat']
#filenames = ['/Volumes/backup/2017/Michael/Axopatch/20170824/30B_10mMKClCis100mMTrans_PosElecToCis_IV.dat']
expname = 'All'

filenames = askopenfilenames(filetypes = (("chimera files", "*.log"), ("all files", "*.*")))

for filename in filenames:
    print(filename)
    out=uf.ImportChimeraRaw(filename, 200*1e3)
    coefficients = {'a': 0.999, 'E': 0, 'S': 5, 'eventlengthLimit': 10e-3 * out['samplerate']}
    fig1, ax = plt.subplots(1)
    ax.plot(np.arange(0, len(out['i1']))/out['sampleratelp'], out['i1'] * 1e9, 'r')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Current [nA]')
    fig1.savefig(filename[:-4] + '.png', dpi=150)
    fig1.clear()
    plt.close()