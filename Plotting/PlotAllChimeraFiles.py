from ObsoleteFunctions import UsefulFunctions as uf
import matplotlib.pyplot as plt
import numpy as np
import os
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

LowPassCutOFF = 50 #kHz

filenames = askopenfilenames(filetypes = (("chimera files", "*.log"), ("all files", "*.*")))

for filename in filenames:
    print(filename)
    out=uf.ImportChimeraRaw(filename, LowPassCutOFF*1e3)
    coefficients = {'a': 0.999, 'E': 0, 'S': 5, 'maxEventLength': 10e-3 * out['samplerate']}
    fig1, ax = plt.subplots(1)
    ax.plot(np.arange(0, len(out['i1']))/out['sampleratelp'], out['i1'] * 1e9, 'g')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Current [nA]')
    plt.show()
    #fig1.savefig(filename[:-4] + '.png', dpi=150)
    #fig1.clear()