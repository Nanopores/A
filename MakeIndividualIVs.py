import numpy as np
import scipy
import scipy.signal as sig
import Functions as uf
import pyqtgraph as pg
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

Tk().withdraw()
os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')

#filenames = ['/Volumes/backup/2017/Michael/Axopatch/20170512/10mMKClInFlowCellORingPore1mm.dat']
#filenames = ['/Volumes/Michael/Axopatch/20180430/R17_ph74_KCl_1M_1M_0mW_473nm_IV_3.dat']
expname = 'All'
reversePolarity = 0

filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file

for filename in filenames:
    print(filename)

    #Make Dir to save images
    output = uf.OpenFile(filename)
    if reversePolarity:
        print('Polarity Reversed!!!!')
        output['i1']= -output['i1']
        output['v1']= -output['v1']

    directory = (str(os.path.split(filename)[0]) + os.sep + expname + '_SavedImages')
    if not os.path.exists(directory):
        os.makedirs(directory)

    AllData = uf.MakeIVData(output, delay = 9)#, UseVoltageRange = [-0.4, 0.4])
    if AllData == 0:
        print('!!!! No Sweep in: ' + filename)
        continue

    #Plot Considered Part
    #figExtracteParts = plt.figure(1)
    #ax1 = figExtracteParts.add_subplot(211)
    #ax2 = figExtracteParts.add_subplot(212, sharex=ax1)
    #(ax1, ax2) = uf.PlotExtractedPart(output, AllData, current = 'i1', unit=1e9, axis = ax1, axis2=ax2)
    #plt.show()
    #figExtracteParts.savefig(directory + os.sep + 'PlotExtracted_' + str(os.path.split(filename)[1])[:-4] + '.eps')
    #figExtracteParts.savefig(directory + os.sep + 'PlotExtracted_' + str(os.path.split(filename)[1])[:-4] + '.png', dpi=150)


    # Plot IV
    if output['graphene']:
        figIV2 = plt.figure(3)
        figIV2.clear()
        ax2IV = figIV2.add_subplot(111)
        ax2IV = uf.PlotIV(output, AllData, current='i2', unit=1e9, axis=ax2IV, WithFit=1)
        figIV2.tight_layout()
        figIV2.savefig(directory + os.sep + str(os.path.split(filename)[1]) + '_IV_i2.png', dpi=300)
        figIV2.savefig(directory + os.sep + str(os.path.split(filename)[1]) + '_IV_i2.eps')

    figIV = plt.figure(2)
    ax1IV = figIV.add_subplot(111)
    #ax1IV = uf.PlotIV(output, AllData, current='i1', unit=1e9, axis = ax1IV, WithFit = 1, PoreSize=[10, 0.7e-9], useEXP = 0, color ='b',
    #                  labeltxt='Exponential Fit', Polynomial=0)
    ax1IV = uf.PlotIV(output, AllData, current='i1', unit=1e9, axis = ax1IV, WithFit = 1, useEXP = 0, color ='y',
                      labeltxt='MeanFit', PoreSize=[10, 1e-9], title=str(os.path.split(filename)[1]))

    ax1IV.legend(loc='upper center', bbox_to_anchor=(0.8, 0.2),
                 ncol=1, fancybox=True, shadow=True, prop=fontP)

    figIV.tight_layout()
    #Save Figures
    #plt.show()
    figIV.savefig(directory + os.sep + str(os.path.split(filename)[1]) + 'IV_i1.png', dpi=300)
    figIV.savefig(directory + os.sep + str(os.path.split(filename)[1]) + 'IV_i1.eps')

    #plt.show()
    figIV.clear()