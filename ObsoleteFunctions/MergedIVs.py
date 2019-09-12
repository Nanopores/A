
import UsefulFunctions as uf
import os

import matplotlib.pyplot as plt
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

Tk().withdraw()
os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')
#
#filenames = ['/Volumes/backup/2017/Michael/Axopatch/20170512/10mMKClInFlowCellORingPore1mm.dat']
#filenames = ['/Volumes/Michael/Axopatch_DR/20171025/03A_GroundOnTrans_OvernightAfterChimeraMeasurements_IV.dat','/Volumes/Michael/Axopatch_DR/20171025/03A_AfterWash_1MKClpH7_GroundToTrans_IV_2.dat','/Volumes/Michael/Axopatch_DR/20171025/03A_AfterWash_1MKClpH115_GroundToTrans_IV_2.dat','/Volumes/Michael/Axopatch_DR/20171025/03A_AfterWash_1MKClpH7_GroundToTrans_500mV_IV_6.dat','/Volumes/Michael/Axopatch_DR/20171025/03A_AfterWash_1MKClpH341_GroundToTrans_500mV_IV_3.dat', '/Volumes/Michael/Axopatch_DR/20171025/03A_AfterWash_1MKClpH7_GroundToTrans_500mV_IV_7.dat']

#labels = ['1) After Poly-D-Lysine, pH7', '2) 1M KCl, pH 7','3) 1M KCl, pH 11', '4) 1M KCl, pH 7','5) 1M KCl, pH 3''6) 1M KCl, pH 7']
#
expname = 'All'

filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file

figIV = plt.figure(2)
ax1IV = figIV.add_subplot(111)

for filename in filenames:
    print(filename)

    #Make Dir to save images
    output = uf.OpenFile(filename)
    directory = (str(os.path.split(filename)[0]) + os.sep + expname + '_SavedImages')
    if not os.path.exists(directory):
        os.makedirs(directory)

    AllData = uf.MakeIVData(output, delay = 2)
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
        ax2IV = figIV2.add_subplot(111)
        ax2IV = uf.PlotIV(output, AllData, current='i2', unit=1e9, axis=ax2IV, WithFit=1)
        figIV2.tight_layout()
        figIV2.savefig(directory + os.sep + str(os.path.split(filename)[1]) + '_IV_i2.png', dpi=150)
        figIV2.savefig(directory + os.sep + str(os.path.split(filename)[1]) + '_IV_i2.eps')
        figIV2.clear()

    ax1IV = uf.PlotIV(output, AllData, current = 'i1', unit=1e9, axis = ax1IV, WithFit = 0)
    figIV.tight_layout()

    #Save Figures
ax1IV.legend()
figIV.savefig(directory + os.sep + str(os.path.split(filename)[1]) + 'Merged_IV_i1.png', dpi=300)
figIV.savefig(directory + os.sep + str(os.path.split(filename)[1]) + 'Merged_IV_i1.eps')
figIV.clear()


# print(AllData['i2']['YorkFitValues']['Yintercept'])

#plt.show()