import numpy as np
import scipy
import scipy.signal as sig
import Functions as uf
import pyqtgraph as pg
import os
import matplotlib.pyplot as plt
import matplotlib
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
from matplotlib.font_manager import FontProperties
import platform
import csv

fontP = FontProperties()
fontP.set_size('small')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib.ticker import EngFormatter
Res = EngFormatter(unit='Ω', places=2)
Cond = EngFormatter(unit='S', places=2)
SpesCond = EngFormatter(unit='S/m', places=2)
size = EngFormatter(unit='m', places=2)

Tk().withdraw()

if (platform.system()=='Darwin'):
    os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')

expname = 'All'
Type='Nanocapillary' #Nanopore, Nanocapillary, NanocapillaryShrunken
reversePolarity = 0
specificConductance=10.5 #10.5 S/m for 1M KCl

delay=0.5 #seconds for reading current

#Nanopore
poreLength =  1e-9

#Nanocapillary
taperLength =  3.3e-3
innerDiameter = 0.2e-3
taperLengthShaft = 543e-9
innerDiameterShaft = 514e-9

filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file
#filenames={'/mnt/lben-archive/2018 - CURRENT/Jochem/Chimera/2018/2018-08-27/NCC3_1MKCl_1/IV/IV_NCC_1MKCl_1_20180827_084204.log'}

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

    AllData = uf.MakeIVData(output, delay = delay)#, UseVoltageRange = [-0.4, 0.4])
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
        figIV2 = plt.figure(3, figsize=(10, 10))
        figIV2.clear()
        ax2IV = figIV2.add_subplot(111)
        ax2IV = uf.PlotIV(output, AllData, current='i2', unit=1e9, axis=ax2IV, WithFit=1)
        figIV2.tight_layout()
        figIV2.savefig(directory + os.sep + str(os.path.split(filename)[1]) + '_IV_i2.png', dpi=300)
        figIV2.savefig(directory + os.sep + str(os.path.split(filename)[1]) + '_IV_i2.eps')

    figIV = plt.figure(2, figsize=(10, 7))
    ax1IV = figIV.add_subplot(111)
    current = 'i1'

    #ax1IV = uf.PlotIV(output, AllData, current='i1', unit=1, axis = ax1IV, WithFit = 1, useEXP = 0, color ='y',
    #                labeltxt='MeanFit', PoreSize=[10, 1e-9], title=str(os.path.split(filename)[1]))
    if Type=='Nanopore':
        textstr = 'Nanopore Size\n\nSpecific Conductance: {}\nLenght: {}\n\nConductance: {}\nDiameter: {}'\
            .format(SpesCond.format_data(specificConductance),size.format_data(poreLength), Cond.format_data(AllData[current]['YorkFitValues']['Slope']),
                    size.format_data(uf.CalculatePoreSize(AllData[current]['YorkFitValues']['Slope'], poreLength, specificConductance)))
    elif Type=='Nanocapillary':
        textstr = 'Nanocapillary Size\n\nSpecific Conductance: {}\nTaper lenght: {}:\nInner diameter: {}:\n\nConductance: {}\nDiameter: {}'.\
            format(SpesCond.format_data(specificConductance),size.format_data(taperLength),size.format_data(innerDiameter),Cond.format_data(AllData[current]['YorkFitValues']['Slope']),
                   size.format_data(uf.CalculateCapillarySize(AllData[current]['YorkFitValues']['Slope'], innerDiameter, taperLength, specificConductance)))
    elif Type=='NanocapillaryShrunken':
        NCSize=uf.CalculateShrunkenCapillarySize(AllData[current]['YorkFitValues']['Slope'],innerDiameter, taperLength,specificConductance,taperLengthShaft,innerDiameterShaft)
        textstr = 'Shrunken Nanocapillary Size\n\nSpecific Conductance: {}\nTaper lenght: {}\nInner diameter: {}\nTaper length at shaft: {}' \
                  '\nInner Diameter at shaft: {}:\n\nConductance: {}\nDiameter: {}'.\
            format(SpesCond.format_data(specificConductance),size.format_data(taperLength),size.format_data(innerDiameter),
                   size.format_data(taperLengthShaft),size.format_data(innerDiameterShaft), Cond.format_data(AllData[current]['YorkFitValues']['Slope']),
                   size.format_data(NCSize))


    ax1IV.text(0.05, 0.95, textstr, transform=ax1IV.transAxes, fontsize=12,
              verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ind = np.argsort(AllData[current]['Voltage'])
    p = np.polyfit(AllData[current]['Voltage'][ind], AllData[current]['Mean'][ind], 1)
    ax1IV.errorbar(AllData[current]['Voltage'][ind], AllData[current]['Mean'][ind],
                  yerr=AllData[current]['STD'][ind], fmt='o', color='b')
    ax1IV.plot(AllData[current]['Voltage'][ind], np.polyval(p, AllData[current]['Voltage'][ind]), color='r')
    ax1IV.set_title(str(os.path.split(filename)[1])+ '\nR=' + Res.format_data(1/p[0]) + ', G=' + Cond.format_data(p[0]))
    ax1IV.set_ylabel('Current')
    ax1IV.set_xlabel('Voltage')
    ax1IV.xaxis.set_major_formatter(EngFormatter(unit='V'))
    ax1IV.yaxis.set_major_formatter(EngFormatter(unit='A'))



    figIV.savefig(directory + os.sep + str(os.path.split(filename)[1]) + Type+'IV_i1.pdf', transparent=True)



    x=AllData[current]['Voltage'][ind]
    y=AllData[current]['Mean'][ind]

    csvfile=directory + os.sep + str(os.path.split(filename)[1]) + Type+'IV_i1.csv'
    with open(csvfile, 'w') as output:
        writer=csv.writer(output, lineterminator='\n')
        for i in range(len(x)):
            writer.writerow([x[i] , y[i]])

    plt.show()
    figIV.clear()