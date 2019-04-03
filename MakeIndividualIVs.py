# -*- coding: utf-8 -*-

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
import argparse
from pprint import pprint
import LoadData

#Set font properties
fontP = FontProperties()
fontP.set_size('small')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

#Set units formatting
from matplotlib.ticker import EngFormatter
Res = EngFormatter(unit='Ω', places=2)
Cond = EngFormatter(unit='S', places=2)
SpesCond = EngFormatter(unit='S/m', places=2)
size = EngFormatter(unit='m', places=2)

Tk().withdraw()

if (platform.system()=='Darwin'):
    os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')

#Define default parameters for size fitting
expname = 'All'
Parameters = {
  'Type': 'Nanopore', #Nanopore, Nanocapillary, NanocapillaryShrunken
  'reversePolarity':  0,
  'specificConductance': 10.5, #10.5 S/m for 1M KCl
  'delay':3, #seconds for reading current
  #Nanopore
  'poreLength' :  1e-9,
  #Nanocapillary
  'taperLength' :  3.3e-3,
  'innerDiameter' : 0.2e-3,
  'taperLengthShaft' : 543e-9,
  'innerDiameterShaft' : 514e-9,

  #Fit option
  'CurveFit' : 'PolyFit', #PolyFit YorkFit

   'fittingMethod' : 'expFit'  #expFit
}

def GetParameters():
    print("Usage:")
    print("run(filenames,newParameters={},verbose=False)")
    print()
    print("Default Parameters:")
    pprint(Parameters)

def run(filenames,newParameters={},verbose=False):
    """
    Function used to call all the necessary other functions to make I-V curves.
    It takes a list of filenames as necessary argument. For each data file, a corresponding
    I-V curve will be generated. 
    
    It iterates on the filenames in the list and calls the corresponding MakeIVData from
    Functions module. 
    
    Parameters
    ----------
    filenames : list of strings
        Full paths to data files.
    newParameters : dict 
        Contains the default parameters for the analysis.
    verbose : bool, optional
        False by default. If True, it will display a simple figure with the shape of the signal.

    """

    for key in newParameters:
        Parameters[key]=newParameters[key]
    Type=Parameters['Type']
    CurveFit=Parameters['CurveFit']
    specificConductance=Parameters['specificConductance']

    for filename in filenames:
        os.chdir(os.path.dirname(filename))
        print(filename)
        #Make Dir to save images
        output = LoadData.OpenFile(filename, verbose=verbose)
        if Parameters['reversePolarity']:
            print('Polarity Reversed!!!!')
            output['i1']= -output['i1']
            output['v1']= -output['v1']

        directory = (str(os.path.split(filename)[0]) + os.sep + expname + '_SavedImages')
        if not os.path.exists(directory):
            os.makedirs(directory)

        AllData = uf.MakeIVData(output, approach=Parameters['fittingMethod'], delay = Parameters['delay'])#, UseVoltageRange = [-0.4, 0.4])
        print('uf.MakeIVData')
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
        Slope=np.abs(AllData[current][CurveFit]['Slope'])
        Yintercept=AllData[current][CurveFit]['Yintercept']

        if Type=='Nanopore':
            poreLength=Parameters['poreLength']
            textstr = 'Nanopore Size\n\nSpecific Conductance: {}\nLength: {}\n\nConductance: {}\nDiameter: {}'\
                .format(SpesCond.format_data(specificConductance),size.format_data(poreLength), Cond.format_data(Slope),
                        size.format_data(uf.CalculatePoreSize(Slope, poreLength, specificConductance)))
        elif Type=='Nanocapillary':
            taperLength=Parameters['taperLength']
            innerDiameter=Parameters['innerDiameter']
            textstr = 'Nanocapillary Size\n\nSpecific Conductance: {}\nTaper lenghth {}:\nInner diameter: {}:\n\nConductance: {}\nDiameter: {}'.\
                format(SpesCond.format_data(specificConductance),size.format_data(taperLength),size.format_data(innerDiameter),Cond.format_data(Slope),
                       size.format_data(uf.CalculateCapillarySize(Slope, innerDiameter, taperLength, specificConductance)))
        elif Type=='NanocapillaryShrunken':
            taperLength=Parameters['taperLength']
            innerDiameter=Parameters['innerDiameter']
            taperLengthShaft=Parameters['taperLengthShaft']
            innerDiameterShaft=Parameters['innerDiameterShaft']
            NCSize=uf.CalculateShrunkenCapillarySize(Slope,innerDiameter, taperLength,specificConductance,taperLengthShaft,innerDiameterShaft)
            textstr = 'Shrunken Nanocapillary Size\n\nSpecific Conductance: {}\nTaper length: {}\nInner diameter: {}\nTaper length at shaft: {}' \
                      '\nInner Diameter at shaft: {}:\n\nConductance: {}\nDiameter: {}'.\
                format(SpesCond.format_data(specificConductance),size.format_data(taperLength),size.format_data(innerDiameter),
                       size.format_data(taperLengthShaft),size.format_data(innerDiameterShaft), Cond.format_data(Slope),
                       size.format_data(NCSize))


        ax1IV.text(0.05, 0.95, textstr, transform=ax1IV.transAxes, fontsize=12,
                  verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))


        ind = np.argsort(AllData[current]['Voltage'])

        ax1IV.errorbar(AllData[current]['Voltage'][ind], AllData[current]['Mean'][ind],
                      yerr=AllData[current]['SE'][ind], fmt='o', color='b')
        ax1IV.plot(AllData[current]['Voltage'][ind], np.polyval([Slope,Yintercept], AllData[current]['Voltage'][ind]), color='r')
        ax1IV.set_title(str(os.path.split(filename)[1])+ '\nR=' + Res.format_data(1/Slope) + ', G=' + Cond.format_data(Slope))
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


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input file')
    parser.add_argument('-c', '--coeff',
                        help='Coefficients for calculating pore size',
                        nargs='+')
    parser
    args = parser.parse_args()

    inputData=args.input
    if inputData==None:
        inputData=askopenfilenames()
    else:
        inputData={inputData}

    newParameters={}
    if args.coeff is not None:
        if len(args.coeff) % 2 == 0:
            for i in range(0, len(args.coeff), 2):
                if i <= len(args.coeff):
                    newParameters[str(args.coeff[i])]=args.coeff[i+1]


    run(inputData,newParameters)
Tk().withdraw()
root=Tk()
root.update()