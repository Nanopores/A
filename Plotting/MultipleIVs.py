# -*- coding: utf-8 -*-
from sys import platform as sys_pf
import matplotlib
from sys import platform as sys_pf

if sys_pf == 'darwin':
    import matplotlib

    matplotlib.use("TKAgg")
import numpy as np
import Functions as uf
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

# Set font properties
fontP = FontProperties()
fontP.set_size('small')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Set units formatting
from matplotlib.ticker import EngFormatter

Res = EngFormatter(unit='Ω', places=2)
Cond = EngFormatter(unit='S', places=2)
SpesCond = EngFormatter(unit='S/m', places=2)
size = EngFormatter(unit='m', places=2)

if (platform.system() == 'Darwin'):
    Tk().withdraw()
    os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')

# Define default parameters for size fitting

expname = 'All'
Parameters = {
    'Type': 'Nanopore',  # Nanopore, Nanocapillary, NanocapillaryShrunken
    'reversePolarity': 0,
    'specificConductance': 21,  # 10.5 S/m for 1M KCl
    'delay': 5,  # seconds for reading current
    # Nanopore
    'poreLength': 1e-9,

    # Fit option
    'CurveFit': 'PolyFit',  # PolyFit YorkFit
    'fittingMethod': 'mean'  # expFit, mean
}


def GetParameters():
    print("Usage:")
    print("run(filenames,newParameters={},verbose=False, noPlot=False")
    print('Returns poresize in m')
    print()
    print("Default Parameters:")
    print(Parameters)


def run(filenames, newParameters={}, verbose=False, noPlot=False):

    for key in newParameters:
        Parameters[key] = newParameters[key]
    Type = Parameters['Type']
    CurveFit = Parameters['CurveFit']
    specificConductance = Parameters['specificConductance']

    for filename in filenames:
        os.chdir(os.path.dirname(filename))
        print(filename)
        # Make Dir to save images
        output = LoadData.OpenFile(filename, verbose=verbose)
        if Parameters['reversePolarity']:
            print('Polarity Reversed!!!!')
            output['i1'] = -output['i1']
            output['v1'] = -output['v1']

        directory = (str(os.path.split(filename)[0]) + os.sep + expname + '_SavedImages')
        if not os.path.exists(directory):
            os.makedirs(directory)

        AllData = uf.MakeIVData(output, approach=Parameters['fittingMethod'],
                                delay=Parameters['delay'])  # , UseVoltageRange = [-0.51, 0.51])
        if AllData == 0:
            print('!!!! No Sweep in: ' + filename)
            continue

        if output['graphene']:
            currents = ['i1', 'i2']
        else:
            currents = ['i1']

        for current in currents:
            Slope = AllData[current][CurveFit]['Slope']
            Yintercept = AllData[current][CurveFit]['Yintercept']

            if Type == 'Nanopore':
                poreLength = Parameters['poreLength']
                poresize = uf.CalculatePoreSize(np.abs(Slope), poreLength, specificConductance)
            if not noPlot:
                figIV = plt.figure(2, figsize=(10, 7))
                ax1IV = figIV.add_subplot(111)
                ind = np.argsort(AllData[current]['Voltage'])

        ax1IV.errorbar(AllData[current]['Voltage'][ind], AllData[current]['Mean'][ind],
                       yerr=AllData[current]['SE'][ind], fmt='o')
        ax1IV.plot(AllData[current]['Voltage'][ind], np.polyval([Slope, Yintercept], AllData[current]['Voltage'][ind]), label=str(os.path.split(filename)[1])+ '\nR=' + Res.format_data(1/Slope) + ', G=' + Cond.format_data(Slope))
        ax1IV.set_ylabel('Current ' + current)
        ax1IV.set_xlabel('Voltage')
        ax1IV.legend()
        ax1IV.xaxis.set_major_formatter(EngFormatter(unit='V'))
        ax1IV.yaxis.set_major_formatter(EngFormatter(unit='A'))
        # ax1IV.set_xlim(-0.6, 0.6)
        # ax1IV.set_ylim(-15e-9, 15e-9)
        ax1IV.grid(True)
        figIV.savefig(directory + os.sep + str(os.path.split(filename)[1]) + Type + '_' + current + 'IV.pdf',
                      transparent=True)
        plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input file')
    parser.add_argument('-c', '--coeff',
                        help='Coefficients for calculating pore size',
                        nargs='+')
    parser
    args = parser.parse_args()

    inputData = args.input
    if inputData == None:
        inputData = askopenfilenames()
    else:
        inputData = {inputData}

    newParameters = {}
    if args.coeff is not None:
        if len(args.coeff) % 2 == 0:
            for i in range(0, len(args.coeff), 2):
                if i <= len(args.coeff):
                    newParameters[str(args.coeff[i])] = args.coeff[i + 1]

    run(inputData, newParameters)

if (platform.system() == 'Darwin'):
    Tk().withdraw()
    root = Tk()
    root.destroy()