import numpy as np
from matplotlib.ticker import FuncFormatter
from scipy import constants as cst
import Functions as uf
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

path ='/Users/migraf/SWITCHdrive/PhD/Comsol Conductance Model/OsmoVsConc'
xls = pd.ExcelFile('/Users/migraf/SWITCHdrive/PhD/Comsol Conductance Model/OsmoVsConc/AllData.xlsx')
data = xls.parse(0)
ThingsToPlot = data.keys()

whichX = 3
whichLoop = 1
whichConst = 0


poresizes = np.unique(data[ThingsToPlot[0]])
SC = np.unique(data[ThingsToPlot[1]])
Volt = np.unique(data[ThingsToPlot[2]])
c_cis = np.unique(data[ThingsToPlot[3]])

## Several Surface Charges
SC = [-0.01, -50, -100]
fig2 = plt.figure(2, figsize=(10, 10))
ax2 = fig2.add_subplot(2, 1, 1)
ax3 = fig2.add_subplot(2, 1, 2)

var={}
## Make All IVs
var['Poresize'] = np.array([])
var['SurfaceCharge'] = np.array([])
var['C_Cis'] = np.array([])
var['OsmoticCurrent'] = np.array([])
var['OsmoticVoltage'] = np.array([])
var['IVData'] = np.array([])
for po in poresizes:
    for sc in SC:
        yVolt = []
        yCurr = []
        for cisc in c_cis:
            ind1 = np.argwhere(np.isclose(data[ThingsToPlot[0]].values, po))
            ind2 = np.intersect1d(ind1, np.argwhere(np.isclose(data[ThingsToPlot[1]].values, sc)))
            ind = np.intersect1d(ind2, np.argwhere(np.isclose(data[ThingsToPlot[3]].values, cisc)))
            x = -data[ThingsToPlot[2]].values[ind]*1e-3
            y = -(data[ThingsToPlot[4]].values[ind]-data[ThingsToPlot[5]].values[ind]) * cst.N_A * cst.e
            p = np.polyfit(x, y, 1)
            var['Poresize'] = np.append(var['Poresize'], po)
            var['SurfaceCharge'] = np.append(var['SurfaceCharge'], sc)
            var['C_Cis'] = np.append(var['C_Cis'], cisc)
            var['OsmoticCurrent'] = np.append(var['OsmoticCurrent'], p[1])
            var['OsmoticVoltage'] = np.append(var['OsmoticVoltage'], p[1] / p[0])
            var['IVData'] = np.append(var['IVData'], [x, y])
            yVolt.append(p[1] / p[0])
            yCurr.append(p[1])
        ax2.plot(1000/c_cis, yVolt, ls='--', marker = 'o', linewidth=2, markersize=10, label =EngFormatter(unit=r'$\frac{C}{m^2}$', places=0).format_data(sc*1e-3))
        ax3.plot(1000/c_cis, yCurr, ls='--', marker = 'o', linewidth=2, markersize=10, label =EngFormatter(unit=r'$\frac{C}{m^2}$', places=0).format_data(sc*1e-3))

    ax2.xaxis.set_major_formatter(EngFormatter(unit=r'$\frac{c_{max}}{c_{min}}$'))  # r'$\frac{K^+}{Cl^-}$'))
    ax2.yaxis.set_major_formatter(EngFormatter(unit='V'))
    ax3.yaxis.set_major_formatter(EngFormatter(unit='A'))
    ax3.xaxis.set_major_formatter(EngFormatter(unit=r'$\frac{c_{max}}{c_{min}}$'))  # r'$\frac{K^+}{Cl^-}$'))

    ax2.set_ylabel('Osmotic Voltage')
    ax3.set_ylabel('Osmotic Current')
    ax2.set_xlabel('Concentration Gradient')
    ax3.set_xlabel('Concentration Gradient')
    ax2.set_xscale('log')
    ax3.set_xscale('log')
    ax2.legend()
    ax2.set_title('Pore Size: {}nm'.format(np.round(po)))
    fig2.savefig(path + os.sep + 'OsmoticVsGradientFor{}nm.pdf'.format(int(np.round(po))), transparent=True)
    ax2.clear()
    ax3.clear()