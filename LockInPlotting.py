import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib
from matplotlib.font_manager import FontProperties
from tkinter import Tk
from matplotlib.ticker import EngFormatter

matplotlib.rcParams['agg.path.chunksize'] = 100000

fontP = FontProperties()
fontP.set_size('small')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
Tk().withdraw()
os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')

filebyfile = 1
expname = 'All'
#filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file
filenames={'/Volumes/lben/lben-commun/2018 User Data/Michael/LockIn/modelcell_05-02-2018_14_40.txt', '/Volumes/lben/lben-commun/2018 User Data/Michael/LockIn/2modelcells, resistor_05-02-2018_14_55.txt'}
range = [10, 10e3]
voltage = 100e-3

for file in filenames:
    filehandle = open(file, 'r')
    (f, V, theta, I, R, C) = np.loadtxt(filehandle, unpack=True)
    V = V*1e-3
    I = I*1e-9
    C = C*1e-9
    theta = -theta
    print(I)
    if range[1]:
        ind = np.argwhere((f >= range[0]) & (f <= range[1])).flatten()
        print(ind)
        f = f[ind]
        V = V[ind]
        theta = theta[ind]
        I = I[ind]
        R = R[ind]
        C = C[ind]
        print(I)

    filehandle.close()
    figIV = plt.figure(1)
    ImpR = voltage/I*np.cos(np.deg2rad(theta))
    ImpIm = voltage/I*np.sin(np.deg2rad(theta))

    axR = figIV.add_subplot(222)
    axR.set_title('Resistance')
    axR.set_xlabel('Frequency')
    axR.set_ylabel('Resistance')
    axR.set_xscale('log')
    axR.set_yscale('log')
    axR.yaxis.set_major_formatter(EngFormatter(unit='$\Omega$', places=0))
    axR.xaxis.set_major_formatter(EngFormatter(unit='Hz', places=0))
    axR.plot(f, R)

    axC = figIV.add_subplot(224)
    axC.set_title('Capacitance')
    axC.set_xlabel('Frequency')
    axC.set_ylabel('Capacitance')
    axC.set_xscale('log')
    axC.set_yscale('linear')
    axC.yaxis.set_major_formatter(EngFormatter(unit='C', places=0))
    axC.xaxis.set_major_formatter(EngFormatter(unit='Hz', places=0))
    axC.plot(f, C)

    axI = figIV.add_subplot(221)
    axI.set_title('Current')
    axI.set_xlabel('Frequency')
    axI.set_ylabel('Current')
    axI.set_xscale('log')
    axI.set_yscale('log')
    axI.yaxis.set_major_formatter(EngFormatter(unit='A', places=0))
    axI.xaxis.set_major_formatter(EngFormatter(unit='Hz', places=0))
    axI.plot(f, I)

    axT = figIV.add_subplot(223)
    axT.set_title('Phase')
    axT.set_xlabel('Frequency')
    axT.set_ylabel('Phase')
    axT.set_xscale('log')
    axT.set_yscale('linear')
    axT.yaxis.set_major_formatter(EngFormatter(unit='°', places=0))
    axT.xaxis.set_major_formatter(EngFormatter(unit='Hz', places=0))
    axT.plot(f, theta)

    figImp = plt.figure(2)
    axNyq = figImp.add_subplot(111)
    axNyq.set_title('Nyquist Plot')
    axNyq.set_xlabel('Real Part [Re(Z)]')
    axNyq.set_ylabel('Imaginary Part [Im(Z)]')
    axNyq.plot(ImpR, ImpIm, ls = 'None', marker = 'o')

    figIV.savefig(file[:-4] + '10kHz.png', dpi=300)
    figImp.savefig(file[:-4] + '_Nyquist10kHz.png', dpi=300)
   # plt.show()
    figIV.clear()
    figImp.clear()

    ## Save Img, Real and Freq to txt file
    np.savetxt(file[:-4] + '_Data10kHz.txt', np.transpose((ImpR, ImpIm, f)), fmt='%.14E', delimiter=' ', newline='\n', header=str(len(ImpIm)), comments='')
