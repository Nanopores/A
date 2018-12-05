#Input Stuff
import AnalysisParameters as pm
import numpy as np
import scipy
import scipy.signal as sig
import Functions as f
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
from matplotlib.ticker import EngFormatter
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
import pyqtgraph as pg

fontP.set_size('small')
pm.init()
root = Tk()
root.withdraw()
os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')
SaveAsPDF=0
foldername='PSD'

# Load File, empty string prompts a pop.up window for file selection. Else a file-path can be given
root.update()
filenames = ['/Users/migraf/Desktop/Temp/BigBerta_Noise_2.dat']
#filenames = askopenfilenames()
root.destroy()

# Make Dir to save images
directory = (str(os.path.split(filenames[0])[0]) + os.sep + foldername)
if not os.path.exists(directory):
    os.makedirs(directory)
if SaveAsPDF:
    pp = PdfPages(directory + os.sep + '_ALLPSD.pdf')

for filename in filenames:
    print(filename)
    #filename = '/Users/migraf/Desktop/TestTrace/07B_10mMKClBoth_1kBDN_BothChannels12.dat'
    inp = f.OpenFile(filename)
    folder = str(os.path.split(filename)[0]) + os.sep +'PSD'
    file = os.sep + str(os.path.split(filename)[1][:-4])
    if not os.path.exists(folder):
        os.makedirs(folder)

    fig1 = plt.figure()
    ax = fig1.add_subplot(211)
    ax2 = fig1.add_subplot(212, sharex = ax)
    fr, Pxx_den = scipy.signal.periodogram(inp['i1'], inp['samplerate'])
    #f, Pxx_den = scipy.signal.welch(input, samplerate, nperseg=10*256, scaling='spectrum')
    ax.set_ylabel(r'PSD of Ionic [$\frac{pA^2}{Hz}$]')
    ax.set_xlabel('Frequency')
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.set_ylim(1e4, 1e17)
    #ax.set_xlim(100e-3, 100e3)

    ax.xaxis.set_major_formatter(EngFormatter(unit='Hz'))
    ax.plot(fr, Pxx_den*1e24, 'b')
    ax.grid(1)
    ax.autoscale()

    if inp['graphene']:
        print(np.mean(inp['v2']))
        fr, Pxx_den = scipy.signal.periodogram(inp['i2'], inp['samplerate'])
        ax2.plot(fr, Pxx_den * 1e24, 'r')
        ax2.set_ylabel(r'PSD of Transverse [$\frac{pA^2}{Hz}$]')
        ax2.set_xlabel('Frequency')
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        ax2.xaxis.set_major_formatter(EngFormatter(unit='Hz'))
        ax2.plot(fr, Pxx_den * 1e24, 'r')
        ax2.grid(1)
        ax2.autoscale()
        #ax.legend(['Ion', 'Transverse'])
        textstr = 'STD ionic: {}\nSTD trans: {}'.format(pg.siFormat(np.std(inp['i1'])), pg.siFormat(np.std(inp['i2'])))
        ax2.set_title('STD trans: {}'.format(pg.siFormat(np.std(inp['i2']))))
        ax.set_title('STD ionic: {}'.format(pg.siFormat(np.std(inp['i1']))))
    else:
        textstr = 'STD ionic: {}'.format(pg.siFormat(np.std(inp['i1'])))
        ax.set_title(textstr)
    #ax2.text(0.75, 0.1, textstr, transform=ax.transAxes, fontsize=8,
    #            verticalalignment='bottom', bbox=props)
    fig1.tight_layout()
    if SaveAsPDF:
        pp.savefig(fig1)
    else:
        fig1.savefig(directory + os.sep + str(os.path.split(filename)[1][:-4]) + '_PSD.png', dpi=300)

    fig1.clear()
    ax.clear()
    plt.close(fig1)

if SaveAsPDF:
    pp.close()
