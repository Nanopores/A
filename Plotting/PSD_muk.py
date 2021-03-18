# Plots PSD, fits 1/f range and has notch filter to exclude power frequencies

import numpy as np
import scipy
from scipy import signal
import scipy.signal as sig
import scipy.signal as sig
import Functions as uf
import LoadData
import os
import matplotlib
import matplotlib.pyplot as plt
import platform

if platform.system() == 'Darwin':
    matplotlib.use('MacOSX')

from tkinter import Tk
from tkinter.filedialog import askopenfilenames
from matplotlib.font_manager import FontProperties

fontP = FontProperties()
fontP.set_size('small')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib.ticker import EngFormatter

root = Tk()
root.withdraw()

#file = '/Volumes/lben/lben-commun/2018 User Data/Michael/Axopatch/20181003/NIPm10_5nm_1MKCl_pH75_Noise_100kHz_0mV_1.dat'
filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file
#colors=['k','r','b']

root.update()

cutoff_list = [50, 149.9, 249.8, 349.4, 449.37, 450.1, 549.4, 550, 649.3, 748.8, 848.4]
Q_list = [100, 300, 200, 400, 1000, 1000, 2000, 1200, 500, 1200, 2000]

# def func(ln_x, alpha, ln_beta):
#     #return beta * x**(-alpha) # beta is I2A, x is frequency (A is the the low-frequency noise amplitude)
#     return ln_beta - (alpha * ln_x) # Y = beta * x^-alpha, take log of this equation to linear fit and then back to real space

for i,file in enumerate(filenames):
    fig = plt.figure(1, figsize=(10, 8))
    ax = fig.add_subplot(111)
    dat = LoadData.OpenFile(file)
    filename = str(os.path.split(file)[1][:-4])
    os.chdir(os.path.dirname(file))
    directory = (str(os.path.split(file)[0]) + os.sep + 'PSD' + '_SavedImages')

    if not os.path.exists(directory):
        os.makedirs(directory)
    # f, Pxx = sig.welch(dat['i1'], dat['samplerate'], nperseg=2**18)
    # ax.plot(f, Pxx*1e24, label = filename)

    #Notch filtering
    sample_rate = dat['samplerate']
    filtered_data = dat['i1']
    for cutoff, Q in zip(cutoff_list, Q_list):
        b, a = signal.iirnotch(cutoff, Q=Q, fs=sample_rate)
        filtered_data = signal.filtfilt(b, a, filtered_data)

    ff, Pxxf = sig.welch(filtered_data, fs=sample_rate, nperseg=2 ** 18)
    ax.plot(ff, Pxxf * 1e24, label=filename)  # , color = colors[i])
    plt.show()

    # Current Noise RMS calculation
    rms = np.sqrt(np.mean(dat['i1'] ** 2) * 1e12)
    print(rms)


    #1/f curve fitting

    index_restrict = np.where((ff < 100) & (ff>0))
    x_data = ff[index_restrict]
    y_data = Pxxf[index_restrict] * 1e24


    # y_data = y_data[:last_data]
    # popt = np.polyfit(np.log(x_data), np.log(y_data), 1)


    # popt, pcov = scipy.optimize.curve_fit(func, np.log(x_data), np.log(y_data)) # linear fitting with : log Y = log Beta - alpha (log x)
    # y_fit = np.exp(func(np.log(x_data), *popt)) #Y-data for line plot (returns Y from log Y)

    popt = np.polyfit(np.log(x_data), np.log(y_data), 1) #returns sorted according to degree
    polynomial_lin = np.poly1d(popt)
    y_fit = np.exp(polynomial_lin(np.log(x_data)))

    ax.plot(x_data, y_fit, label = "alpha = {}, beta(I2A) = {}".format(popt[0],np.exp(popt[1])))

    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel(r'PSD ($\frac{pA^2}{Hz}$)')
    ax.set_yscale('log')
    ax.set_xscale('log')
    # plt.xlim(10, 1e3)
    ax.legend()

plt.show()

# fig.savefig(directory + os.sep + filename + '_PSD.pdf', transparent=True)
# fig.savefig(directory + os.sep + filename + '_PSD.png', dpi=300)

import numpy as np
import scipy
from scipy import signal
import scipy.signal as sig
import scipy.signal as sig
import Functions as uf
import LoadData
import os
import matplotlib
import matplotlib.pyplot as plt
import platform

if platform.system() == 'Darwin':
    matplotlib.use('MacOSX')

from tkinter import Tk
from tkinter.filedialog import askopenfilenames
from matplotlib.font_manager import FontProperties

fontP = FontProperties()
fontP.set_size('small')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib.ticker import EngFormatter

root = Tk()
root.withdraw()

#file = '/Volumes/lben/lben-commun/2018 User Data/Michael/Axopatch/20181003/NIPm10_5nm_1MKCl_pH75_Noise_100kHz_0mV_1.dat'
filenames = askopenfilenames() # show an "Open" dialog box and return the path to the selected file


root.update()

### For normalizing the traces to zero:

for i,file in enumerate(filenames):
    fig = plt.figure(1, figsize=(10, 8))
    ax = fig.add_subplot(111)
    dat = LoadData.OpenFile(file)

    filename = str(os.path.split(file)[1][:-4])
    os.chdir(os.path.dirname(file))
    directory = (str(os.path.split(file)[0]) + os.sep + 'PSD' + '_SavedImages')

    if not os.path.exists(directory):
        os.makedirs(directory)

    current_array = np.array(dat['i1'] * 1e9)
    current_mean = np.mean(current_array)
    current_array_normalized = current_array - current_mean

    #
    # print (len(current_array))
    # print ((current_mean))
    # print (len(current_array_normalized))

ax.plot(current_array_normalized)

