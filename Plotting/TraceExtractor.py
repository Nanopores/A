
### Extracts a part of the trace:

from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")

import NanoporeClasses
import LoadData
import numpy as np
import math
import scipy
from scipy.optimize import curve_fit
from scipy import constants as cst
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import h5py
import pandas as pd
from matplotlib.ticker import EngFormatter
from matplotlib import rc
import LoadData
#rc('text', usetex=True)
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
import pickle as pkl
import scipy.signal as sig
import os
root = Tk()
root.withdraw()


filenames = askopenfilenames()
root.update()

for i,file in enumerate(filenames):
    dat = LoadData.OpenFile(file)


Data = LoadData.OpenFile(file)
samplerate = Data['samplerate']
trace= Data['i1']

x = np.fromfile(file, np.dtype('>f4'))

f = open(file, 'rb')
for i in range(10):
    a=str(f.readline())
    print(a)


end = len(x)
i1 = np.array(x[250:end-1:2])
v1 = np.array(x[251:end:2])

savename2 = 'Extracted_3.dat'


def round_up_to_even(f):
    return math.ceil(f / 2.) * 2

startT  = round_up_to_even(int(188.39*samplerate)*2)  #Specify starting time stamp
endT = round_up_to_even(int(189.61*samplerate)*2)     #Specify ending time stamp
print(startT)
print(endT)
header = x[:250]
data = x[startT:endT]

f2 = open(savename2,'wb')
header.tofile(f2)
data.tofile(f2)
f2.close()

f=open(file, 'rb')
f2 = open(savename2,'rb')
for i in range(11):
    print("    Lines:")
    a=f.readline()
    a2=f2.readline()
    print(str(a))
    print( " hahaha ")
    print(str(a2))

f2.close()
f.close()
print('Done')