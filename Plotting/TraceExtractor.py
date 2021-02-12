
### Extracts a part of the trace:

from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import os
import numpy as np
import math
import LoadData
#rc('text', usetex=True)
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
root = Tk()
root.withdraw()

filenames = askopenfilenames()
root.update()

for i,file in enumerate(filenames):
    dat = LoadData.OpenFile(file)

Data = LoadData.OpenFile(file)
samplerate = Data['samplerate']
trace = Data['i1']

x = np.fromfile(file, np.dtype('>f4'))

f = open(file, 'rb')  #rb : Opens the file as read-only in binary format and starts reading from the beginning of the file.
for i in range(10): #range() function. To loop through a set of code a specified number of times, we can use the range() function, The range() function returns a sequence of numbers, starting from 0 by default, and increments by 1 (by default), and ends at a specified number.
    a = str(f.readline())
    print(a)

end = len(x)
i1 = np.array(x[250:end-1:2])
v1 = np.array(x[251:end:2])

savename2 = 'multiple_trappings_extracted_.dat'

#Save options
os.chdir(os.path.dirname(file))
directory = (str(os.path.split(file)[0]) + os.sep + 'Extracted' + '_traces')
if not os.path.exists(directory):
    os.makedirs(directory)

def round_up_to_even(f):
    return math.ceil(f / 2.) * 2

startT = round_up_to_even (int(134.2*samplerate)*2)  #Specify starting time stamp
endT = round_up_to_even (int(260*samplerate)*2)     #Specify ending time stamp
print(startT)
print(endT)

header = x[:250]
data = x[startT:endT]

f2 = open(savename2,'wb') # wb indicates that the file is opened for writing only in binary mode
header.tofile(f2)
data.tofile(f2)
f2.close()

f = open(file, 'rb') # Opens a file for reading only in binary format
f2 = open(savename2,'rb')

for i in range(11):
    print("    Lines:")
    a = f.readline()
    a2 = f2.readline()
    print(str(a))
    print( " hahaha ")
    print(str(a2))

f2.close()
f.close()
print('Done')