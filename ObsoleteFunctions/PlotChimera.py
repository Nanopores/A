#Input Stuff
import MiscParameters as pm
import numpy as np
import Functions as f
import os
import matplotlib.pyplot as plt
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
from tkinter.filedialog import askopenfilename
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
pm.init(LoadFiles = 0)
root = Tk()
root.withdraw()
os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')

# Load File, empty string prompts a pop.up window for file selection. Else a file-path can be given

LP = 10e3

root.update()
filenames = askopenfilenames()
root.destroy()

for filename in filenames:
    print(filename)
    inp = f.OpenFile(filename, ChimeraLowPass=LP, approxImpulseResponse=True)
    folder = pm.OutputFolder
    file = os.sep + str(os.path.split(filename)[1][:-4])
    print('Number of samples in file: {}'.format(len(inp['i1'])))
    if not os.path.exists(folder):
        os.makedirs(folder)
    fig = plt.figure(1, figsize=(16, 5))
    ax = fig.add_subplot(111)
    ax.plot(np.arange(0, len(inp['i1'][100:]))/inp['samplerate'], inp['i1'][100:]*1e9)
    ax.set_ylabel('Current (nA)')
    ax.set_xlabel('Time (s)')
    plt.show()
    #fig.savefig('/Volumes/lben/lben-commun/2018 User Data/Publications/Ke-ssdsDNA/raw materials/raw chimera traces/'+ file+'.pdf', transparent=True)
    ax.clear()
    fig.clear()