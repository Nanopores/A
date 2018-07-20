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

root.update()
#filenames={'/Volumes/lben/lben-commun/2018 User Data/Martina/Axonpatch/20180518_13A/13A_1MKCl_events250mV3.dat'}
filenames = askopenfilenames()
root.destroy()

#filenames=['/Volumes/lben/lben-commun/2018 User Data/Mike/mike 27 03 2018/Events 3/Events 3_20180327_195201.log']

for filename in filenames:
    print(filename)
    inp = f.OpenFile(filename, ChimeraLowPass=pm.ChimeraLowPass)
    folder = pm.OutputFolder
    file = os.sep + str(os.path.split(filename)[1][:-4])
    print('Number of samples in file: {}'.format(len(inp['i1'])))
    if not os.path.exists(folder):
        os.makedirs(folder)

    #Low Pass Event Detection
    AnalysisResults = {}

    chan = 'i1'
    AnalysisResults[chan] = {}
    AnalysisResults[chan]['RoughEventLocations'] = f.RecursiveLowPassFast(inp[chan], pm.coefficients[chan], inp['samplerate'])
    print('Found {} events'.format(len(AnalysisResults[chan]['RoughEventLocations'])))

    ##Event time limit -> deletes events that are too long
    ind = np.array([])
    for (k,i) in enumerate(AnalysisResults[chan]['RoughEventLocations'][:]):
        if i[1] - i[0] <= pm.minmalEventLengthPoints:
            ind = np.append(ind, k)
    AnalysisResults[chan]['RoughEventLocations'] = np.delete(
    AnalysisResults[chan]['RoughEventLocations'][:], ind, axis=0)

    if pm.UpwardsOn: # Upwards detection can be turned on or off
        AnalysisResults[chan + '_Up'] = {}
        AnalysisResults[chan + '_Up']['RoughEventLocations'] = f.RecursiveLowPassFastUp(inp[chan], pm.coefficients[chan], inp['samplerate'])
        ind = np.array([])
        for (k, i) in enumerate(AnalysisResults[chan + '_Up']['RoughEventLocations'][:]):
            if i[1] - i[0] <= pm.minmalEventLengthPoints:
                ind = np.append(ind, k)
            AnalysisResults[chan + '_Up']['RoughEventLocations'] = np.delete(
                AnalysisResults[chan + '_Up']['RoughEventLocations'][:], ind, axis=0)

    print('Deleted All Events shorter than {:0.2e}s'.format(pm.minmalEventLengthPoints*inp['samplerate']))


    ############Plot the Lowpass Detections
    if pm.PlotTheLowPassDetection:
        fig = plt.figure(1, figsize=(16,5) )
        ax = fig.add_subplot(111)

        ax.plot(np.arange(0, len(inp['i1']), 1)/inp['samplerate'], inp['i1'], 'b')

        for i in AnalysisResults['i1']['RoughEventLocations']:
                ax.plot(np.arange(np.uint64(i[0]), np.uint64(i[1]), 1)/inp['samplerate'], inp['i1'][np.uint64(i[0]):np.uint64(i[1])], 'r')

        if pm.UpwardsOn:
            for i in AnalysisResults['i1_Up']['RoughEventLocations']:
                    ax.plot(y=inp['i1'][np.uint64(i[0]):np.uint64(i[1])], x=np.arange(np.uint64(i[0]), np.uint64(i[1]), 1)/inp['samplerate'], pen='r')

        ax.set
        ax.set_ylabel("Current [A]")
        ax.set_xlabel("Time [s]")



    # Refine the Rough Event Detection done by the LP filter and Add event infos
    if inp['graphene']:
        AnalysisResults = f.RefinedEventDetection(inp, AnalysisResults, signals=['i1', 'i2'], limit=pm.MinimalFittingLimit*inp['samplerate'])
    else:
        AnalysisResults = f.RefinedEventDetection(inp, AnalysisResults, signals=['i1'], limit=pm.MinimalFittingLimit*inp['samplerate'])

    ## Print Lenghts
    f.SaveToHDF5(inp, AnalysisResults, pm.coefficients, folder)

    if pm.PlotTheLowPassDetection:
        plt.show()