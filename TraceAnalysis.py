import numpy as np
from pprint import pprint
import EventDetectionJD as ED
import matplotlib.pyplot as plt
from Plotting.NanoporePlots import PlotI_tau
from tkinter.filedialog import askopenfilenames,askdirectory

#Parameters
coefficients = {}
coefficients = {'a': 0.999, 'E': 0, 'S': 4, 'eventlengthLimit': 0.5,'minEventLength': 0}

#folder = '/mnt/lben-archive/2018 - CURRENT/Jochem/Chimera/2018-07-02/30nmPDMA-MTase'
#file='30nmPDMA_20180702_152228.log'

#folder='Z:\lben-commun\2018 User Data\Mike\Raw data\7.August\20180809\P23_1MKClboth_Pyrene_TL_200mV'
folder = askdirectory()

#TranslocationEvents=ED.batcheventdetection(folder,'*.log',coefficients)  #to re-run ananlysis uncomment these
#TranslocationEvents.SaveEvents(folder) #to re-run ananlysis uncomment these

TranslocationEvents=ED.LoadEvents(folder) #to load the already run analysis fles, uncomment this


LengthList=TranslocationEvents.GetAllLengths()
meanIDrop=TranslocationEvents.GetAllIdrops()
PlotI_tau(meanIDrop,LengthList)

minCurrent=-400e-12
maxCurrent=400e-12
minLength=0
maxLength=1
Newevents=TranslocationEvents.GetEventsMinCondition(minCurrent,maxCurrent,minLength,maxLength)

LengthList=Newevents.GetAllLengths()
meanIDrop=Newevents.GetAllIdrops()
PlotI_tau(meanIDrop,LengthList)


Newevents.PlotAllEvents()
