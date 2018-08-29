import numpy as np
from pprint import pprint
import EventDetectionJD as ED
import matplotlib.pyplot as plt
from Plotting.NanoporePlots import PlotI_tau
from tkinter.filedialog import askopenfilenames,askdirectory
import os


#Parameters
coefficients = {}
coefficients = {'a': 0.99, 'E': 0, 'S': 4, 'eventlengthLimit': 0.5,'minEventLength': 0}

#folder = '/mnt/lben-archive/2018 - CURRENT/Jochem/Chimera/2018-07-02/30nmPDMA-MTase'
#file='30nmPDMA-MTase_20180702_155437.log'

folder='//mnt/lben/lben-commun/2018 User Data/Mike/Raw data/7.August/20180809/P23_1MKClboth_Pyrene_TL_200mV'
#file='P23_1MKClboth_Pyrene_TL_20180809_005537.log'
#file='P23_1MKClboth_Pyrene_TL_20180809_012433.log'
savefilename='/home/jochem/Documents/test3'

#fullfilename=os.path.join(folder,file)
#folder = askdirectory()

#TranslocationEvents=ED.batcheventdetection(folder,'*.log',coefficients)  #to re-run ananlysis uncomment these
#TranslocationEvents=ED.eventdetection(fullfilename,coefficients) # To run analysis on a single file uncomment this
#TranslocationEvents.SaveEvents(savefilename) #to re-run ananlysis uncomment these
TranslocationEvents=ED.LoadEvents(folder) #to load the already run analysis fles, uncomment this

minCurrent=-400e-12
maxCurrent=400e-12
minLength=0
maxLength=0.015
ewevents=TranslocationEvents.GetEventsMinCondition(minCurrent,maxCurrent,minLength,maxLength)

LengthList=ewevents.GetAllLengths()
meanIDrop=ewevents.GetAllIdrops()
PlotI_tau(meanIDrop,LengthList,True)
#TranslocationEvents.PlotAllEvents()
# #
# minCurrent=-400e-12
# maxCurrent=400e-12
# minLength=0.025
# maxLength=1
# Newevents=TranslocationEvents.GetEventsMinCondition(minCurrent,maxCurrent,minLength,maxLength)
# #
# # LengthList=Newevents.GetAllLengths()
# # meanIDrop=Newevents.GetAllIdrops()
# #PlotI_tau(meanIDrop,LengthList)
#
#
# Newevents.PlotAllEvents()
