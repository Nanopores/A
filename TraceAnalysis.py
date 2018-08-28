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
#file='30nmPDMA_20180702_152228.log'

folder='//mnt/lben/lben-commun/2018 User Data/Mike/Raw data/7.August/20180809/P23_1MKClboth_Pyrene_TL_200mV'
file='P23_1MKClboth_Pyrene_TL_20180809_005537.log'
#file='P23_1MKClboth_Pyrene_TL_20180809_012433.log'
#folder = askdirectory()

#TranslocationEvents=ED.batcheventdetection(folder,'*.log',coefficients)
#TranslocationEvents.SaveEvents(folder)

test=ED.eventdetection(os.path.join(folder, file), coefficients)
test.PlotAllEvents()

#TranslocationEvents=ED.LoadEvents(folder)


#LengthList=TranslocationEvents.GetAllLengths()
#meanIDrop=TranslocationEvents.GetAllIdrops()
#PlotI_tau(meanIDrop,LengthList)

#Newevents=TranslocationEvents.GetEventsMinCondition(-1,1.5e-11,0.01,[])

#Newevents.PlotAllEvents()
