import numpy as np
from pprint import pprint
import EventDetection as ED
import matplotlib.pyplot as plt
from Plotting.EventPlots import PlotI_tau
from tkinter.filedialog import askopenfilenames,askdirectory
import os


savefilename='/mnt/lben-archive/2018 - CURRENT/Jochem/Chimera/2017/2017-12-01/120802B_4MLiCl_5kbpcontrol/200mv_1/Data/Data20180912'
savefilename = '/mnt/lben-archive/2018 - CURRENT/Jochem/Chimera/2018/2018-08-27/NCC3_1MKCl_1/Data/Data20180925'
#fullfilename=os.path.join(folder,file)
#folder = askdirectory()

#TranslocationEvents=ED.batcheventdetection(folder,'*.log',coefficients)  #to re-run ananlysis uncomment these
#TranslocationEvents=ED.eventdetection(fullfilename,coefficients) # To run analysis on a single file uncomment this
#TranslocationEvents.SaveEvents(savefilename) #to re-run ananlysis uncomment these
TranslocationEvents=ED.LoadEvents(savefilename) #to load the already run analysis fles, uncomment this

TranslocationEvents.PlotHistogram()
TranslocationEvents.PlotG_tau()


#
# LengthList=TranslocationEvents.GetAllLengths()
# meanIDrop=TranslocationEvents.GetAllIdrops()
# PlotI_tau(meanIDrop,LengthList)
#
# TranslocationEvents.PlotAllEvents()

# minCurrent=0
# maxCurrent=1 #2.5e-9
# minLength=0 #5e-3
# maxLength=50e-3 #
# newevents=TranslocationEvents.GetEventsMinCondition(minCurrent,maxCurrent,minLength,maxLength)
# newevents.PlotI_tau()
# newevents.PlotAllEvents()
#
# LengthList=newevents.GetAllLengths()
# meanIDrop=newevents.GetAllIdrops()
# PlotI_tau(meanIDrop,LengthList)
#
# newevents.PlotAllEvents()
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
