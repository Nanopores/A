import pandas as pd
import datetime
import os
from tkinter.filedialog import askopenfilenames
import numpy as np
title = '20180711_RW18_Temp'

filenames = list(askopenfilenames()) # show an "Open" dialog box and return the path to the selected file

timestrings=[]
for fi in filenames:
    timestrings.append(datetime.datetime.fromtimestamp(os.path.getmtime(fi)).isoformat())

# Create a Pandas dataframe from some data.
df = pd.DataFrame({'Filename': filenames,
                   'Cis Conc [M]': np.ones(len(filenames)),
                   'Trans Conc [M]': np.ones(len(filenames)),
                   'Laser Power [mW]': np.zeros(len(filenames)),
                   'Lens [mm]': 0*np.ones(len(filenames)),
                   'IV': np.ones(len(filenames)),
                   'UseForFigure1': np.ones(len(filenames)),
                   'Time': timestrings,
                   'UseForCondEvolution': np.ones(len(filenames)),
                   'Wavelength': 640*np.ones(len(filenames)),
                   'SurfaceChargePlot': np.ones(len(filenames)),
                   'UseForDecay': np.ones(len(filenames)),
                   'pH': 7.4*np.ones(len(filenames))})

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter(str(os.path.split(filenames[0])[0])+os.sep+title+'.xlsx', engine='xlsxwriter')

# Convert the dataframe to an XlsxWriter Excel object.
df.to_excel(writer, sheet_name='Sheet1')

# Close the Pandas Excel writer and output the Excel file.
writer.save()