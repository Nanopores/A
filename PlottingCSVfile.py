import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from tkinter.filedialog import askopenfilenames


#Updated 23-11-2018

filenames=askopenfilenames()
for filename in filenames: #Update your column names
    data_file = pd.read_csv(filename) #For Excel files or csv files change yourself
    column1 = data_file ['x000'] #x-axis data
    column2 = data_file ['y000)'] #y-axis data 1
    #column3 = data_file ['Put column 3 name here'] #y-axis data 2 #Uncomment if you have multiple y-axis data e.g. column 3 etc

    # Aesthetics: comment the part that starts from marker stuff if not required
    plt.plot (column1, column2, marker='.', markerfacecolor='blue', markersize=12, color='red', linewidth=4)
    plt.xscale('log')
    plt.title('MoS2 nanopore ionic current, 300 mV, 1M KCl, pH 7.5')
    plt.grid(True)
    #plt.plot (column1, column3, marker='', color='olive', linewidth=2, linestyle='dashed', label="toto")
    plt.xlabel('Time (ks)') and plt.ylabel('Current (nA)')
    plt.show()
#print (data_file) just in case if u want to see wht you are actually printing, sometimes the plot is reversed!!!

