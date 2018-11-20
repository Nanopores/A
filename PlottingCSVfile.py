import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from tkinter.filedialog import askopenfilenames


#ONE CAN ADAPT THIS CODE TO MAKE A SIMPLE 2D GRAPH Updated 29-08-2018

filenames=askopenfilenames()
for filename in filenames: #Update your column names
    data_file = pd.read_csv(filename) #For Excel files
    column1 = data_file ['Time (ks)'] #x-axis data
    column2 = data_file ['Current (nA)'] #y-axis data 1
#column3 = data_file ['Put column 3 name here'] #y-axis data 2 #Uncomment if you have multiple y-axis data e.g. column 3 etc

    # Aesthetics: comment the part that starts from marker stuff if not required
    plt.plot (column1, column2, marker='o', markerfacecolor='blue', markersize=5, color='skyblue', linewidth=4)
    #plt.plot (column1, column3, marker='', color='olive', linewidth=2, linestyle='dashed', label="toto")

    plt.xlabel('Time (ks)') and plt.ylabel('Current (nA)')
    plt.show()
#print (data_file) just in case if u want to see wht you are actually printing, sometimes the plot is reversed!!!