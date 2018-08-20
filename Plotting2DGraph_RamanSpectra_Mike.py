import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


#ONE CAN ADAPT THIS CODE TO MAKE A SIMPLE 2D GRAPH 19-08-2018

data_file = pd.read_excel('raman_spectra_1.xlsx') #Makesure your readble file is in the current directory like in my case M Drive, Michael/Jochem can we make this openfile like thingy?
column1 = data_file ['Raman shift (cm-1)'] #x-axis
column2 = data_file ['MoS2 monolayer traingles'] #y-axis 1
column3 = data_file ['MoS2 monolayer traingles 2 test'] #y-axis 2

plt.plot (column1, column2, marker='o', markerfacecolor='blue', markersize=12, color='skyblue', linewidth=4) #comment the part from marker stuff if not required
plt.plot (column1, column3, marker='', color='olive', linewidth=2, linestyle='dashed', label="toto") #comment the part from marker stuff if not required

plt.xlabel('Raman Shift (cm-1)') and plt.ylabel('Counts (a.u.)')
plt.show()
print (data_file) #just in case if u want to see wht you are actually printing, sometimes the plot is reversed!!!