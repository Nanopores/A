import pandas as pd
from scipy.optimize import curve_fit
import Functions as uf
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
from scipy import constants as cst
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def Exp(x, a, b):
    return np.exp(x*a)+b

conc = [1.0, 0.1, 0.01, 0.001, 0.0001]
cond = [10.554, 12.9e-1, 1.41e-1, 0.147e-1, 0.015e-1]
p = np.polyfit(np.log(conc), np.log(cond), 3)
print(p)

concrange = np.linspace(1e-6, 1, 1000)

fig = plt.figure(1, figsize=(5, 5))
ax = fig.add_subplot(111)
ax.plot(np.log(conc), np.log(cond), linestyle = 'None', marker='o', color='b')
ax.plot(np.log(concrange), np.polyval(p,np.log(concrange)), linestyle = '--', color='r')
#ax.set_xscale('log')
#ax.set_yscale('log')

plt.show()