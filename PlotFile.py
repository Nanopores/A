import numpy as np
import scipy
from scipy.optimize import curve_fit
from scipy import constants as cst
import Functions as uf
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import h5py
import pandas as pd
from matplotlib.ticker import EngFormatter
from matplotlib import rc
#rc('text', usetex=True)
formCurrent = EngFormatter(unit='A', places=3)
formA = EngFormatter(unit='A', places=3)
formB = EngFormatter(unit='s', places=3)

file = '/Volumes/lben/lben-commun/2018 User Data/Michael/Axopatch/20180504/W11R149_ph74_1M_1mM_KCl_100mW_relaxation_473.dat'
dat = uf.OpenFile(file)

fig1 = plt.figure(1, figsize=(8, 4))
ax_part1 = fig1.add_subplot(1, 1, 1)

ax_part1.plot(dat['i1'])

plt.show()