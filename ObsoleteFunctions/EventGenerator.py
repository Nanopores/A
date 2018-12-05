##SHOULD NOT BE USED FOR MINI-MEETINGS!!! PAPERS ONLY

import numpy as np
from matplotlib.ticker import FuncFormatter
from scipy import constants as cst
import Functions as uf
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

currForm = EngFormatter(unit='A', places=1)
VoltForm = EngFormatter(unit='V', places=1)

sr = 10e3
baselinetime = 50e-3
baselineCurrent = 5e-9
dropPercent = 0.0
dt = 25e-3
std = 100e-12

signal = np.array([])
signal = np.append(signal, baselineCurrent*np.ones(int(baselinetime*sr)))
signal = np.append(signal, (baselineCurrent-dropPercent*baselineCurrent)*np.ones(int(dt*sr)))
signal = np.append(signal, baselineCurrent*np.ones(int(baselinetime*sr)))
noise = np.random.normal(0, std, len(signal))
signal = signal + noise
t = np.arange(len(signal))/sr

fig2 = plt.figure(2, figsize=(8, 4))
ax2 = fig2.add_subplot(1, 1, 1)
ax2.plot(t, signal)
ax2.yaxis.set_major_formatter(EngFormatter(unit='A'))
ax2.xaxis.set_major_formatter(EngFormatter(unit='s'))
ax2.set_ylabel('Ionic Current')
ax2.set_xlabel('Time')
fig2.savefig('FakeEvent_.pdf', transparent=True)
plt.close()