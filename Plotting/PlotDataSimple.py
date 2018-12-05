from pyqtgraph.Qt import QtGui, QtCore
from pyqtgraph import FileDialog
from pyqtgraph import PlotDataItem
from pyqtgraph import PlotItem
import types
from types import MethodType
import scipy.signal as sig
import numpy as np
import scipy
import pyqtgraph as pg
import Functions as f
import LoadData


inp = LoadData.OpenFile('/Volumes/Michael/Axopatch/20180316/1msSignal_AppliedToModelCellThroughAxopatch_1.dat')

#Downsampling Methods
ds_factor = 10.0

#QtGui.QApplication.setGraphicsSystem('raster')
app = QtGui.QApplication([])
#mw = QtGui.QMainWindow()
#mw.resize(800,800)

win = pg.GraphicsWindow(title="Basic plotting examples")
win.resize(1000, 600)

win.setWindowTitle('pyqtgraph example: Plotting')
"""
p0 = win.addPlot(y=inp['i1'], x=np.arange(0, len(inp['i1']))/inp['samplerate'])
lr = pg.LinearRegionItem([1000, len(inp['i1']-1000)])
lr.setZValue(-10)
p0.addItem(lr)
p0.hideButtons()
#p0.enableAutoRange(y=False, x=False)
p0.setMenuEnabled(enableMenu=False, enableViewBoxMenu='same')
pg.setConfigOptions(antialias=True)
p0.setLabel('left', "Current", units='A')
p0.setLabel('bottom', "Time", units='s')
p0.setMouseEnabled(x=False, y=False)
win.nextRow()
"""
p1 = win.addPlot(title='Ionic Current')
p1.setDownsampling(ds=10, auto=True, mode='subsample')
p1.setClipToView(True)
pg.setConfigOptions(antialias=True)
p1.plot(y=inp['i1'], x=np.arange(0, len(inp['i1']))/inp['samplerate'])
p1.setLabel('left', "Current", units='A')
p1.setLabel('bottom', "Time", units='s')

win.nextRow()

p2 = win.addPlot(title='Transverse Current')
p2.setDownsampling(ds=10, auto=True, mode='subsample')
p2.setClipToView(True)
pg.setConfigOptions(antialias=True)
p2.plot(y=inp['i2'], x=np.arange(0, len(inp['i2']))/inp['samplerate'], pen='r')
p2.setXLink(p1)
p2.setLabel('left', "Current", units='A')
p2.setLabel('bottom', "Time", units='s')

win.nextRow()

p3 = win.addPlot(title='Voltages')
p3.setDownsampling(ds=100, auto=True, mode='subsample')
p3.setClipToView(True)
pg.setConfigOptions(antialias=True)
p3.plot(y=inp['v2'], x=np.arange(0, len(inp['i2']))/inp['samplerate'], pen='r')
p3.plot(y=inp['v1'], x=np.arange(0, len(inp['i2']))/inp['samplerate'])
p3.setXLink(p1)
p3.setLabel('left', "Current", units='A')
p3.setLabel('bottom', "Time", units='s')
"""
def updatePlot():
    p1.setXRange(*lr.getRegion(), padding=0)
def updateRegion():
    lr.setRegion(p1.getViewBox().viewRange()[0])
lr.sigRegionChanged.connect(updatePlot)
p1.sigXRangeChanged.connect(updateRegion)
updatePlot()
"""


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
