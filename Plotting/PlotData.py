from pyqtgraph.Qt import QtGui, QtCore
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
# Downsample data
downsampled_data = sig.resample(np.float64(inp['i1']), np.uint64(len(inp['i1'])/ds_factor))
downsampled_data_x = np.arange(0, len(downsampled_data))/(inp['samplerate']/ds_factor)

# Replacement for the default getData function
def getData2(obj):
    print('In Function')
    # Calculate the visible range
    range = obj.viewRect()
    print(range)
    if range is not None:
        dx = float(inp['i1'][1] - inp['i1'][0]) / (inp['i1'].size[0] - 1)
        x0 = (range.left() - inp['i1'][0]) / dx
        x1 = (range.right() - inp['i1'][0]) / dx
    # Decide whether to use downsampled or original data
    if (x1 - x0) > 20000:
        print('Is Downsampled')
        obj.xData = downsampled_data_x
        obj.yData = downsampled_data
    else:
        print('Not Downsampled')
        obj.xData = np.arange(0, len(inp['i1']))/inp['samplerate']
        obj.yData = inp['i1']
    # Run the original getData of PlotDataItem
    return PlotDataItem.getData(obj)


#QtGui.QApplication.setGraphicsSystem('raster')
app = QtGui.QApplication([])
#mw = QtGui.QMainWindow()
#mw.resize(800,800)

win = pg.GraphicsWindow(title="Basic plotting examples")
win.resize(1000, 600)
win.setWindowTitle('pyqtgraph example: Plotting')

p1 = PlotItem(y=downsampled_data, x=np.arange(0, len(downsampled_data))/(inp['samplerate']/ds_factor))
p1.getData = types.MethodType(getData2, p1)

# Enable antialiasing for prettier plots
pg.setConfigOptions(antialias=True)
#p1.addItem(win)
win.addItem(p1)
#p1 = win.addPlot(title="Ionic Current Trace", y=inp['i1'], x=np.arange(0,len(inp['i1']))/inp['samplerate'])
p1.setLabel('left', "Current", units='A')
p1.setLabel('bottom', "Time", units='s')

## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
