from PyQt5.QtWidgets import QMainWindow, QApplication
import pyqtgraph.exporters
from PyQt5 import QtGui
from plotgui import Ui_Plotting
import pyqtgraph as pg
import Functions as f
import sys
import numpy as np
import os
import glob
import LoadData


#   Draw UI
app = QApplication(sys.argv)
window = QMainWindow()
ui = Ui_Plotting()
ui.setupUi(window)

folder = ''
currentFile = ''
ui.graphicsView.setBackground('w')
ui.graphicsView_2.setBackground('w')
ui.graphicsView_3.setBackground('w')


def ChangedFolder():
    global folder
    ui.listWidget.clear()
    folder = QtGui.QFileDialog.getExistingDirectory(window, 'Select in which Folder to save your data', os.getcwd())
    if folder:
        print('Opening ' + folder)
        os.chdir(folder)
        file_list = glob.glob(folder + os.sep + '*.dat')
        file_list.extend(glob.glob(folder + os.sep + '*.log'))
        file_list.extend(glob.glob(folder + os.sep + '*.abf'))
        file_list.sort(key = os.path.getmtime)
        for i in file_list:
            ui.listWidget.addItem(os.path.split(i)[1])

def ItemChanged(it):
    print(it.text())
    global currentFile
    currentFile = it.text()
    # Update Plots
    UpdatePlots(currentFile)
    ui.label.setText(currentFile)

def UpdatePlots(currentFile):
    inp = LoadData.OpenFile(folder + os.sep + currentFile, None, True)
    ui.graphicsView.plotItem.clear()
    ui.graphicsView_2.plotItem.clear()
    ui.graphicsView_3.plotItem.clear()
    if 'i1' in inp:
        ui.graphicsView.plot(y=inp['i1'], x=np.arange(0, len(inp['i1']))/ inp['samplerate'], pen = 'k')
    if 'i2' in inp:
        ui.graphicsView_3.plot(y=inp['i2'], x=np.arange(0, len(inp['i2']))/inp['samplerate'], pen = 'b')
    if 'v2' in inp:
        ui.graphicsView_2.plot(y=inp['v2'], x=np.arange(0, len(inp['i2']))/inp['samplerate'], pen = 'b')
    if 'v1' in inp:
        ui.graphicsView_2.plot(y=inp['v1'], x=np.arange(0, len(inp['i1']))/inp['samplerate'], pen = 'k')

def SaveWindow():
    win = pg.GraphicsWindow()
    #win.addItem(ui.graphicsView.plotItem)
    win.show()
    exporter = pg.exporters.ImageExporter(ui.graphicsView.plotItem)
    exporter.parameters()['width'] = 10000  # (note this also affects height parameter)
    exporter.export('fileName.png')
    exporter2 = pg.exporters.ImageExporter(ui.graphicsView_3.plotItem)
    exporter2.parameters()['width'] = 10000  # (note this also affects height parameter)
    exporter2.export('fileName_1.png')


ui.listWidget.itemDoubleClicked.connect(ItemChanged)
ui.pushButton.pressed.connect(ChangedFolder)
ui.pushSave.pressed.connect(SaveWindow)

## Plotting Settings
ui.graphicsView.plotItem.setDownsampling(ds=10, auto=True, mode='subsample')
ui.graphicsView.plotItem.setClipToView(True)
pg.setConfigOptions(antialias=True)
pg.setConfigOptions(background=None)
ui.graphicsView.plotItem.setLabel('left', "Current", units='A')
ui.graphicsView.plotItem.setLabel('bottom', "Time", units='s')

ui.graphicsView_2.plotItem.setDownsampling(ds=10, auto=True, mode='subsample')
ui.graphicsView_2.plotItem.setClipToView(True)
ui.graphicsView_2.plotItem.setXLink(ui.graphicsView.plotItem)
ui.graphicsView_2.plotItem.setLabel('left', "Voltage", units='V')
ui.graphicsView_2.plotItem.setLabel('bottom', "Time", units='s')

ui.graphicsView_3.plotItem.setDownsampling(ds=100, auto=True, mode='subsample')
ui.graphicsView_3.plotItem.setClipToView(True)
ui.graphicsView_3.plotItem.setXLink(ui.graphicsView.plotItem)
ui.graphicsView_3.plotItem.setLabel('left', "Current", units='A')
ui.graphicsView_3.plotItem.setLabel('bottom', "Time", units='s')

window.show()
sys.exit(app.exec_())