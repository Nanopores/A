import sys
from PyQt5.QtWidgets import *

import os
import numpy as np
import shelve

import matplotlib
matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from Plotting import EventPlots
import NanoporeClasses as NC

class AnalysisUI(QWidget):

    def __init__(self):
        super().__init__()
        # self.AllVariable()


    #def AllVariable(self):
        # Let's dump all the variable we will need here...
        self.importedFiles = []
        self.selectedFiles = []
        self.analyzedFiles = []

        # Defining the plots
        self.TurnToolbarsOn = False
        self.pltFig = plt.figure(1, figsize=(12, 10))
        self.pltFig.__attached = True
        self.figurePlotting = FigureCanvas(figure=self.pltFig)

        self.initUI()

    def initUI(self):
        # Main Window
        self.setGeometry(300, 200, 1200, 850)
        self.setWindowTitle('Translocation event analysis')

        # Assemble the different layout pieces
        self.MakeFileImportLayout()
        self.DisplayAnalysis()
        windowLayout = QHBoxLayout()
        windowLayout.addWidget(self.FileImportLayout)
        windowLayout.addWidget(self.AnalysisResults)
        self.setLayout(windowLayout)

        # All Actions go here:
        self.button_fileimport.clicked.connect(self.ImportButtonPushed)
        self.button_clearfilelist.clicked.connect(self.ClearListButtonPushed)
        self.list_filelist.clicked.connect(self.SelectionInFileListChanged)
        # self.button_startAnalysis.clicked.connect(self.AnalysisButtonPressed)

        self.plotCurrent.clicked.connect(self.PlotSelected)
        self.plotVoltage.clicked.connect(self.PlotSelected)
        self.showLog.clicked.connect(self.PlotSelected)
        self.show()

    def MakeFileImportLayout(self):
        self.FileImportLayout = QGroupBox("File Import")
        self.FileImportLayout.setFixedWidth(400)
        layout = QGridLayout()
        self.button_fileimport = QPushButton('Open Files')
        self.button_fileimport.setToolTip('Select the files you want to import into the list. You can do multiple selections.')
        self.button_clearfilelist = QPushButton('Clear List')
        self.button_clearfilelist.setToolTip('This deletes your current list of files.')
        self.list_filelist = QListWidget()
        self.list_filelist.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.list_filelist.setToolTip('This is the list of current files. selected files will be plotted.')
        layout.addWidget(self.button_fileimport, 0, 0)
        layout.addWidget(self.list_filelist, 1, 0)
        layout.addWidget(self.button_clearfilelist, 2, 0)
        self.FileImportLayout.setLayout(layout)

    def DisplayAnalysis(self):
        self.AnalysisResults = QGroupBox("Analysis results")
        self.Plotsettings = QGroupBox("Plot Settings")
        self.Plotsettings.setSizePolicy(QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum))
        self.Plotting = QGroupBox("Plotting")

        # layout.setColumnStretch(1, 4) # stretch factors
        # layout.setColumnStretch(2, 4)

        self.plotCurrent = QCheckBox('Show Current')
        self.plotVoltage = QCheckBox('Seperate by voltage')
        self.showLog = QCheckBox('Show log')

        layoutsettings = QGridLayout()
        layoutsettings.addWidget(self.plotCurrent, 0, 0)
        layoutsettings.addWidget(self.plotVoltage, 0, 1)
        layoutsettings.addWidget(self.showLog, 0, 2)
        self.Plotsettings.setLayout(layoutsettings)

        layoutPlotting = QGridLayout()
        layoutPlotting.addWidget(self.figurePlotting)
        self.Plotting.setLayout(layoutPlotting)

        # Final Assembly
        layout = QGridLayout()
        layout.addWidget(self.Plotsettings, 0, 0)
        layout.addWidget(self.Plotting, 1, 0)
        self.AnalysisResults.setLayout(layout)

    # Button Actions
    def ImportButtonPushed(self, event):
        files = QFileDialog.getOpenFileNames(self, 'Select Files to Add to list', filter="data files (*.dat)")
        # Add to file list if unique:
        for file in files[0]:
            if file not in self.importedFiles:
                self.importedFiles.append(file)
        self.UpdateFileList()

    def ClearListButtonPushed(self, event):
        self.importedFiles = []
        self.UpdateFileList()

    def UpdateFileList(self):
        self.list_filelist.clear()
        i = 0
        for file in self.importedFiles:
            item = QListWidgetItem(os.path.split(file)[1], type=i)
            self.list_filelist.addItem(item)
            i += 1

    def SelectionInFileListChanged(self, event):
        items = self.list_filelist.selectedItems()
        selectedFiles = []
        for item in items:
            idx = item.type()
            selectedFiles.append(str(self.importedFiles[idx]))
        self.selectedFiles = selectedFiles
        print('The following files are selected in the list:')
        print(selectedFiles)
        self.PlotSelected()

    def ClearFigure(self):
        self.pltFig.clf()

    def PlotSelected(self):
        self.ClearFigure()
        plotCurrent = self.plotCurrent.isChecked()
        plotVoltage = self.plotVoltage.isChecked()
        showLog = self.showLog.isChecked()

        if len(self.selectedFiles)>0:
            translocationEvents = NC.AllEvents()
            for filename in self.selectedFiles:
                shelfFile = shelve.open(os.path.splitext(filename)[0])
                translocationEventstemp = shelfFile['TranslocationEvents']
                shelfFile.close()
                translocationEvents.AddEvent(translocationEventstemp)

            EventPlots.PlotG_tau(translocationEvents, fig=self.pltFig, showCurrent=plotCurrent, showLog=showLog)

            self.figurePlotting.draw()


    def closeEvent(self, event):
        reply = QMessageBox.question(self, 'Message',
                                     "Are you sure to quit?", QMessageBox.Yes |
                                     QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()


# Execution
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = AnalysisUI()
    sys.exit(app.exec_())
