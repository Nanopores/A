import sys
import os
from PyQt5.QtWidgets import *
from PyQt5 import QtCore
from PyQt5.QtGui import QColor
import matplotlib
matplotlib.use('QT5Agg')
import numpy as np
import matplotlib.pylab as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
font = {'family' : 'monospace',
        'weight' : 'regular',
        'size'   : 4}
matplotlib.rc('font', **font)  # pass in the font dict as kwargs
from EventDetection import *


class AnalysisUI(QWidget):

    def __init__(self):
        super().__init__()
        self.initUI()

    def AllVariable(self):
        # Let's dump all the variable we will need here...
        self.importedFiles = []
        self.selectedFiles = []
        self.analyzedFiles = []

        # Defining the plots
        self.TurnToolbarsOn = False
        self.fig_singleEvent = plt.figure(1, figsize=(10, 10))
        self.ax_singleEvent = self.fig_singleEvent.add_subplot(111)
        self.ax_singleEvent.plot(np.arange(1, 10), np.arange(1, 10))

        self.fig_wholeTrace = plt.figure(2, figsize=(16, 9))
        self.ax_wholeTrace = self.fig_wholeTrace.add_subplot(111)
        self.ax_wholeTrace.plot(np.arange(1, 10), -np.arange(1, 10))

        self.figure_singleEvent = FigureCanvas(figure=self.fig_singleEvent)
        self.figure_wholeTrace = FigureCanvas(figure=self.fig_wholeTrace)
        if self.TurnToolbarsOn:
            self.toolbarWholeTrace = NavigationToolbar(self.figure_wholeTrace, self)
            self.toolbarSingleEvent = NavigationToolbar(self.figure_singleEvent, self)

    def initUI(self):
        # Main Window
        self.AllVariable()
        self.setGeometry(0, 0, 1200, 750)
        self.setWindowTitle('Do not use this tool for shitty science...')

        # Assemble the different layout pieces
        self.MakeFileImportLayout()
        self.MakeAnalysisParametersAndRunLayout()
        self.MakeEventNavigatorLayout()
        windowLayout = QHBoxLayout()
        windowLayout.addWidget(self.FileImportLayout, 10)
        windowLayout.addWidget(self.AnalysisParameters)
        windowLayout.addWidget(self.EventNavigatorLayout)
        self.setLayout(windowLayout)

        ## All Actions go here:
        self.button_fileimport.clicked.connect(self.ImportButtonPushed)
        self.button_clearfilelist.clicked.connect(self.ClearListButtonPushed)
        self.list_filelist.clicked.connect(self.SelectionInFileListChanged)
        self.button_startAnalysis.clicked.connect(self.AnalysisButtonPressed)
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
        self.list_filelist.setToolTip('This is the list of current files. Entries in green are already analyzed. Your selection defines which files are used in the analysis.')
        layout.addWidget(self.button_fileimport, 0, 0)
        layout.addWidget(self.list_filelist, 1, 0)
        layout.addWidget(self.button_clearfilelist, 2, 0)
        self.FileImportLayout.setLayout(layout)

    def MakeAnalysisParametersAndRunLayout(self):
        self.AnalysisParameters = QGroupBox("Analysis")
        self.LowPassSettings = QGroupBox("Low Pass Settings")
        self.CusumSettings = QGroupBox("Cusum Settings")
        self.CusumSettings.setSizePolicy(QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum))     #This allows this part to remain small in vertical position when the UI resizes.
        self.OtherSettings = QGroupBox("Other Settings")
        self.button_startAnalysis = QPushButton('Start Analysis')
        self.button_startAnalysis.setToolTip('This button will provoke your computer to sel-destruct. Use at your own risk!!!!!')

        ## Low Pass Settings
        self.settings_LP_a = QComboBox()
        availableoptions = [0.99, 0.999, 0.9999]
        for i in availableoptions:
            self.settings_LP_a.addItem(str(i))
        self.settings_LP_a.setToolTip('How strong should the filter be? Higher number is stronger filtering, which means less found events...')


        ## EndValueE
        self.settings_LP_E = QDoubleSpinBox()
        self.settings_LP_E.setValue(0)
        self.settings_LP_E.setSingleStep(0.1)
        self.settings_LP_E.setToolTip('The current needs to go above mean + E*std to terminate the event')

        ## StartValueS
        self.settings_LP_S = QDoubleSpinBox()
        self.settings_LP_S.setValue(5)
        self.settings_LP_S.setSingleStep(0.1)
        self.settings_LP_S.setToolTip('The current needs to go below mean + S*std to initiate an event')

        layout1 = QGridLayout()
        layout1.addWidget(QLabel('Low Pass Filter Coefficient (a)'), 0, 0)
        layout1.addWidget(self.settings_LP_a, 0, 1)
        layout1.addWidget(QLabel('Start Coefficient (S)'), 1, 0)
        layout1.addWidget(self.settings_LP_S, 1, 1)
        layout1.addWidget(QLabel('End Coefficient (E)'), 3, 0)
        layout1.addWidget(self.settings_LP_E, 3, 1)
        self.LowPassSettings.setLayout(layout1)


        ##Cusum Settings
        self.settings_cusum_hBook = QDoubleSpinBox()
        self.settings_cusum_hBook.setValue(1)
        self.settings_cusum_hBook.setSingleStep(0.1)
        self.settings_cusum_hBook.setToolTip('Some random stupid parameter')


        self.settings_cusum_dI = QDoubleSpinBox()
        self.settings_cusum_dI.setValue(1)
        self.settings_cusum_dI.setSuffix(' nA')
        self.settings_cusum_dI.setSingleStep(0.1)
        self.settings_cusum_dI.setToolTip('The most probable current drop')


        layoutCusum = QGridLayout()
        layoutCusum.addWidget(QLabel('h-Book'), 0, 0)
        layoutCusum.addWidget(self.settings_cusum_hBook, 0, 1)
        layoutCusum.addWidget(QLabel('Delta I'), 1, 0)
        layoutCusum.addWidget(self.settings_cusum_dI, 1, 1)
        self.CusumSettings.setLayout(layoutCusum)


        ## Other Settings
        self.settings_other_ChimeraLP = QDoubleSpinBox()
        self.settings_other_ChimeraLP.setValue(10)
        self.settings_other_ChimeraLP.setSingleStep(10)
        self.settings_other_ChimeraLP.setSuffix(' kHz')
        self.settings_other_ChimeraLP.setToolTip('Low pass setting for high-bandwidth Chimera Data')


        self.settings_other_MaxEventLength = QDoubleSpinBox()
        self.settings_other_MaxEventLength.setValue(0.2)
        self.settings_other_MaxEventLength.setSingleStep(0.05)
        self.settings_other_MaxEventLength.setSuffix(' s')
        self.settings_other_MaxEventLength.setToolTip('maximal event length to be considered an event')


        self.settings_other_MinEventLength = QDoubleSpinBox()
        self.settings_other_MinEventLength.setValue(600)
        self.settings_other_MinEventLength.setSingleStep(100)
        self.settings_other_MinEventLength.setSuffix(' us')
        self.settings_other_MinEventLength.setToolTip('maximal event length to be considered an impulse')


        self.settings_other_fitLength = QDoubleSpinBox()
        self.settings_other_fitLength.setValue(3)
        self.settings_other_fitLength.setSingleStep(0.5)
        self.settings_other_fitLength.setSuffix(' ms')
        self.settings_other_fitLength.setToolTip('minimal event length to be fitted for impulses')

        layoutOther = QGridLayout()
        layoutOther.addWidget(QLabel('Chimera Low Pass'), 0, 0)
        layoutOther.addWidget(self.settings_other_ChimeraLP, 0, 1)
        layoutOther.addWidget(QLabel('Max event length'), 1, 0)
        layoutOther.addWidget(self.settings_other_MaxEventLength, 1, 1)
        layoutOther.addWidget(QLabel('Min event length'), 2, 0)
        layoutOther.addWidget(self.settings_other_MinEventLength, 2, 1)
        layoutOther.addWidget(QLabel('Fit length'), 3, 0)
        layoutOther.addWidget(self.settings_other_fitLength, 3, 1)
        self.OtherSettings.setLayout(layoutOther)

        ##Final Assembly
        layout = QGridLayout()
        layout.addWidget(self.LowPassSettings, 0, 0)
        layout.addWidget(self.CusumSettings, 1, 0)
        layout.addWidget(self.OtherSettings, 2, 0)
        layout.addWidget(self.button_startAnalysis, 3, 0)
        self.AnalysisParameters.setLayout(layout)

    def MakeEventNavigatorLayout(self):
        self.EventNavigatorLayout = QGroupBox("Event Navigator")
        self.navigation_buttons =QButtonGroup()
        self.button_eventForward = QPushButton('Forward')
        self.button_eventForward.setToolTip('Advance one event forward')
        self.button_eventBackward = QPushButton('Backward')
        self.button_eventBackward.setToolTip('Go one event back')
        self.text_eventNumber = QLineEdit('Current Event Number')
        self.text_eventNumber.setToolTip('This is the number of the current event. You can edit this field to jump directly to an event.')
        self.text_eventNumber.setAlignment(QtCore.Qt.AlignCenter)
        layout = QGridLayout()
        layout.addWidget(self.figure_singleEvent, 0, 0, 1, 3)
        layout.addWidget(self.button_eventForward, 2, 2)
        layout.addWidget(self.text_eventNumber, 2, 1)
        layout.addWidget(self.button_eventBackward, 2, 0)
        layout.addWidget(self.figure_wholeTrace, 3, 0, 1, 3)
        if self.TurnToolbarsOn:
            layout.addWidget(self.toolbarSingleEvent, 1, 0, 1, 3)
            layout.addWidget(self.toolbarWholeTrace, 4, 0, 1, 3)
        self.EventNavigatorLayout.setLayout(layout)


## Button Actions

    def ImportButtonPushed(self, event):
        files = QFileDialog.getOpenFileNames(self, 'Select Files to Add to list', filter = "Images (*.dat *.log)")
        ## Add to file list if unique:
        for i in files[0]:
            if i not in self.importedFiles:
                self.importedFiles.append(i)
        self.UpdateFileList()

    def ClearListButtonPushed(self, event):
        self.importedFiles = []
        self.UpdateFileList()

    def UpdateFileList(self):
        self.list_filelist.clear()
        for ind, i in enumerate(self.importedFiles):
            self.list_filelist.addItem(os.path.split(i)[1])
            folder = os.path.dirname(i) + os.sep + 'analysisfiles'
            filename, file_extension = os.path.splitext(os.path.basename(i))
            potentialanalysisfile = folder + os.sep + filename + 'data.db'
            print(potentialanalysisfile)
            print(os.path.exists(potentialanalysisfile))
            if os.path.isfile(potentialanalysisfile):
            ## If file is present, make the row green. This means analysis was done on it.
                self.list_filelist.item(ind).setBackground(QColor('green'))

    def SelectionInFileListChanged(self, event):
        items = self.list_filelist.selectedItems()
        x = []
        for i in range(len(items)):
            x.append(str(self.list_filelist.selectedItems()[i].text()))
        self.selectedFiles = x
        print('The following files are selected in the list:')
        print(x)

    def AnalysisButtonPressed(self):
        # Get Full File Paths:
        fullfilepaths = []
        for i in self.selectedFiles:
            for j in self.importedFiles:
                if i in j:
                    fullfilepaths.append(j)
        print(fullfilepaths)
        # Fill coefficient dictionary from the user inputs
        coefficients = {}
        coefficients['maxEventLength'] = self.settings_other_MaxEventLength.value()
        coefficients['minEventLength'] = self.settings_other_MinEventLength.value()
        coefficients['fitLength'] = self.settings_other_fitLength.value()
        coefficients['delta'] = self.settings_cusum_dI.value()
        coefficients['hbook'] = self.settings_cusum_hBook.value()
        coefficients['a'] = np.float(self.settings_LP_a.currentText())
        coefficients['S'] = self.settings_LP_S.value()
        coefficients['E'] = self.settings_LP_E.value()
        coefficients['ChimeraLowPass'] = self.settings_other_ChimeraLP.value()
        coefficients['dt'] = 25
        coefficients['deltaRel'] = None
        for i in fullfilepaths:
            eventdetectionwithfilegeneration(i, coefficients, forceRun=True, verboseLevel=1)
        self.UpdateFileList()
    def closeEvent(self, event):
        '''
        reply = QMessageBox.question(self, 'Message',
                                     "Are you sure to quit? Don't make me sad... :( :( :(", QMessageBox.Yes |
                                     QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()
        '''

## Execution
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = AnalysisUI()
    sys.exit(app.exec_())
