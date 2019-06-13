import sys
from PyQt5.QtWidgets import *


class AnalysisUI(QWidget):

    def __init__(self):
        super().__init__()
        self.AllVariable()
        self.initUI()

    def AllVariable(self):
        # Let's dump all the variable we will need here...
        self.importedFiles = []
        self.selectedFiles = []
        self.analyzedFiles = []

    def initUI(self):
        # Main Window
        self.setGeometry(300, 200, 1200, 750)
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

    def DisplayAnalysis(self):
        self.AnalysisResults = QGroupBox("Analysis results")
        self.Plotsettings = QGroupBox("Plot Settings")
        self.Plotsettings.setSizePolicy(QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum))
        self.Plotting = QGroupBox("Plotting")

        # layout.setColumnStretch(1, 4) # stretch factors
        # layout.setColumnStretch(2, 4)

        self.plotCurrent = QCheckBox('Show Current')
        self.plotVoltage = QCheckBox('Seperate events by voltage')

        self.showImpulse = QCheckBox('Show impulse events')
        self.showImpulse.setChecked(True)
        self.showCUSUM = QCheckBox('Show CUSUM fitted events')
        self.showCUSUM.setChecked(True)
        self.showRough = QCheckBox('Show non-fitted events')

        layoutsettings = QGridLayout()
        layoutsettings.addWidget(self.plotCurrent, 0, 0)
        layoutsettings.addWidget(self.plotVoltage, 0, 1)
        layoutsettings.addWidget(self.showImpulse, 0, 2)
        layoutsettings.addWidget(self.showCUSUM, 0, 3)
        layoutsettings.addWidget(self.showRough, 0, 4)
        self.Plotsettings.setLayout(layoutsettings)

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
        for i in self.importedFiles:
            self.list_filelist.addItem(os.path.split(i)[1])
        ## If file is present, make the row green. This means analysis was done on it.

    def SelectionInFileListChanged(self, event):
        items = self.list_filelist.selectedItems()
        x = []
        for i in range(len(items)):
            x.append(str(self.list_filelist.selectedItems()[i].text()))
        self.selectedFiles = x
        print('The following files are selected in the list:')
        print(x)


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
