import sys
import os
from PyQt5.QtWidgets import *



class AnalysisUI(QWidget):

    def __init__(self):
        super().__init__()
        self.initUI()

    def AllVariable(self):
        # Let's dump all the variable we will need here...
        self.importedFiles = []
        self.selectedFiles = []

    def initUI(self):
        # Main Window
        self.AllVariable()
        self.setGeometry(300, 300, 500, 500)
        self.setWindowTitle('This will be our interface for the analysis')

        # Assemble the different layout pieces
        self.MakeFileImportLayout()
        self.MakeAnalysisParametersAndRunLayout()
        windowLayout = QHBoxLayout()
        windowLayout.addWidget(self.FileImportLayout)
        windowLayout.addWidget(self.AnalysisParameters)
        self.setLayout(windowLayout)

        ## All Actions go here:
        self.button_fileimport.clicked.connect(self.ImportButtonPushed)
        self.button_clearfilelist.clicked.connect(self.ClearListButtonPushed)
        self.list_filelist.clicked.connect(self.SelectionInFileListChanged)
        self.show()

    def MakeFileImportLayout(self):
        self.FileImportLayout = QGroupBox("File Import")
        layout = QGridLayout()
        self.button_fileimport = QPushButton('Open Files')
        self.button_clearfilelist = QPushButton('Clear List')
        self.list_filelist = QListWidget()
        self.list_filelist.setSelectionMode(QAbstractItemView.ExtendedSelection)
        layout.addWidget(self.button_fileimport, 0, 0)
        layout.addWidget(self.list_filelist, 1, 0)
        layout.addWidget(self.button_clearfilelist, 2, 0)
        self.FileImportLayout.setLayout(layout)


    def MakeAnalysisParametersAndRunLayout(self):
        self.AnalysisParameters = QGroupBox("Analysis Parameters")
        layout = QGridLayout()
        layout.setColumnStretch(1, 4) # stretch factors
        layout.setColumnStretch(2, 4)
        layout.addWidget(QPushButton('1'), 0, 0)
        layout.addWidget(QPushButton('2'), 0, 1)
        layout.addWidget(QPushButton('3'), 0, 2)
        layout.addWidget(QPushButton('4'), 1, 0)
        layout.addWidget(QPushButton('5'), 1, 1)
        layout.addWidget(QPushButton('6'), 1, 2)
        layout.addWidget(QPushButton('7'), 2, 0)
        layout.addWidget(QPushButton('8'), 2, 1)
        layout.addWidget(QPushButton('9'), 2, 2)
        self.AnalysisParameters.setLayout(layout)

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

    def closeEvent(self, event):
        reply = QMessageBox.question(self, 'Message',
                                     "Are you sure to quit? Don't make me sad... :( :( :(", QMessageBox.Yes |
                                     QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def UpdateFileList(self):
        self.list_filelist.clear()
        for i in self.importedFiles:
            self.list_filelist.addItem(os.path.split(i)[1])

    def SelectionInFileListChanged(self, event):
        items = self.list_filelist.selectedItems()
        x = []
        for i in range(len(items)):
            x.append(str(self.list_filelist.selectedItems()[i].text()))
        self.selectedFiles = x
        print('The following files are selected in the list:')
        print(x)


## Execution
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = AnalysisUI()
    sys.exit(app.exec_())
