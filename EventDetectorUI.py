import sys
from PyQt5.QtWidgets import *



class AnalysisUI(QWidget):

    def __init__(self):
        super().__init__()
        self.initUI()

    def AllVariable(self):
        # Let's dump all the variable we will need here...
        self.importedFiles = []

    def initUI(self):
        # Main Window
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



        self.show()

    def MakeFileImportLayout(self):
        self.FileImportLayout = QGroupBox("File Import")
        layout = QGridLayout()
        self.button_fileimport = QPushButton('Open Files')
        layout.addWidget(self.button_fileimport, 0, 0)
        layout.addWidget(QPushButton('2'), 1, 0)
        layout.addWidget(QPushButton('3'), 2, 0)
        layout.addWidget(QPushButton('4'), 3, 0)
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
        QMessageBox.information(self, 'Really?', 'Fuck Off!!')

    def closeEvent(self, event):
        reply = QMessageBox.question(self, 'Message',
                                     "Are you sure to quit? Don't make me sad... :( :( :(", QMessageBox.Yes |
                                     QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

## Execution
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = AnalysisUI()
    sys.exit(app.exec_())
