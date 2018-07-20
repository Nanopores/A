import AnalysisParameters as pm
import time
import MiscParameters as pm
import h5py
from openpyxl import Workbook
from tkinter.filedialog import askopenfilenames
pm.init(LoadFiles = 0)

filenames = askopenfilenames(filetypes=[('HDF5 files', '*.hdf5')])
ExperimentName = pm.OutputFolder + '02B_10mM_1M_1kb_pH10.xlsx'
wb = Workbook()
ws1 = wb.active

for (i,k) in enumerate(filenames):
    ws1['A{}'.format(i+1)] = k
    f = h5py.File(k, 'r')
    ws1['B{}'.format(i + 1)] = f['General/TransverseRecorded'].value
    ws1['C{}'.format(i + 1)] = time.ctime(f['General/TimeFileWritten'].value)

wb.save(filename = ExperimentName)