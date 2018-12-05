import matplotlib
matplotlib.use('module://kivy.garden.matplotlib.backend_kivy')
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import Functions as f
import os
import kivy
from kivy.app import App
from kivy.lang import Builder
from kivy.uix.boxlayout import BoxLayout
from kivy.garden.matplotlib.backend_kivyagg import FigureCanvas
import kivy.uix.filechooser as fc
from kivy.uix.button import Button


#inp = f.OpenFile(filename, 10e3)
#file = os.sep + str(os.path.split(filename)[1][:-4])
fig = plt.figure(1, figsize=(10, 4))
ax = fig.add_subplot(111)
canvas = FigureCanvas(fig)

def FileSelected(a,b,c):
    fileAdress = b[0]
    plt.clf()
    print(fileAdress)
    inp = f.OpenFile(fileAdress, 10e3)
    ax.plot(np.arange(0, len(inp['i1']))/inp['samplerate'], inp['i1'])
    canvas.draw_idle()

class MyApp(App):

    def build(self):
        box = BoxLayout(orientation = 'vertical')
        brows = fc.FileChooserListView(path = '/Volumes/lben/lben-commun/2018 User Data/Michael/')
        brows.bind(on_submit=FileSelected)
        box.add_widget(brows)
        box.add_widget(canvas)
        btnlayout = BoxLayout(orientation = 'horizontal')
        btn = Button(text='Ok')
        btn.bind()
        btnlayout.add_widget(btn)
        btn = Button(text='Cancel')
        btn.bind()
        btnlayout.add_widget(btn)
        box.add_widget(btnlayout)
        return box

MyApp().run()