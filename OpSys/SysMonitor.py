

__all__ = ["getTemprature", "GUISysMonitor"]


N_MAX_CORE = 8


import subprocess

def getTemprature(vbose=False):
    p = subprocess.Popen('sensors', stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    out, err = p.communicate()
    if out==None:
        return None
    
    out_str = out.decode('utf-8')
    if vbose:
        print(out_str)

    Temp_Main = None
    str_start = 'Physical id 0:'
    str_end = chr(176)+'C'
    if str_start in out_str:
        ind_start = out_str.find(str_start)+len(str_start)
        ind_end   = out_str.find(str_end, ind_start)
        temp = out_str[ind_start:ind_end].strip()
        
        Temp_Main = float(temp)
        if vbose:
            print('Temp = ', temp)

    Temp_Core = []
    for i in range(N_MAX_CORE):
        str_start = 'Core {}:'.format(i)
        str_end = chr(176)+'C'
        if str_start in out_str:
            ind_start = out_str.find(str_start)+len(str_start)
            ind_end   = out_str.find(str_end, ind_start)
            temp = out_str[ind_start:ind_end].strip()
            
            Temp_Core.append(float(temp))
            if vbose:
                print('Core 0 => ', temp)
        else:
            break
            
    return [Temp_Main, Temp_Core]





#from multiprocessing import Process, Queue
from threading import Thread
import queue
from tkinter import *
import time
from random import randint

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure


class GUISysMonitor(Thread):
    ## TODO: resolve system crash
    
    def __init__(self, que=None):
        self.queue = que
        self.quit = False
        Thread.__init__(self)
        #self.daemon = True
        self.start()
        return
        
    def SetCommChannel(self, que):
        self.queue = que
        return
                
    def run(self):
        self.root = Tk()

        label = Label(self.root, text="Start Page")#, font=LARGE_FONT)
        label.pack(pady=10,padx=10)
        #self.label = label

        
        fig = Figure(figsize=(5,5), dpi=100)
        ax = fig.add_subplot(111)
        ax.plot([1,2,3,4,5,6,7,8],[5,6,1,3,8,9,3,5])
        #self.fig = fig
        #self.ax = ax
        
        figcanvas = FigureCanvasTkAgg(fig, self.root)
        figcanvas.show()
        figcanvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=True)

        #self.figcanvas = figcanvas
        
        
        self.root.protocol("WM_DELETE_WINDOW", self.TerminateRoot)
        
        self.root.after(1000, self.CommandLoop)
        self.root.mainloop()
        print('thread exiting..!')
        return
        
    def TerminateRoot(self):
        print('stopping thread!')
        self.quit = True
        self.root.quit()
        self.root.destroy()
        #del self.figcanvas
        del self.root
        return

    def CommandLoop(self):
        if not self.quit:
            if self.queue != None:
                if not self.queue.empty():                
                    comm = self.queue.get()
                    self.ProcessCommand(comm)
                    self.queue.task_done()
            #self.figcanvas.draw_idle()
            self.root.after(1000, self.CommandLoop)
        return

    def ProcessCommand(self, comm):
        if comm[0]=="plot":
            self.ax.clear()
            self.ax.plot(*tuple(comm[1]))


 


