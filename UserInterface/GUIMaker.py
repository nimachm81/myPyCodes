

__init__ = ["GUIMakerND"]
                
                    
##-----------------------------------------------------------------------
##---------------------- GUI  ------------

import numpy as np
import scipy as sp
import random
import math

import sympy
from sympy import Symbol, symbols, lambdify
from sympy.parsing.sympy_parser import parse_expr

from enum import Enum, IntEnum

from threading import Thread
import queue
from tkinter import *
import time
from random import randint


import wx
import wx.grid
import sys
from wx import glcanvas
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *


class NameGenerator():
    def __init__(self):
        self.names = []
        
    def NameExists(self, name):
        for n in self.names:
            if n==name:
                return True
        return False
        
    def GenerateNewName(self, name):
        i = 1
        name_next = name
        while self.NameExists(name_next):
            name_next = name + str(i)
            i += 1
        return name_next

    def AddName(self, name):
        self.names.append(name)


nameGeneraor = NameGenerator()


class GUIMakerND(Thread):
    
    def __init__(self, W_3D=np.array([1.0, 1.0, 1.0]), que=None):
        self.queue = que
        self.app = None
        self.W_3D = W_3D
        
        Thread.__init__(self)
        #self.daemon = True
        self.start()
        return
        
    def SetCommChannel(self, que):
        self.queue = que
        if self.app!=None:
            self.app.SetQueue(self.queue)
        return

    def run(self):
        self.app = InteractiveGLApp(W_3D=self.W_3D)
        self.app.SetQueue(self.queue)
        #print('GUI queue: ', self.queue)

        self.app.SetCommandLoop()
        self.app.MainLoop()
        self.app.Destroy()

        self.app.deleteGUIItems()
        del self.app




class InteractiveGLApp(wx.App):
    
    def __init__(self, W_3D=np.array([1.0, 1.0, 1.0]), que=None):
        wx.App.__init__(self, redirect=False)
        
        sys.excepthook = MyExceptionHook
        
        self.W_3D = W_3D        # world dimensions
        self.dW_3D = W_3D/10.0
        self.V_3D = np.array([0.0, 0.0, 0.0])        # View origin
        self.V_phi = 20.0       # view angle phi
        self.V_theta = 30.0     # view angle theta
        
        self.SetFrameDimensions()
        self.queue = que
        self.frame.queue = self.queue
        
        self.modif_state = "select"
        
        
    def SetQueue(self, que):
        self.queue = que
        self.frame.queue = self.queue
        #print('app queue: ', self.queue)

    def OnInit(self, W_3D=np.array([1.0, 1.0, 1.0])):
        frame = GLFrame(None, -1, "ND Space", pos=(0,0), size=(800, 600),
                        style=wx.DEFAULT_FRAME_STYLE, name="ND Space Program")
        frame.CreateStatusBar()

        ##------ menubar
        menuBar = wx.MenuBar()
        menu = wx.Menu()
        item = menu.Append(wx.ID_NEW, "N&ew\tCtrl-N", "New")
        self.Bind(wx.EVT_MENU, self.OnNew, item)
        item = menu.Append(wx.ID_EXIT, "E&xit\tCtrl-Q", "Exit")
        self.Bind(wx.EVT_MENU, self.OnExitApp, item)
        menuBar.Append(menu, "&File")
        
        menu = wx.Menu()
        item = menu.Append(wx.NewId(), "Add Dipole", "Add dipole ...")
        self.Bind(wx.EVT_MENU, self.OnAddDipole, item)
        menuBar.Append(menu, "&Sources")

        menu = wx.Menu()
        item = menu.Append(wx.NewId(), "Add Diagonal Medium", "Add Diagonal Medium ...")
        self.Bind(wx.EVT_MENU, self.OnAddDiagonalMedium, item)
        menuBar.Append(menu, "&Media")

        menu = wx.Menu()
        item = menu.Append(wx.NewId(), "Add Boundary Condition", "Add Boundary Condition ...")
        self.Bind(wx.EVT_MENU, self.OnAddBoundaryCondition, item)
        menuBar.Append(menu, "&Boundary")

        menu = wx.Menu()
        item = menu.Append(wx.NewId(), "Add Output", "Add Output ...")
        self.Bind(wx.EVT_MENU, self.OnAddOutputPlane, item)
        menuBar.Append(menu, "&Output")


        menu = wx.Menu()
        item = menu.Append(wx.NewId(), "Simulation Parameters", "Simulation Parameters ...")
        self.Bind(wx.EVT_MENU, self.OnSetSimulationParams, item)
        item = menu.Append(wx.NewId(), "Run Simulation", "Run simulattion ...")
        self.Bind(wx.EVT_MENU, self.OnRunSimulation, item)
        item = menu.Append(wx.NewId(), "Global Variables", "Global Variables ...")
        self.Bind(wx.EVT_MENU, self.OnSetGlobalVariables, item)
        menuBar.Append(menu, "&Simulation")

        menu = wx.Menu()
        item = menu.Append(wx.NewId(), "Cube", "Draw Cube ...")
        self.Bind(wx.EVT_MENU, self.OnDrawCube, item)
        menuBar.Append(menu, "&Draw")


        frame.SetMenuBar(menuBar)

        ##----- toolbar
        toolbar = wx.ToolBar(frame, id=-1)
        toolbar.SetToolBitmapSize( (17, 17) ) # Square spacer equired for non-standard size buttons on MSW
        
        tool_save = toolbar.AddTool(wx.ID_ANY, 'Save', wx.Bitmap('./images/icons/save-16x16.png'))
        toolbar.AddSeparator()
        tool_select = toolbar.AddTool(wx.ID_ANY, 'Select', wx.Bitmap('./images/icons/select-16x16.png'), kind=wx.ITEM_RADIO)
        tool_move = toolbar.AddTool(wx.ID_ANY, 'Move', wx.Bitmap('./images/icons/move-16x16.png'), kind=wx.ITEM_RADIO)
        tool_rotate = toolbar.AddTool(wx.ID_ANY, 'Rotate', wx.Bitmap('./images/icons/rotate-16x16.png'), kind=wx.ITEM_RADIO)
        tool_zoom = toolbar.AddTool(wx.ID_ANY, 'Zoom', wx.Bitmap('./images/icons/zoom-16x16.png'), kind=wx.ITEM_RADIO)
        toolbar.AddSeparator()
        tool_toggle_grid = toolbar.AddTool(wx.ID_ANY, 'gridOn', wx.Bitmap('./images/icons/grid-16x16.png'), kind=wx.ITEM_CHECK)
        tool_toggle_axis = toolbar.AddTool(wx.ID_ANY, 'gridOn', wx.Bitmap('./images/icons/axis-16x16.png'), kind=wx.ITEM_CHECK)
        
        
        toolbar.Realize()
        
        frame.SetToolBar(toolbar)
        
        frame.Bind(wx.EVT_TOOL, self.OnToolSave, tool_save)
        frame.Bind(wx.EVT_TOOL, self.OnToolSelect, tool_select)
        frame.Bind(wx.EVT_TOOL, self.OnToolMove, tool_move)
        frame.Bind(wx.EVT_TOOL, self.OnToolRotate, tool_rotate)
        frame.Bind(wx.EVT_TOOL, self.OnToolZoom, tool_zoom)
        frame.Bind(wx.EVT_TOOL, self.OnToolToggleGrid, tool_toggle_grid)
        frame.Bind(wx.EVT_TOOL, self.OnToolToggleAxis, tool_toggle_axis)
               
        
        ##-------
        
        frame.Show(True)
        frame.Bind(wx.EVT_CLOSE, self.OnCloseFrame)

        
        self.SetTopWindow(frame)
        self.frame = frame
        
        self.quit = False
        return True
        
    def OnNew(self, evt):
        newDlg = NewSceneDialog(None, title="")
        newDlg.ShowModal()        
        
        if newDlg.simulation_params!=None:
            n_dims = newDlg.simulation_params["NumDimensions"]
            self.frame.SetNumDimensions(n_dims)
        
            self.frame.simulationParameters["Global"] = newDlg.simulation_params
        
        #newDlg.Destroy()        
        evt.Skip()

    def OnExitApp(self, evt):
        self.frame.Close(True)
        self.quit = True
        print('OnExitApp')

    def OnCloseFrame(self, evt):
        if hasattr(self, "window") and hasattr(self.window, "ShutdownDemo"):
            self.window.ShutdownDemo()
        evt.Skip()
        self.quit = True
        self.autoCaller.Stop()
        print('OnCloseFrame')
        
    def OnToolSave(self, evt):
        self.frame.LogText('Save!\n')
        evt.Skip()
        
    def OnToolSelect(self, evt):
        self.frame.selstate = SelStates.select
        self.frame.LogText('Select! '+repr(self.frame.selstate)+"\n")
        evt.Skip()

    def OnToolMove(self, evt):
        self.frame.selstate = SelStates.move
        self.frame.LogText('Move! '+repr(self.frame.selstate)+"\n")
        evt.Skip()

    def OnToolRotate(self, evt):
        self.frame.selstate = SelStates.rotate
        self.frame.LogText('Rotate! '+repr(self.frame.selstate)+"\n")
        evt.Skip()

    def OnToolZoom(self, evt):
        self.frame.selstate = SelStates.zoom
        self.frame.LogText('Zoom! '+repr(self.frame.selstate)+"\n")
        evt.Skip()
        
    def OnToolToggleGrid(self, evt):
        self.frame.gridOn = not self.frame.gridOn
        self.frame.UpdateCanvas()
        evt.Skip()

    def OnToolToggleAxis(self, evt):
        self.frame.axisOn = not self.frame.axisOn
        self.frame.UpdateCanvas()
        evt.Skip()
        
    def OnAddDipole(self, evt):
        self.frame.LogText('Add dipole!\n')
        dlg = DipoleParamsDialog(None, title="")
        dlg.ShowModal()
 
        if dlg.dipole_params!=None:
            if "Sources" not in self.frame.simulationParameters:
                self.frame.simulationParameters["Sources"] = [dlg.dipole_params]
            else:
                self.frame.simulationParameters["Sources"].append(dlg.dipole_params)
  
            x = self.StrExpressionToFloat(dlg.dipole_params['x'])
            y = self.StrExpressionToFloat(dlg.dipole_params['y'])
            z = self.StrExpressionToFloat(dlg.dipole_params['z'])
            r = np.array([x, y, z])
            pol = None
            elem = {'type':'dipole', 'v_arr': r, 'pol':pol, 'name':dlg.dipole_params['name']}
            self.frame.AddNewElement(elem)
            
            self.frame.UpdateCanvas()
            
        evt.Skip()
        
    
    def StrExpressionToNumpyFunc(self, expr, vars=sympy.Symbol('t')):
        eq = parse_expr(expr)
        for var in self.frame.globalvariables:
            eq = eq.subs(Symbol(var), self.frame.globalvariables[var])
            
        np_func = sympy.lambdify(vars, eq, modules='numpy')  
        return np_func
        

    def StrExpressionToFloat(self, expr):
        eq = parse_expr(expr)
        for var in self.frame.globalvariables:
            eq = eq.subs(Symbol(var), self.frame.globalvariables[var])
        
        return float(eq.evalf())
        

    def StrExpressionToComplex(self, expr):
        eq = parse_expr(expr)
        for var in self.frame.globalvariables:
            eq = eq.subs(Symbol(var), self.frame.globalvariables[var])
            
        return complex(eq.evalf())
        
    def OnAddDiagonalMedium(self, evt):
        self.frame.LogText('Add diagonal medium!\n')
        mediumDlg = MediumParamsDialog(None, title="")
        mediumDlg.ShowModal()
 
        
        if mediumDlg.medium_params!=None:
            if "Media" not in self.frame.simulationParameters:
                self.frame.simulationParameters["Media"] = [mediumDlg.medium_params]
            else:
                self.frame.simulationParameters["Media"].append(mediumDlg.medium_params)
        
        
        #mediumDlg.Destroy()        
        evt.Skip()

    def OnAddBoundaryCondition(self, evt):
        self.frame.LogText('Add Boundary Condition!\n')
        boundaryDlg = BoundaryParamsDialog(None, title="")
        boundaryDlg.ShowModal()
        
        if boundaryDlg.boundary_params!=None:
            if "Boundaries" not in self.frame.simulationParameters:
                self.frame.simulationParameters["Boundaries"] = boundaryDlg.boundary_params

        evt.Skip()

    def OnAddOutputPlane(self, evt):
        self.frame.LogText('Add Output Plane!\n')
        outputDlg = OutputParamsDialog(None, title="")
        outputDlg.ShowModal()

        if outputDlg.params!=None:
            if "Outputs" not in self.frame.simulationParameters:
                self.frame.simulationParameters["Outputs"] = [outputDlg.params]
            else:
                self.frame.simulationParameters["Outputs"].append(outputDlg.params)

        evt.Skip()



    def OnSetSimulationParams(self, evt):
        self.frame.LogText('Set Simulation Parameters!\n')
        simulationDlg = SimulationParamsDialog(None, title="")
        simulationDlg.ShowModal()

        if simulationDlg.simulation_params!=None:
            self.frame.simulationParameters["Simulation"] = simulationDlg.simulation_params

        evt.Skip()
        
    def OnRunSimulation(self, evt):
        self.frame.LogText('Run Simulation!\n')

        print("simultion parameters: \n", self.frame.simulationParameters)
        print("elements: \n", self.frame.elements)
        print("global variables: \n", self.frame.globalvariables)

        from Electromagnetics.FDTD import SimSpecToFDTDModel
        
        fdtdModeler = SimSpecToFDTDModel(self.frame.simulationParameters, self.frame.elements,
                            self.frame.globalvariables)
        
        fdtdModeler.validateSpecs()

        evt.Skip()
        

    def OnSetGlobalVariables(self, evt):
        self.frame.LogText('Set Global Variables!\n')

        dlg = GlobalVarsDialog(self.frame.globalvariables, None)
        dlg.ShowModal()

        if dlg.params!=None:
            var_name = dlg.params["name"]
            var_value = dlg.params["value"]
            assert var_name not in self.frame.globalvariables
            self.frame.globalvariables[var_name] = var_value

        evt.Skip()


    def OnDrawCube(self, evt):
        self.frame.LogText('Draw Cube!\n')
    
        cubeDlg = CubeParamsDialog(None, title="")
        cubeDlg.ShowModal()

        if cubeDlg.cube_params!=None:
            if "Objects" not in self.frame.simulationParameters:
                self.frame.simulationParameters["Objects"] = [cubeDlg.cube_params]
            else:
                self.frame.simulationParameters["Objects"].append(cubeDlg.cube_params)

            r0_str = cubeDlg.cube_params['r0']
            r0 = np.array([self.StrExpressionToFloat(r0_str[i]) for i in range(3)])
            r1_str = cubeDlg.cube_params['r1']
            r1 = np.array([self.StrExpressionToFloat(r1_str[i]) for i in range(3)])
            elem = {'type':'cube', 'name':cubeDlg.cube_params['name'], 'v_arr': [r0, r1], \
                        'color':cubeDlg.cube_params['color']}
            self.frame.AddNewElement(elem)
            
            self.frame.UpdateCanvas()

        evt.Skip()

    

    def SetViewOrigin(self, V_3D):
        self.V_3D = V_3D
        self.SetFrameDimensions()
        
    def SetViewAngles(self, theta, phi):
        self.V_theta = theta
        self.V_phi = phi
        self.SetFrameDimensions()
        
    def SetGuideGridSize(self, dW):
        self.dW_3D = dw
        self.SetFrameDimensions()
        
    def SetFrameDimensions(self):
        self.frame.UpdateDimensions(W_3D=self.W_3D, dW_3D=self.dW_3D, 
            V_3D=self.V_3D, V_theta=self.V_theta, V_phi=self.V_phi)
        
    def deleteGUIItems(self):
        ## deletes window components (frame,...)
        self.quit = True
        del self.frame
        
    def SetCommandLoop(self, time_ms=1000):
        self.time_ms = time_ms
        self.autoCaller = wx.CallLater(time_ms, self.CommandLoop)
        self.autoCaller.Start()
        
        
    def CommandLoop(self):
        if not self.quit:
            if self.queue != None:
                #print('queue: ', self.queue)
                if not self.queue.empty():
                    comm = self.queue.get()
                    self.ProcessCommand(comm)
                    self.queue.task_done()
            self.autoCaller.Restart()
        return

    def ProcessCommand(self, comm):
        if comm[0]=='deleteAllElems':
            self.frame.elements = []
        elif comm[0]=='add elem':
            if comm[1] == 'cube':
                elem = {'type':'cube', 'v_arr': comm[2]}
                self.frame.AddNewElement(elem)
            elif comm[1] == 'triangle':
                elem = {'type':'triangle', 'v_arr': comm[2]}
                self.frame.AddNewElement(elem)
            elif comm[1] == 'quad':
                elem = {'type':'quad', 'v_arr': comm[2]}
                self.frame.AddNewElement(elem)
            self.frame.UpdateCanvas()
        elif comm[0]=="logSimulationParameters":
            self.LogSimulationParams()
        elif comm[0]=="logtext":
            text = comm[1]
            self.frame.LogText(text)
            #print("test command received!")
        return
        
        
    def LogSimulationParams(self):
        self.frame.LogText('Simulation parameters: {} \n '.format(self.frame.simulationParameters))
        self.frame.LogText('Elements: {} \n '.format(self.frame.elements))
    


class SelStates(Enum): 
    select = 1
    move = 2
    rotate = 3
    zoom = 4
    
    
class GLFrame(wx.Frame):

    def __init__(self, parent, id, title, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=wx.DEFAULT_FRAME_STYLE,
                 name='frame', W_3D=np.array([1.0, 1.0, 1.0])):

        sys.excepthook = MyExceptionHook

        ##-- stencil index
        self.stencil_index_next = 1
        self.elements = []

        ##-- dimensions and units
        self.W_3D = W_3D          # world dimensions
        self.dW_3D = W_3D/10.0
        self.V_3D = np.array([0.0, 0.0, 0.0])        # View origin
        self.V_phi = 20.0       # view angle phi
        self.V_theta = 30.0     # view angle theta

        self.queue = None       ## can also be used to send commands to the parent window

        self.selstate = SelStates.select
        self.move_started = False
        self.rotate_started = False
        self.zoom_started = False
        
        self.gridOn = True
        self.axisOn = True
        
        self.SetInitPerspParams()

        #
        # Forcing a specific style on the window.
        #   Should this include styles passed?
        style = wx.DEFAULT_FRAME_STYLE | wx.NO_FULL_REPAINT_ON_RESIZE

        super(GLFrame, self).__init__(parent, id, title, pos, size, style, name)

        self.GLinitialized = False
        attribList = (glcanvas.WX_GL_RGBA, # RGBA
                      glcanvas.WX_GL_DOUBLEBUFFER, # Double Buffered
                      glcanvas.WX_GL_DEPTH_SIZE, 24) # 24 bit

        #
        # Create the canvas
        self.canvas = glcanvas.GLCanvas(self, attribList=attribList)
        self.context = glcanvas.GLContext(self.canvas)
        
        # create cursors
        self.SetupCursers()
        
        #
        # Set the event handlers.
        self.canvas.Bind(wx.EVT_ERASE_BACKGROUND, self.processEraseBackgroundEvent)
        self.canvas.Bind(wx.EVT_SIZE, self.processSizeEvent)
        self.canvas.Bind(wx.EVT_PAINT, self.processPaintEvent)
        self.canvas.Bind(wx.EVT_LEFT_DOWN, self.OnCanvasLeftMouseDown)
        self.canvas.Bind(wx.EVT_LEFT_UP, self.OnCanvasLeftMouseUp)
        self.canvas.Bind(wx.EVT_MOTION, self.OnCanvasMouseMove)
        self.canvas.Bind(wx.EVT_RIGHT_DOWN, self.OnCanvasRightMouseDown)
        
        ## Sizer
        box = wx.BoxSizer(wx.VERTICAL)
        #box.Add((20, 30))

        self.canvas.SetMinSize((500, 400))
        box.Add(self.canvas, 1, wx.ALIGN_CENTER|wx.ALL|wx.EXPAND, 2)

        ##--- textbox
        self.textbox = wx.TextCtrl(self, size=(300, 50), style=wx.TE_MULTILINE)        
        box.Add(self.textbox, 0, wx.ALIGN_CENTER|wx.BOTTOM|wx.EXPAND, 5)

        self.SetAutoLayout(True)
        self.SetSizer(box)
                
                
        ##--- guide grid
        el_guidegrid = {'type':'guidegrid'}
        self.AddNewElement(el_guidegrid)
        
        ##--- default colors
        self.color_ambient_default = (0.2, 0.2, 0.2, 1.0)
        self.color_emission_default = (0.0, 0.0, 0.0, 1.0)
        self.color_specular_default = (1.0, 1.0, 1.0, 1.0)
        self.shininess_default = 60
        
        ##---
        self.n_dims = 3
        self.projection = "P"   ## "P":perspective   "O":orthogonal
  
        ##-- simulation parameters
        self.simulationParameters = {}
        
        self.globalvariables = {}
        
        
        
    def AddNewElement(self, elem, assign_stencil_index=True):
        elem['STI'] = 0
        if assign_stencil_index:
            elem['STI'] = self.stencil_index_next
            self.stencil_index_next += 1
        
        elem['is_selected'] = False
        self.elements.append(elem)
    
  
    def SetNumDimensions(self, n_dims):
        assert 1<=n_dims<=3
        self.n_dims = n_dims     
        if n_dims==3:
           self.projection = "P"
        elif n_dims<=2:
           self.projection = "O"
           
        self.LogText('Number of dimensions Set to {} \n'.format(n_dims))
        
        size_GL = self.GetGLExtents()
        self.SetProjection(size_GL.width, size_GL.height)
        self.UpdateCanvas()
    
    
    def LogText(self, text):
        self.textbox.AppendText(text)
    
    #
    # Canvas Proxy Methods

    def GetGLExtents(self):
        """Get the extents of the OpenGL canvas."""
        return self.canvas.GetClientSize()

    def SwapBuffers(self):
        """Swap the OpenGL buffers."""
        self.canvas.SwapBuffers()

    #
    # wxPython Window Handlers

    def processEraseBackgroundEvent(self, event):
        """Process the erase background event."""
        pass # Do nothing, to avoid flashing on MSWin

    def processSizeEvent(self, event):
        #size_frame = self.GetSize()
        #X_frame, Y_frame = size_frame
        ##---
        #self.textbox.SetWidth(X_frame)
        
        """Process the resize event."""
        size_GL = self.GetGLExtents()
        #self.canvas.SetCurrent(self.context)
        self.OnReshape(size_GL.width, size_GL.height)
        self.canvas.Refresh(False)
        self.canvas.Update()
        event.Skip()

    def processPaintEvent(self, event):
        """Process the drawing event."""
        #self.Show()
        self.canvas.SetCurrent(self.context)

        # initialize OpenGL
        if not self.GLinitialized:
            self.OnInitGL()
            self.GLinitialized = True

        self.OnDraw()
        self.canvas.Refresh(False)
        self.canvas.Update()
        event.Skip()
        
    def OnCanvasLeftMouseDown(self, event):
        size_GL = self.GetGLExtents()
        x = event.GetX()
        y = size_GL.height - 1 - event.GetY()
        #self.LogText('Canvas: Left mouse down! (x={}, y={})\n'.format(x, y))
        
        if self.selstate==SelStates.select:
            stencil_index = glReadPixels(x, y, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT)
            self.LogText('Stencil index: {}\n'.format(stencil_index))
            stencil_index = stencil_index[0][0]
            if stencil_index>0:
                elem = self.FindTheClickedObject(stencil_index)
                if elem['is_selected']==False:
                    self.UnselectAll(False)
                    elem['is_selected'] = True
                    self.UpdateCanvas()
            else:
                self.UnselectAll()
        elif self.selstate==SelStates.move:
            self.x_sel = x
            self.y_sel = y
            self.move_started = True
            self.ChangeCurser('hand')
        elif self.selstate==SelStates.rotate:
            self.x_sel = x
            self.y_sel = y
            self.rotate_started = True
            self.ChangeCurser('rotate')
        elif self.selstate==SelStates.zoom:
            self.x_sel = x
            self.y_sel = y
            self.zoom_started = True
            self.ChangeCurser('magnifier')
        
        event.Skip()

    def OnCanvasLeftMouseUp(self, event):
        if self.selstate==SelStates.select:
            pass
        elif self.selstate==SelStates.move:
            self.move_started = False
            self.UpdateCanvas()
            self.ChangeCurser('arrow')
        elif self.selstate==SelStates.rotate:
            self.rotate_started = False
            self.UpdateCanvas()
            self.ChangeCurser('arrow')
        elif self.selstate==SelStates.zoom:
            self.zoom_started = False
            self.UpdateCanvas()
            self.ChangeCurser('arrow')
        
        event.Skip()

    def OnCanvasMouseMove(self, event):
        size_GL = self.GetGLExtents()
        x = event.GetX()
        y = size_GL.height - 1 - event.GetY()

        if self.selstate==SelStates.select:
            pass
        elif self.selstate==SelStates.move:
            if self.move_started:
                dx = x - self.x_sel
                dy = y - self.y_sel
                X, Y, Z = self.W_3D[0], self.W_3D[1], self.W_3D[2]
                X_ratio = float(X)/size_GL.width
                Y_ratio = float(Y)/size_GL.height
                self.P_Tx += dx * X_ratio
                self.P_Ty += dy * Y_ratio 
                #self.LogText('dx={}, dy={} '.format(dx, dy))
                #self.queue.put(['logtext', 'dx:{} dy:{} '.format(dx, dy)])
                self.SetProjection(size_GL.width, size_GL.height)
                self.UpdateCanvas()
                self.x_sel = x
                self.y_sel = y
        elif self.selstate==SelStates.rotate:
            if self.rotate_started:
                dx = x - self.x_sel
                dy = y - self.y_sel
                X, Y, Z = self.W_3D[0], self.W_3D[1], self.W_3D[2]
                X_ratio = float(X)/size_GL.width*25
                Y_ratio = float(Y)/size_GL.height*25
                self.P_Rz += dx * X_ratio
                self.P_Rx -= dy * Y_ratio 
                self.SetProjection(size_GL.width, size_GL.height)
                self.UpdateCanvas()
                self.x_sel = x
                self.y_sel = y
        elif self.selstate==SelStates.zoom:
            if self.zoom_started:
                dx = x - self.x_sel
                dy = y - self.y_sel
                X, Y, Z = self.W_3D[0], self.W_3D[1], self.W_3D[2]
                #X_ratio = float(X)/size_GL.width*25
                Y_ratio = float(Y)/size_GL.height*15
                self.P_angle -= dy * Y_ratio
                if self.P_angle<10.0:
                    self.P_angle = 10.0
                if self.P_angle>120.0:
                    self.P_angle = 120.0
                self.SetProjection(size_GL.width, size_GL.height)
                self.UpdateCanvas()
                self.x_sel = x
                self.y_sel = y
        
        event.Skip()



    def OnCanvasRightMouseDown(self, event):
        size_GL = self.GetGLExtents()
        x = event.GetX()
        y = size_GL.height - 1 - event.GetY()
        #self.LogText('Canvas: Left mouse down! (x={}, y={})\n'.format(x, y))
        
        if self.selstate==SelStates.select:
            stencil_index = glReadPixels(x, y, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT)[0][0]
            #self.LogText('Stencil index: {}\n'.format(stencil_index))
            if stencil_index>0:
                elem = self.FindTheClickedObject(stencil_index)
                if elem['is_selected']==True:
                    if elem['type'] in ["cube", "sphere"]:
                        obj_type = "VOL_OBJ"
                        if elem['type']=="cube":
                            ## assert it is not rotated
                            is_rotated = False
                            if not is_rotated:
                                obj_type = "CUBE_NOROT"
                        menu = PopupMenu_Object(obj_type=obj_type)
                        self.PopupMenu(menu, event.GetPosition())
                        menu.Destroy()
                        
                        if menu.request == menu.request_types.volobj_assign_material:
                            self.LogText("Assign material\n")
                            medium_list = self.GetMediumList()
                            assignMedDlg = AssignMediumDialog(medium_list, None, title="")
                            assignMedDlg.ShowModal()
                            
                            params = assignMedDlg.params
                            if params != None:
                                elem_name = elem['name']
                                assert "Objects" in self.simulationParameters
                                objects = self.simulationParameters["Objects"]
                                for obj in objects:
                                    if obj["name"] == elem_name:
                                        obj['material'] = params['name']
                        elif menu.request == menu.request_types.volobj_show_hide:
                            self.LogText("Flip Show/Hide\n")
                            if 'show' in elem:
                                elem['show'] = not elem['show']
                            else:
                                elem['show'] = False
                            self.UpdateCanvas()
                        elif menu.request == menu.request_types.cube_set_as_simulation_box:
                            self.LogText("Set as simulation box.\n")
                            self.simulationParameters["SimBox"] = elem["name"]
                            
                    if elem['type'] in ["dipole"]:
                        menu = PopupMenu_Source(obj_type="dipole")
                        self.PopupMenu(menu, event.GetPosition())
                        menu.Destroy()
                        
            else:
                self.UnselectAll()

                menu = PopupMenu_Object()
                self.PopupMenu(menu, event.GetPosition())
                menu.Destroy()

                if menu.request == menu.request_types.none_set_units:
                    self.LogText("Set Units\n")
                
        event.Skip()



    ## --- 
    def GetMediumList(self):
        medium_list = []
        if "Media" in self.simulationParameters:
            media = self.simulationParameters['Media']
            for med in media:
                medium_list.append(med["name"])
        return medium_list


    ## cursors 

    def SetupCursers(self, model=None):
        if model==None or model=='arrow':
            self.cursor_arrow = wx.Cursor(wx.CURSOR_ARROW)
        if model==None or model=='hand':
            image = wx.Image(r'./images/cursors/hand-16x16.png', wx.BITMAP_TYPE_PNG)
            self.cursor_hand = wx.Cursor(image)
        if model==None or model=='magnifier':
            image = wx.Image(r'./images/cursors/magnifier-16x16.png', wx.BITMAP_TYPE_PNG)
            self.cursor_magnifier = wx.Cursor(image)
        if model==None or model=='move':
            image = wx.Image(r'./images/cursors/move-16x16.png', wx.BITMAP_TYPE_PNG)
            self.cursor_move = wx.Cursor(image)
        if model==None or model=='rotate':
            image = wx.Image(r'./images/cursors/rotate-16x16.png', wx.BITMAP_TYPE_PNG)
            self.cursor_rotate = wx.Cursor(image)

    def ChangeCurser(self, model='arrow'):
        if model=='arrow':
            self.canvas.SetCursor(self.cursor_arrow)
        elif model=='hand':
            self.canvas.SetCursor(self.cursor_hand)
        elif model=='magnifier':
            self.canvas.SetCursor(self.cursor_magnifier)
        elif model=='move':
            self.canvas.SetCursor(self.cursor_move)
        elif model=='rotate':
            self.canvas.SetCursor(self.cursor_rotate)

    #
    # GLFrame OpenGL Event Handlers
    
    def UpdateCanvas(self):
        self.canvas.Refresh(True)
        self.canvas.Update()
        #self.canvas.SwapBuffers()
        

    def OnInitGL(self):
        """Initialize OpenGL for use in the window."""
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
        glEnable(GL_DEPTH_TEST)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT)
        glClearColor(1, 1, 1, 1)

        glEnable(GL_STENCIL_TEST)
        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE)

        glEnable(GL_LIGHTING)

        glEnable(GL_LIGHT0)
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE)
        
        ambient = [ 0.2, 0.2, 0.2, 1.0 ]
        diffuse = [ 1.0, 1.0, 1.0, 1.0 ]
        position = [ 10.0, 10.0, 10.0, 1.0]
        glLightfv(GL_LIGHT0, GL_AMBIENT, ambient)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse)
        glLightfv(GL_LIGHT0, GL_POSITION, position)
        
        glEnable(GL_LIGHT1)
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE)
        
        ambient = [ 0.2, 0.2, 0.2, 1.0 ]
        diffuse = [ 1.0, 1.0, 1.0, 1.0 ]
        position = [ -10.0, -10.0, -10.0, 1.0]
        glLightfv(GL_LIGHT1, GL_AMBIENT, ambient)
        glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse)
        glLightfv(GL_LIGHT1, GL_POSITION, position)


    def OnReshape(self, width, height):
        """Reshape the OpenGL viewport based on the dimensions of the window."""
        glViewport(0, 0, width, height)

        self.SetProjection(width, height)
        
        
    def SetInitPerspParams(self):
        X, Y, Z = self.W_3D[0], self.W_3D[1], self.W_3D[2]
        XYZ_Max = max(X, Y, Z)
        XYZ_Min = min(X, Y, Z)

        self.P_angle = 60.0
        self.P_near = 0.1*XYZ_Min
        self.P_far = 6.0*XYZ_Max
        
        self.P_Tx = 0.0
        self.P_Ty = 0.0
        self.P_Tz = -2.0*XYZ_Max
        
        self.P_Rx = -60.0
        self.P_Ry = 0.0
        self.P_Rz = -20.0
                

    def SetPerspective(self, width, height):

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(self.P_angle, float(width)/float(height), self.P_near, self.P_far)
        glTranslated(self.P_Tx, self.P_Ty, self.P_Tz)
        glRotated(self.P_Rx, 1.0, 0.0, 0.0)
        glRotated(self.P_Ry, 0.0, 1.0, 0.0)
        glRotated(self.P_Rz, 0.0, 0.0, 1.0)
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

    
    def SetParallelProjection(self, width, height):
        X, Y , Z = self.W_3D
    
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if X/Y <= width/height:
            ratio = width/height
            gluOrtho2D(-Y/2*ratio, Y/2*ratio, -Y/2, Y/2)
        else:
            ratio = height/width
            gluOrtho2D(-X/2, X/2, -X/2*ratio, X/2*ratio)
        glTranslated(self.P_Tx, self.P_Ty, 0.0)

        """
        glTranslated(self.P_Tx, self.P_Ty, self.P_Tz)
        glRotated(self.P_Rx, 1.0, 0.0, 0.0)
        glRotated(self.P_Ry, 0.0, 1.0, 0.0)
        glRotated(self.P_Rz, 0.0, 0.0, 1.0)
        """
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
    
    
    def SetProjection(self, width, height):
        if self.projection=="P":
            self.SetPerspective(width, height)
        else:
            self.SetParallelProjection(width, height)
        

    def OnDraw(self, *args, **kwargs):
        """Draw the window.
        """
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT)
        
        for elem in self.elements:
            glStencilFunc(GL_ALWAYS, elem['STI'], -1)
            if elem['type']=='guidegrid':
                self.DrawGuideGrid()
            elif elem['type']=='triangle':
                r0, r1, r2 = elem['v_arr']
                self.DrawTriangle(r0, r1, r2)
            elif elem['type']=='cube':
                r0, r1 = elem['v_arr']
                color = (1.0, 0.0, 0.0, 1.0)
                if 'color' in elem:
                    color = elem['color']
                is_selected = elem['is_selected']
                show = True
                if 'show' in elem:
                    show = elem['show']
                self.DrawCube(r0, r1, color=color, is_selected=is_selected, show=show)
            elif elem['type']=='dipole':
                r0 = elem['v_arr']
                pol = elem['pol']
                color = (1.0, 0.0, 0.0, 1.0)
                if 'color' in elem:
                    color = elem['color']
                is_selected = elem['is_selected']
                self.DrawDipolePoint(r0, color=color, is_selected=is_selected)
            elif elem['type']=='quad':
                r0, r1, r2, r3 = elem['v_arr']
                self.DrawQuad(r0, r1, r2, r3)
        
        self.SwapBuffers()
                    

    def UpdateDimensions(self, W_3D, dW_3D, V_3D, V_theta, V_phi):
        self.W_3D = W_3D
        self.dW_3D = dW_3D
        self.V_3D = V_3D
        self.V_theta = V_theta
        self.V_phi = V_phi
        

    def DrawGuideGrid(self):
        dx = self.dW_3D[0]
        dy = self.dW_3D[1]
        X = self.W_3D[0]
        Y = self.W_3D[1]
        Z = self.W_3D[2]
        
        glLineWidth(1.5)
        glBegin(GL_LINES)

        if self.axisOn:
            color_emission = [1.0, 0.0, 0.0, 1.0]
            glColor4d(*tuple(color_emission))
            glMaterialfv(GL_FRONT, GL_EMISSION, color_emission)
            glMaterialfv(GL_FRONT, GL_DIFFUSE, color_emission)
            glMaterialfv(GL_FRONT, GL_AMBIENT, color_emission)
            glMaterialfv(GL_FRONT, GL_SPECULAR, color_emission)
            glMaterialf(GL_FRONT, GL_SHININESS, 120)
            glVertex3d(0.0, 0.0, 0.0)
            glVertex3d(+X/2, 0.0, 0.0)
            
            color_emission = [0.0, 1.0, 0.0, 1.0]
            glColor4d(*tuple(color_emission))
            glMaterialfv(GL_FRONT, GL_EMISSION, color_emission)
            glMaterialfv(GL_FRONT, GL_DIFFUSE, color_emission)
            glMaterialfv(GL_FRONT, GL_AMBIENT, color_emission)
            glMaterialfv(GL_FRONT, GL_SPECULAR, color_emission)
            glMaterialf(GL_FRONT, GL_SHININESS, 120)
            glVertex3d(0.0, 0.0, 0.0)
            glVertex3d(0.0, Y/2, 0.0)

            color_emission = [0.0, 0.0, 1.0, 1.0]
            glMaterialfv(GL_FRONT, GL_EMISSION, color_emission)
            glMaterialfv(GL_FRONT, GL_DIFFUSE, color_emission)
            glMaterialfv(GL_FRONT, GL_AMBIENT, color_emission)
            glMaterialfv(GL_FRONT, GL_SPECULAR, color_emission)
            glMaterialf(GL_FRONT, GL_SHININESS, 120)
            glVertex3d(0.0, 0.0, 0.0)
            glVertex3d(0.0, 0.0, Z/2)
        glEnd()


        glLineWidth(1.0)
        glBegin(GL_LINES)
        if self.gridOn:
            color = [0.0, 0.0, 0.0, 1.0]
            glColor4d(*tuple(color))
            glMaterialfv(GL_FRONT, GL_AMBIENT, color)
            glMaterialfv(GL_FRONT, GL_DIFFUSE, color)
            glMaterialfv(GL_FRONT, GL_SPECULAR, color)
            glMaterialfv(GL_FRONT, GL_EMISSION, color)
            glMaterialf(GL_FRONT, GL_SHININESS, 0.0)
            for i in range(-5, 6):
                glVertex3d(i*dx, -Y/2, 0.0)
                glVertex3d(i*dx, +Y/2, 0.0)
            for i in range(-5, 6):
                glVertex3d(-X/2, i*dy, 0.0)
                glVertex3d(+X/2, i*dy, 0.0)
        
        glEnd()


    def DrawTriangle(self, r0, r1, r2, color=None, color_ambient=None, color_diffuse=None, color_specular=None, color_emission=None, shininess=None):
        # Drawing an example triangle in the middle of the screen
        glBegin(GL_TRIANGLES)
        if color==None:
            color = (random.random(), random.random(), random.random(), 1.0)
        if color_ambient==None:
            color_ambient = self.color_ambient_default
        if color_diffuse==None:
            color_diffuse = color
        if color_specular==None:
            color_specular = self.color_specular_default
        if color_emission==None:
            color_emission = self.color_emission_default
        if shininess==None:
            shininess = self.shininess_default


        glColor4d(*tuple(color))
        glMaterialfv(GL_FRONT, GL_AMBIENT, color_ambient)
        glMaterialfv(GL_FRONT, GL_DIFFUSE, color_diffuse)
        glMaterialfv(GL_FRONT, GL_SPECULAR, color_specular)
        glMaterialfv(GL_FRONT, GL_EMISSION, color_emission)
        glMaterialf(GL_FRONT, GL_SHININESS, shininess)

        glVertex3d(r0[0], r0[1], r0[2])
        glVertex3d(r1[0], r1[1], r1[2])
        glVertex3d(r2[0], r2[1], r2[2])
        glEnd()
    
    def DrawCube(self, r0, r1, color=[1., 0.0, 0.0, 1.0], color_ambient=None, color_diffuse=None, color_specular=None, color_emission=None, shininess=None, 
                is_selected=False, selection_color=None, show=True):
        x0, y0, z0 = r0
        x1, y1, z1 = r1

        if color==None:
            color = (random.random(), random.random(), random.random(), 1.0)
        if color_ambient==None:
            color_ambient = self.color_ambient_default
        if color_diffuse==None:
            color_diffuse = color
        if color_specular==None:
            color_specular = self.color_specular_default
        if color_emission==None:
            color_emission = self.color_emission_default
        if shininess==None:
            shininess = self.shininess_default

        if show:
            # position viewer
            glMatrixMode(GL_MODELVIEW)
            glLoadIdentity()
            #glRotated(0.0, 1.0, 0.0, 0.0)
            #glTranslated(0.0, 0.0, 0.0)

            # draw six faces of a cube
            glBegin(GL_QUADS)
            glNormal3d(0.0, 0.0, 1.0)
            glColor4d(*tuple(color))
            glMaterialfv(GL_FRONT, GL_AMBIENT, color_ambient)
            glMaterialfv(GL_FRONT, GL_DIFFUSE, color_diffuse)
            glMaterialfv(GL_FRONT, GL_SPECULAR, color_specular)
            glMaterialfv(GL_FRONT, GL_EMISSION, color_emission)
            glMaterialf(GL_FRONT, GL_SHININESS, shininess)
            glVertex3d(x0, y0, z1)
            glVertex3d(x1, y0, z1)
            glVertex3d(x1, y1, z1)
            glVertex3d(x0, y1, z1)

            glNormal3d( 0.0, 0.0, -1.0)
            glVertex3d(x0, y0, z0)
            glVertex3d(x0, y1, z0)
            glVertex3d(x1, y1, z0)
            glVertex3d(x1, y0, z0)

            glNormal3d( 0.0, 1.0, 0.0)
            glVertex3d(x0, y1, z0)
            glVertex3d(x0, y1, z1)
            glVertex3d(x1, y1, z1)
            glVertex3d(x1, y1, z0)

            glNormal3d( 0.0, -1.0, 0.0)
            glVertex3d(x0, y0, z0)
            glVertex3d(x0, y0, z1)
            glVertex3d(x1, y0, z1)
            glVertex3d(x1, y0, z0)
            
            glNormal3d(1.0, 0.0, 0.0)
            glVertex3d(x1, y0, z0)
            glVertex3d(x1, y0, z1)
            glVertex3d(x1, y1, z1)
            glVertex3d(x1, y1, z0)

            glNormal3d(-1.0, 0.0, 0.0)
            glVertex3d(x0, y0, z0)
            glVertex3d(x0, y0, z1)
            glVertex3d(x0, y1, z1)
            glVertex3d(x0, y1, z0)
            glEnd()
        
        if is_selected or not show:
            if selection_color==None:
                selection_color = 1.0 - np.array(color)
                selection_color[3] = 1.0
            if not is_selected:
                selection_color = np.array(color)
                selection_color[3] = 1.0
            
        
            ##self.LogText('selection_color : {} \n'.format(selection_color))
            glLineWidth(8.0)
            if not is_selected:
                glLineWidth(6.0)

            glBegin(GL_LINES)
            glColor4d(*tuple(selection_color))
            glMaterialfv(GL_FRONT, GL_EMISSION, selection_color)
            glMaterialfv(GL_FRONT, GL_DIFFUSE, selection_color)
            glMaterialfv(GL_FRONT, GL_AMBIENT, selection_color)
            glMaterialfv(GL_FRONT, GL_SPECULAR, selection_color)
            glMaterialf(GL_FRONT, GL_SHININESS, 120)

            glVertex3d(x0, y0, z0)
            glVertex3d(x0, y0, z1)
            
            glVertex3d(x0, y0, z0)
            glVertex3d(x0, y1, z0)

            glVertex3d(x0, y1, z0)
            glVertex3d(x0, y1, z1)

            glVertex3d(x0, y0, z1)
            glVertex3d(x0, y1, z1)

            glVertex3d(x1, y0, z0)
            glVertex3d(x1, y0, z1)
            
            glVertex3d(x1, y0, z0)
            glVertex3d(x1, y1, z0)

            glVertex3d(x1, y1, z0)
            glVertex3d(x1, y1, z1)

            glVertex3d(x1, y0, z1)
            glVertex3d(x1, y1, z1)

            glVertex3d(x0, y0, z0)
            glVertex3d(x1, y0, z0)
            
            glVertex3d(x0, y0, z1)
            glVertex3d(x1, y0, z1)
            
            glVertex3d(x0, y1, z0)
            glVertex3d(x1, y1, z0)

            glVertex3d(x0, y1, z1)
            glVertex3d(x1, y1, z1)
            
            glEnd()

        glFlush()

    
    def DrawCylinder(self, base_center, base_radius, top_radius, height, slices, stacks, 
                color=[1., 0.0, 0.0, 1.0], color_ambient=None, color_diffuse=None, 
                color_specular=None, color_emission=None,      shininess=None, is_selected=False, selection_color=None):

        pass    
    

    def DrawDipoleVec(self, r, pol, color):
        """ r: position 
            pol: polarization vector (3 vector)
        """
        pol_norm = np.linalg.norm(pol)
        assert pol_norm != 0
        pol /= pol_norm*10
        x0, y0, z0 = r - pol
        x1, y1, z1 = r + pol
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        
        glLineWidth(3.0)

        glBegin(GL_LINES)
        glColor4d(*tuple(color))

        glVertex3d(x0, y0, z0)
        glVertex3d(x1, y1, z1)
        glEnd()

        glFlush()
        
        
    def DrawDipolePoint(self, r, color, is_selected=False, selection_color=None):
        """ r: position 
            pol: polarization vector (3 vector)
        """
        x0, y0, z0 = r 
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        
        glPointSize(10.0)

        glBegin(GL_POINTS)
        glColor4d(*tuple(color))
        glMaterialfv(GL_FRONT, GL_EMISSION, color)
        glMaterialfv(GL_FRONT, GL_DIFFUSE, color)
        glMaterialfv(GL_FRONT, GL_AMBIENT, color)
        glMaterialfv(GL_FRONT, GL_SPECULAR, color)
        glMaterialf(GL_FRONT, GL_SHININESS, 120)
        if is_selected:
            if selection_color==None:
                selection_color = 1.0 - np.array(color)
                selection_color[3] = 1.0
            glColor4d(*tuple(selection_color))
            glMaterialfv(GL_FRONT, GL_EMISSION, selection_color)
            glMaterialfv(GL_FRONT, GL_DIFFUSE, selection_color)
            glMaterialfv(GL_FRONT, GL_AMBIENT, selection_color)
            glMaterialfv(GL_FRONT, GL_SPECULAR, selection_color)
            glMaterialf(GL_FRONT, GL_SHININESS, 120)


        glVertex3d(x0, y0, z0)
        glEnd()

        glFlush()



    def DrawQuad(self, r0, r1, r2, r3, normal=None, color=None):
        x0, y0, z0 = r0
        x1, y1, z1 = r1
        x2, y2, z2 = r2
        x3, y3, z3 = r3

        # position viewer
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

        # draw six faces of a cube
        glBegin(GL_QUADS)
        if color==None:
            glColor3d(random.random(), random.random(), random.random())
        else:
            glColor3d(color[0], color[1], color[2])
        if normal==None:
            glNormal3d(0.0, 0.0, 1.0)
        else:
            glNormal3d(normal[0], normal[1], normal[2])
            
        glVertex3d(x0, y0, z0)
        glVertex3d(x1, y1, z1)
        glVertex3d(x2, y2, z2)
        glVertex3d(x3, y3, z3)
        glEnd()



    ##--- selection
    
    def UnselectAll(self, set_update=True):
        update = False
        for elem in self.elements:
            if elem['is_selected']: 
                elem['is_selected'] = False
                update = True
        if update and set_update:
            self.UpdateCanvas()
        
    
    def FindTheClickedObject(self, stencil_index):
        for elem in self.elements:
            if 'STI' in elem:
                if elem['STI'] == stencil_index:
                    return elem
        return None




##------------------- Dialog Boxes

class NewSceneDialog(wx.Dialog):
    
    def __init__(self, *args, **kw):
        super(NewSceneDialog, self).__init__(*args, **kw) 
            
        self.InitUI()
        self.SetSize((400, 300))
        self.SetTitle("New Simulation Parameters")
        
        self.simulation_params = None
        
        
    def InitUI(self):
        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        toplevbox = wx.StaticBox(pnl, label='Problem Dimensions')
        toplevsizer = wx.StaticBoxSizer(toplevbox, orient=wx.VERTICAL)     
        _3d_ = wx.RadioButton(pnl, label='3D', style=wx.RB_GROUP)   
        _2d_ = wx.RadioButton(pnl, label='2D')   
        _1d_ = wx.RadioButton(pnl, label='1D')   
        toplevsizer.Add(_3d_)
        toplevsizer.Add(_2d_)
        toplevsizer.Add(_1d_)
        
        pnl.SetSizer(toplevsizer)
         
        self._3d_ = _3d_
        self._2d_ = _2d_
        self._1d_ = _1d_
        
        ##buttons
        butsizer = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        butsizer.Add(okButton)
        butsizer.Add(closeButton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, 
            flag=wx.ALL|wx.EXPAND, border=5)
        vbox.Add(butsizer, 
            flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        
        okButton.Bind(wx.EVT_BUTTON, self.OnOk)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCancel)
        
    def OnClose(self, e):
        self.Destroy()
        
    def OnCancel(self, e):
        self.Destroy()

    def OnOk(self, e):
        simulation_params = {}
        if self._3d_.GetValue()==True:
            simulation_params['NumDimensions'] = 3
        elif self._2d_.GetValue()==True:
            simulation_params['NumDimensions'] = 2
        elif self._1d_.GetValue()==True:
            simulation_params['NumDimensions'] = 1
        
        print("simulation_params: ", simulation_params)
        self.simulation_params = simulation_params
        self.Destroy()


class DipoleParamsDialog(wx.Dialog):
    
    def __init__(self, *args, **kw):
        super(DipoleParamsDialog, self).__init__(*args, **kw) 
            
        self.InitUI()
        self.SetSize((400, 400))
        self.SetTitle("Set Dipole Parameters")
        
        self.dipole_params = None
        
        
    def InitUI(self):

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        toplevbox = wx.StaticBox(pnl, label='Dipole Parameters')
        toplevsizer = wx.StaticBoxSizer(toplevbox, orient=wx.VERTICAL)     
        
        ##---- name
        name_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        name_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'name:  '))
        name =  wx.TextCtrl(pnl, wx.NewId(), value=nameGeneraor.GenerateNewName('dipole'), size=wx.Size(200, 30))
        name_sizer.Add(name)
               
        toplevsizer.Add(name_sizer)
        
        self.name = name
        
        
        ##electric/magnetic
        is_electric = wx.RadioButton(pnl, label='Electric', style=wx.RB_GROUP)   
        is_magnetic = wx.RadioButton(pnl, label='Magnetic')   
        toplevsizer.Add(is_electric)
        toplevsizer.Add(is_magnetic)
        
        self.is_electric = is_electric
        self.is_magnetic = is_magnetic
        
        ##location
        locsizer = wx.BoxSizer(wx.HORIZONTAL)        
        locsizer.Add(wx.StaticText(pnl, wx.NewId(), label='Location:      '))
        locsizer_2 = wx.BoxSizer(wx.VERTICAL)        
        loc_x_txt = wx.TextCtrl(pnl, wx.NewId(), value="0.0")
        loc_y_txt = wx.TextCtrl(pnl, wx.NewId(), value="0.0")
        loc_z_txt = wx.TextCtrl(pnl, wx.NewId(), value="0.0")
        locsizer_2.Add(loc_x_txt, flag=wx.LEFT, border=5)
        locsizer_2.Add(loc_y_txt, flag=wx.LEFT, border=5)
        locsizer_2.Add(loc_z_txt, flag=wx.LEFT, border=5)
        locsizer.Add(locsizer_2, flag=wx.LEFT, border=5)
        toplevsizer.Add(locsizer)
        self.loc_x_txt = loc_x_txt
        self.loc_y_txt = loc_y_txt
        self.loc_z_txt = loc_z_txt
        
        ##polarization
        polsizer = wx.BoxSizer(wx.HORIZONTAL)        
        polsizer.Add(wx.StaticText(pnl, wx.NewId(), label='Amplitude:  '))
        polsizer_2 = wx.BoxSizer(wx.VERTICAL) 
        pol_x_txt = wx.TextCtrl(pnl, wx.NewId(), value="1.5*sin(2*pi*t)", size=wx.Size(300, 30))
        pol_y_txt = wx.TextCtrl(pnl, wx.NewId(), value="0.0", size=wx.Size(300, 30))
        pol_z_txt = wx.TextCtrl(pnl, wx.NewId(), value="0.0", size=wx.Size(300, 30))
        polsizer_2.Add(pol_x_txt, flag=wx.LEFT, border=5)
        polsizer_2.Add(pol_y_txt, flag=wx.LEFT, border=5)
        polsizer_2.Add(pol_z_txt, flag=wx.LEFT, border=5)
        polsizer.Add(polsizer_2, flag=wx.LEFT, border=5)
        toplevsizer.Add(polsizer)

        pnl.SetSizer(toplevsizer)
        
        self.pol_x_txt = pol_x_txt
        self.pol_y_txt = pol_y_txt
        self.pol_z_txt = pol_z_txt
       
        ##buttons
        butsizer = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        butsizer.Add(okButton)
        butsizer.Add(closeButton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, 
            flag=wx.ALL|wx.EXPAND, border=5)
        vbox.Add(butsizer, 
            flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        
        okButton.Bind(wx.EVT_BUTTON, self.OnOk)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCancel)
        
    def OnClose(self, e):
        self.Destroy()
        
    def OnCancel(self, e):
        self.Destroy()

    def OnOk(self, e):
        dipole_params = {}

        name = self.name.GetLineText(0)
        name = nameGeneraor.GenerateNewName(name)
        dipole_params['name'] = name

        nameGeneraor.AddName(name)        

        if self.is_electric.GetValue()==True:
            dipole_params['type'] = "electric delta"
        elif self.is_magnetic.GetValue()==True:
            dipole_params['type'] = "magnetic delta"
        x = self.loc_x_txt.GetLineText(0)
        y = self.loc_y_txt.GetLineText(0)
        z = self.loc_z_txt.GetLineText(0)
        dipole_params['x'] = x
        dipole_params['y'] = y
        dipole_params['z'] = z

        Jx = self.pol_x_txt.GetLineText(0)
        Jy = self.pol_y_txt.GetLineText(0)
        Jz = self.pol_z_txt.GetLineText(0)
        dipole_params['Jx'] = Jx
        dipole_params['Jy'] = Jy
        dipole_params['Jz'] = Jz
        
        print("dipole_params: ", dipole_params)
        self.dipole_params = dipole_params
        self.Destroy()





class MediumParamsDialog(wx.Dialog):
    
    def __init__(self, *args, **kw):
        super(MediumParamsDialog, self).__init__(*args, **kw) 
            
        self.InitUI()
        self.SetSize((400, 350))
        self.SetTitle("Set Medium Parameters")
        
        self.medium_params = None
        
        
    def InitUI(self):

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        toplevbox = wx.StaticBox(pnl, label='Medium Parameters')
        toplevsizer = wx.StaticBoxSizer(toplevbox, orient=wx.VERTICAL)     
        
        
        ##---- name
        name_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        name_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'name:  '))
        name =  wx.TextCtrl(pnl, wx.NewId(), value=nameGeneraor.GenerateNewName('material'), size=wx.Size(200, 30))
        name_sizer.Add(name)
               
        toplevsizer.Add(name_sizer)
        
        self.name = name


        ##---- epsilon_r
        toplevsizer.Add(wx.StaticText(pnl, wx.NewId(), label='Permittivity:      '))

        
        eps_r_xx_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        eps_r_xx_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'epsilon_r^xx:  '))
        eps_r_xx =  wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(300, 30))
        eps_r_xx_sizer.Add(eps_r_xx)
        
        eps_r_yy_sizer = wx.BoxSizer(wx.HORIZONTAL)     
        eps_r_yy_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'epsilon_r^yy:  '))
        eps_r_yy =  wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(300, 30))
        eps_r_yy_sizer.Add(eps_r_yy)

        eps_r_zz_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        eps_r_zz_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'epsilon_r^zz:  '))
        eps_r_zz =  wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(300, 30))
        eps_r_zz_sizer.Add(eps_r_zz)
        
        toplevsizer.Add(eps_r_xx_sizer)
        toplevsizer.Add(eps_r_yy_sizer)
        toplevsizer.Add(eps_r_zz_sizer)
        
        self.eps_r_xx = eps_r_xx
        self.eps_r_yy = eps_r_yy
        self.eps_r_zz = eps_r_zz
        
        
        ##----- mu_r 
        toplevsizer.Add(wx.StaticText(pnl, wx.NewId(), label='Permeability:      '))
        mu_r_xx_sizer = wx.BoxSizer(wx.HORIZONTAL)      
        mu_r_xx_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'mu_r^xx:  '))
        mu_r_xx =  wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(300, 30))
        mu_r_xx_sizer.Add(mu_r_xx)

        mu_r_yy_sizer = wx.BoxSizer(wx.HORIZONTAL)      
        mu_r_yy_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'mu_r^yy:  '))
        mu_r_yy =  wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(300, 30))
        mu_r_yy_sizer.Add(mu_r_yy)

        mu_r_zz_sizer = wx.BoxSizer(wx.HORIZONTAL)     
        mu_r_zz_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'mu_r^zz:  '))
        mu_r_zz = wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(300, 30))
        mu_r_zz_sizer.Add(mu_r_zz)
        
        toplevsizer.Add(mu_r_xx_sizer)
        toplevsizer.Add(mu_r_yy_sizer)
        toplevsizer.Add(mu_r_zz_sizer)
        
        self.mu_r_xx = mu_r_xx
        self.mu_r_yy = mu_r_yy
        self.mu_r_zz = mu_r_zz
        

        pnl.SetSizer(toplevsizer)
        
        ##buttons
        butsizer = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        butsizer.Add(okButton)
        butsizer.Add(closeButton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, 
            flag=wx.ALL|wx.EXPAND, border=5)
        vbox.Add(butsizer, 
            flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        
        okButton.Bind(wx.EVT_BUTTON, self.OnOk)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCancel)
        
    def OnClose(self, e):
        self.Destroy()
        
    def OnCancel(self, e):
        self.Destroy()

    def OnOk(self, e):
        medium_params = {}
                
        medium_params['type'] = "diagonal"

        name = self.name.GetLineText(0)
        name = nameGeneraor.GenerateNewName(name)
        medium_params['name'] = name

        nameGeneraor.AddName(name)        
                
        medium_params['eps_r_xx'] = self.eps_r_xx.GetLineText(0)
        medium_params['eps_r_yy'] = self.eps_r_yy.GetLineText(0)
        medium_params['eps_r_zz'] = self.eps_r_zz.GetLineText(0)
        
        medium_params['mu_r_xx'] = self.mu_r_xx.GetLineText(0)
        medium_params['mu_r_yy'] = self.mu_r_yy.GetLineText(0)
        medium_params['mu_r_zz'] = self.mu_r_zz.GetLineText(0)
        
        print("medium_params: ", medium_params)
        self.medium_params = medium_params
        self.Destroy()




class BoundaryParamsDialog(wx.Dialog):
    
    def __init__(self, *args, **kw):
        super(BoundaryParamsDialog, self).__init__(*args, **kw) 
            
        self.InitUI()
        self.SetSize((400, 300))
        self.SetTitle("Set Boundary Parameters")
        
        self.boundary_params = None
        
        
    def InitUI(self):

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        toplevbox = wx.StaticBox(pnl, label='Boundary Parameters')
        toplevsizer = wx.StaticBoxSizer(toplevbox, orient=wx.VERTICAL)     
        is_pec = wx.RadioButton(pnl, label='Perfect Electric Conductor (PEC)', style=wx.RB_GROUP)   
        is_pmc = wx.RadioButton(pnl, label='Perfect Magnetic Condctor (PMC) ')   
        is_pml = wx.RadioButton(pnl, label='Perfect Matched Layer (PML) ')   
        toplevsizer.Add(is_pec)
        toplevsizer.Add(is_pmc)
        toplevsizer.Add(is_pml)
        
        self.is_pec = is_pec
        self.is_pmc = is_pmc
        self.is_pml = is_pml
        
        ##pml parameters
        pmlsizer = wx.BoxSizer(wx.VERTICAL)        
        pmlsizer.Add(wx.StaticText(pnl, wx.NewId(), label='PML Parameters: '))
        
        pml_thickness_sizer = wx.BoxSizer(wx.HORIZONTAL)        
        pml_thickness_sizer.Add(wx.StaticText(pnl, wx.NewId(), label='Thickness:         '))
        pml_thickness = wx.TextCtrl(pnl, wx.NewId(), value="0.1")
        pml_thickness_sizer.Add(pml_thickness, flag=wx.LEFT, border=5)

        pml_stretch_fact_sizer = wx.BoxSizer(wx.HORIZONTAL)        
        pml_stretch_fact_sizer.Add(wx.StaticText(pnl, wx.NewId(), label='Stretch factor: '))
        pml_stretch_fact = wx.TextCtrl(pnl, wx.NewId(), value="1.0-1.0j")
        pml_stretch_fact_sizer.Add(pml_stretch_fact, flag=wx.LEFT, border=5)

        pmlsizer.Add(pml_thickness_sizer, flag=wx.LEFT, border=5)
        pmlsizer.Add(pml_stretch_fact_sizer, flag=wx.LEFT, border=5)
        
        toplevsizer.Add(pmlsizer)
        
        self.pml_thickness = pml_thickness
        self.pml_stretch_fact = pml_stretch_fact
        

        pnl.SetSizer(toplevsizer)
        
        ##buttons
        butsizer = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        butsizer.Add(okButton)
        butsizer.Add(closeButton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, 
            flag=wx.ALL|wx.EXPAND, border=5)
        vbox.Add(butsizer, 
            flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        
        okButton.Bind(wx.EVT_BUTTON, self.OnOk)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCancel)
        
    def OnClose(self, e):
        self.Destroy()
        
    def OnCancel(self, e):
        self.Destroy()

    def OnOk(self, e):
        boundary_params = {}
        if self.is_pec.GetValue()==True:
            boundary_params['type'] = "PEC"
        elif self.is_pmc.GetValue()==True:
            boundary_params['type'] = "PMC"
        elif self.is_pml.GetValue()==True:
            boundary_params['type'] = "PML"

            pml_thickness = self.pml_thickness.GetLineText(0)
            boundary_params['pml_thickness'] = pml_thickness
            pml_stretch_fact = self.pml_stretch_fact.GetLineText(0)
            boundary_params['pml_stretch_fact'] = pml_stretch_fact
        
        print("boundary_params: ", boundary_params)
        self.boundary_params = boundary_params
        self.Destroy()




class OutputParamsDialog(wx.Dialog):
    
    def __init__(self, *args, **kw):
        super(OutputParamsDialog, self).__init__(*args, **kw) 
            
        self.InitUI()
        self.SetSize((400, 300))
        self.SetTitle("Set Boundary Parameters")
        
        self.params = None
        
        
    def InitUI(self):

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        toplevbox = wx.StaticBox(pnl, label='Output')
        toplevsizer = wx.StaticBoxSizer(toplevbox, orient=wx.VERTICAL)     

        ##---- name
        name_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        name_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'name:  '))
        name =  wx.TextCtrl(pnl, wx.NewId(), value=nameGeneraor.GenerateNewName("field"), size=wx.Size(200, 30))
        name_sizer.Add(name)
               
        toplevsizer.Add(name_sizer)
        
        self.name = name


        ##---- type
        is_entire = wx.RadioButton(pnl, label='Entire Domain', style=wx.RB_GROUP)   
        is_cut = wx.RadioButton(pnl, label='Cut Plane')   
        
        r_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        r_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'r:  '))
        x =  wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(100, 30))
        y =  wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(100, 30))
        z =  wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(100, 30))
        #x.Enable(False)
        r_sizer.Add(x)
        r_sizer.Add(y)
        r_sizer.Add(z)
        
        cut_plane_sizer = wx.BoxSizer(wx.HORIZONTAL)        
        is_x = wx.RadioButton(pnl, label='x', style=wx.RB_GROUP)   
        is_y = wx.RadioButton(pnl, label='y')   
        is_z = wx.RadioButton(pnl, label='z')   
        cut_plane_sizer.Add(is_x)
        cut_plane_sizer.Add(is_y)
        cut_plane_sizer.Add(is_z)

        toplevsizer.Add(is_entire)
        toplevsizer.Add(is_cut)
        toplevsizer.Add(r_sizer)
        toplevsizer.Add(cut_plane_sizer)
        
        self.is_entire = is_entire
        self.is_cut = is_cut
        self.x = x
        self.y = y
        self.z = z
        self.is_x = is_x
        self.is_y = is_y
        self.is_z = is_z
        
        ##--- field
        fieldsizer = wx.BoxSizer(wx.VERTICAL)        
        fieldsizer.Add(wx.StaticText(pnl, wx.NewId(), label='Field: '))
        
        field_exyz_sizer_sizer = wx.BoxSizer(wx.HORIZONTAL)        
        is_Ex = wx.RadioButton(pnl, label='E_x', style=wx.RB_GROUP)   
        is_Ey = wx.RadioButton(pnl, label='E_y')   
        is_Ez = wx.RadioButton(pnl, label='E_z')   
        field_exyz_sizer_sizer.Add(is_Ex)
        field_exyz_sizer_sizer.Add(is_Ey)
        field_exyz_sizer_sizer.Add(is_Ez)

        fieldsizer.Add(field_exyz_sizer_sizer, flag=wx.LEFT, border=5)

        toplevsizer.Add(fieldsizer)
        
        self.is_Ex = is_Ex
        self.is_Ey = is_Ey
        self.is_Ez = is_Ez

        pnl.SetSizer(toplevsizer)
        
        ##buttons
        butsizer = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        butsizer.Add(okButton)
        butsizer.Add(closeButton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, 
            flag=wx.ALL|wx.EXPAND, border=5)
        vbox.Add(butsizer, 
            flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        
        okButton.Bind(wx.EVT_BUTTON, self.OnOk)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCancel)
        
    def OnClose(self, e):
        self.Destroy()
        
    def OnCancel(self, e):
        self.Destroy()

    def OnOk(self, e):
        params = {}
        
        name = self.name.GetLineText(0)
        name = nameGeneraor.GenerateNewName(name)
        params['name'] = name
        
        nameGeneraor.AddName(name)
        
        if self.is_entire.GetValue()==True:
            params['type'] = "entire"
        elif self.is_cut.GetValue()==True:
            params['type'] = "cut"
            
        x = self.x.GetLineText(0)
        y = self.y.GetLineText(0)
        z = self.z.GetLineText(0)
        
        params['r'] = [x, y, z]

        if self.is_x.GetValue()==True:
            params['plane'] = "x"
        elif self.is_y.GetValue()==True:
            params['plane'] = "y"
        elif self.is_z.GetValue()==True:
            params['plane'] = "z"
        
        if self.is_Ex.GetValue()==True:
            params['field'] = "Ex"
        elif self.is_Ey.GetValue()==True:
            params['field'] = "Ey"
        elif self.is_Ez.GetValue()==True:
            params['field'] = "Ez"
        
        print("output_params: ", params)
        self.params = params
        self.Destroy()



class SimulationParamsDialog(wx.Dialog):
    
    def __init__(self, *args, **kw):
        super(SimulationParamsDialog, self).__init__(*args, **kw) 
            
        self.InitUI()
        self.SetSize((400, 300))
        self.SetTitle("Set Simulation Parameters")
        
        self.simulation_params = None
        
        
    def InitUI(self):

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        toplevbox = wx.StaticBox(pnl, label='Simulation Parameters')
        toplevsizer = wx.StaticBoxSizer(toplevbox, orient=wx.VERTICAL)     
        
        
        ##---- d_r
        toplevsizer.Add(wx.StaticText(pnl, wx.NewId(), label='Grid      '))
        
        dr_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        dr_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'dr:  '))
        dx =  wx.TextCtrl(pnl, wx.NewId(), value="0.1", size=wx.Size(100, 30))
        dy =  wx.TextCtrl(pnl, wx.NewId(), value="0.1", size=wx.Size(100, 30))
        dz =  wx.TextCtrl(pnl, wx.NewId(), value="0.1", size=wx.Size(100, 30))
        dr_sizer.Add(dx)
        dr_sizer.Add(dy)
        dr_sizer.Add(dz)
                
        toplevsizer.Add(dr_sizer)

        self.dx = dx
        self.dy = dy
        self.dz = dz

        ##---- Stability factor
        toplevsizer.Add(wx.StaticText(pnl, wx.NewId(), label='Stability factor      '))

        S_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        S_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'S:  '))
        S =  wx.TextCtrl(pnl, wx.NewId(), value="0.99", size=wx.Size(100, 30))
        S_sizer.Add(S)
                
        toplevsizer.Add(S_sizer)

        self.S = S
                
        
        ##---- Save outputs every
        toplevsizer.Add(wx.StaticText(pnl, wx.NewId(), label='Save outputs at N*dt samples: '))

        save_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        save_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'N:  '))
        N_save =  wx.TextCtrl(pnl, wx.NewId(), value="1", size=wx.Size(100, 30))
        save_sizer.Add(N_save)
                
        toplevsizer.Add(save_sizer)

        self.N_save = N_save
        
        ##--
                
        pnl.SetSizer(toplevsizer)
        
        ##buttons
        butsizer = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        butsizer.Add(okButton)
        butsizer.Add(closeButton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, 
            flag=wx.ALL|wx.EXPAND, border=5)
        vbox.Add(butsizer, 
            flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        
        okButton.Bind(wx.EVT_BUTTON, self.OnOk)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCancel)
        
    def OnClose(self, e):
        self.Destroy()
        
    def OnCancel(self, e):
        self.Destroy()

    def OnOk(self, e):
        simulation_params = {}

        simulation_params["dx"] = self.dx.GetLineText(0)
        simulation_params["dy"] = self.dy.GetLineText(0)
        simulation_params["dz"] = self.dz.GetLineText(0)
        
        simulation_params["S"] = self.S.GetLineText(0)
        
        simulation_params["N_save"] = self.N_save.GetLineText(0)

        print("simulation_params: ", simulation_params)
        self.simulation_params = simulation_params
        self.Destroy()


class GlobalVarsDialog(wx.Dialog):
    
    def __init__(self, glob_vars, *args, **kw):
        super(GlobalVarsDialog, self).__init__(*args, **kw) 
            
        self.InitUI(glob_vars)
        self.SetSize((400, 300))
        self.SetTitle("Set Global Variables")
        
        self.params = None
        
        
    def InitUI(self, glob_vars):

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        toplevbox = wx.StaticBox(pnl, label='Global Variables')
        toplevsizer = wx.StaticBoxSizer(toplevbox, orient=wx.VERTICAL)     
        
        
        ##---- d_r
        toplevsizer.AddSpacer(20)
        toplevsizer.Add(wx.StaticText(pnl, wx.NewId(), label='New Varaible    '))
        
        name_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        name_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'name:  '))
        name =  wx.TextCtrl(pnl, wx.NewId(), value="alpha", size=wx.Size(100, 30))
        name_sizer.Add(name)
        ## TODO: some variable names such as 'var' are problematic
        
        val_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        val_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'value:  '))
        val =  wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(100, 30))
        val_sizer.Add(val)
                
        toplevsizer.AddSpacer(10)
        toplevsizer.Add(name_sizer)
        toplevsizer.Add(val_sizer)

        self.name = name
        self.val  = val
        
        ##---- table
        
        n_row = min(len(glob_vars), 5)
        
        grid = wx.grid.Grid(pnl, size=(300, 100))
        grid.CreateGrid(n_row, 2)
        grid.EnableEditing(False)
        
        grid.SetRowLabelSize(40)
        grid.SetColSize(0, 130)
        grid.SetColSize(1, 130)
        
        _ind_ = 0
        for g_var in glob_vars:
            grid.SetCellValue(_ind_, 0, g_var)
            grid.SetCellValue(_ind_, 1, str(glob_vars[g_var]))
            _ind_ += 1
        
        toplevsizer.AddSpacer(20)
        toplevsizer.Add(grid, border=10)
                
        pnl.SetSizer(toplevsizer)
        
        ##buttons
        butsizer = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        butsizer.Add(okButton)
        butsizer.Add(closeButton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, 
            flag=wx.ALL|wx.EXPAND, border=5)
        vbox.Add(butsizer, 
            flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        
        okButton.Bind(wx.EVT_BUTTON, self.OnOk)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCancel)
        
    def OnClose(self, e):
        self.Destroy()
        
    def OnCancel(self, e):
        self.Destroy()

    def OnOk(self, e):
        params = {}

        params["name"] = self.name.GetLineText(0)
        params["value"] = float(self.val.GetLineText(0))

        print("variable params: ", params)
        self.params = params
        self.Destroy()



class AssignMediumDialog(wx.Dialog):
    
    def __init__(self, medium_list_str, *args, **kw):
        super(AssignMediumDialog, self).__init__(*args, **kw) 
            
        self.InitUI(medium_list_str)
        self.SetSize((400, 300))
        self.SetTitle("Assign Medium")
        
        self.params = None
        
        
    def InitUI(self, medium_list_str):

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        toplevbox = wx.StaticBox(pnl, label='Select Medium')
        toplevsizer = wx.StaticBoxSizer(toplevbox, orient=wx.VERTICAL)     
        
        ##---- name
        mediumlist_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        mediumlist_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'medium:  '))
        medium_list =  wx.ComboBox(pnl, choices=medium_list_str)
        mediumlist_sizer.Add(medium_list)
               
        toplevsizer.Add(mediumlist_sizer)
        
        self.medium_list = medium_list

                
        pnl.SetSizer(toplevsizer)
        
        ##buttons
        butsizer = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        butsizer.Add(okButton)
        butsizer.Add(closeButton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, 
            flag=wx.ALL|wx.EXPAND, border=5)
        vbox.Add(butsizer, 
            flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        
        okButton.Bind(wx.EVT_BUTTON, self.OnOk)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCancel)
        
    def OnClose(self, e):
        self.Destroy()
        
    def OnCancel(self, e):
        self.Destroy()

    def OnOk(self, e):
        params = {}
        
        sel_ind = self.medium_list.GetSelection()
        if sel_ind != wx.NOT_FOUND:
            params['name'] = self.medium_list.GetString(sel_ind)
            
            print("assign medium params: ", params)
        self.params = params
        self.Destroy()



class CubeParamsDialog(wx.Dialog):
    
    def __init__(self, *args, **kw):
        super(CubeParamsDialog, self).__init__(*args, **kw) 
            
        self.InitUI()
        self.SetSize((400, 300))
        self.SetTitle("Cube Parameters")
        
        self.cube_params = None
        
        
    def InitUI(self):

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        toplevbox = wx.StaticBox(pnl, label='Cube Parameters')
        toplevsizer = wx.StaticBoxSizer(toplevbox, orient=wx.VERTICAL)     
        
        ##---- name
        name_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        name_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'name:  '))
        name =  wx.TextCtrl(pnl, wx.NewId(), value=nameGeneraor.GenerateNewName("cube"), size=wx.Size(200, 30))
        name_sizer.Add(name)
        color = wx.ColourPickerCtrl(pnl)
        random_color = np.random.randint(0, 256, size=4)
        random_color[3] = 255
        color.SetColour(random_color)
        name_sizer.Add(color)
               
        toplevsizer.Add(name_sizer)
        
        self.name = name
        self.color = color

        ##---- r0
        toplevsizer.Add(wx.StaticText(pnl, wx.NewId(), label='Lower left corner '))
        
        r0_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        r0_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'r0:  '))
        r0_x =  wx.TextCtrl(pnl, wx.NewId(), value="0.0", size=wx.Size(100, 30))
        r0_y =  wx.TextCtrl(pnl, wx.NewId(), value="0.0", size=wx.Size(100, 30))
        r0_z =  wx.TextCtrl(pnl, wx.NewId(), value="0.0", size=wx.Size(100, 30))
        r0_sizer.Add(r0_x)
        r0_sizer.Add(r0_y)
        r0_sizer.Add(r0_z)
                
        toplevsizer.Add(r0_sizer)

        self.r0_x = r0_x
        self.r0_y = r0_y
        self.r0_z = r0_z

        ##---- Stability factor
        toplevsizer.Add(wx.StaticText(pnl, wx.NewId(), label='Upper right corner '))
        
        r1_sizer = wx.BoxSizer(wx.HORIZONTAL)    
        r1_sizer.Add(wx.StaticText(pnl, wx.NewId(), label=r'r1:  '))
        r1_x =  wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(100, 30))
        r1_y =  wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(100, 30))
        r1_z =  wx.TextCtrl(pnl, wx.NewId(), value="1.0", size=wx.Size(100, 30))
        r1_sizer.Add(r1_x)
        r1_sizer.Add(r1_y)
        r1_sizer.Add(r1_z)
                
        toplevsizer.Add(r1_sizer)

        self.r1_x = r1_x
        self.r1_y = r1_y
        self.r1_z = r1_z
                
        pnl.SetSizer(toplevsizer)
        
        ##buttons
        butsizer = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        butsizer.Add(okButton)
        butsizer.Add(closeButton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, 
            flag=wx.ALL|wx.EXPAND, border=5)
        vbox.Add(butsizer, 
            flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        
        okButton.Bind(wx.EVT_BUTTON, self.OnOk)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCancel)
        
    def OnClose(self, e):
        self.Destroy()
        
    def OnCancel(self, e):
        self.Destroy()

    def OnOk(self, e):
        cube_params = {}

        cube_params['type'] = 'cube'
        
        name = self.name.GetLineText(0)
        name = nameGeneraor.GenerateNewName(name)
        cube_params['name'] = name
        
        nameGeneraor.AddName(name)
        
        color = self.color.GetColour().Get()
        cube_params['color'] = tuple(np.array(color)/255.0)

        r0_x = self.r0_x.GetLineText(0)
        r0_y = self.r0_y.GetLineText(0)
        r0_z = self.r0_z.GetLineText(0)
        
        cube_params["r0"] = [r0_x, r0_y, r0_z]

        r1_x = self.r1_x.GetLineText(0)
        r1_y = self.r1_y.GetLineText(0)
        r1_z = self.r1_z.GetLineText(0)
        
        cube_params["r1"] = [r1_x, r1_y, r1_z]

        print("cube_params: ", cube_params)
        self.cube_params = cube_params
        self.Destroy()



##--------------------  popup menue

class PopupMenu_Object(wx.Menu):
    def __init__(self, obj_type=None):
        wx.Menu.__init__(self)
        
        self.request = None
        
        class reuqest_types(Enum): 
            volobj_assign_material = 1
            volobj_show_hide = 2
            #-
            none_set_units = 20
            #-
            cube_set_as_simulation_box = 40
        
        self.request_types = reuqest_types

        if obj_type=="VOL_OBJ" or obj_type=="CUBE_NOROT":
            ## "CUBE_NOROT" : non rotated cube
            item = wx.MenuItem(self, wx.NewId(), r"Show/Hide")
            self.Append(item)
            self.Bind(wx.EVT_MENU, self.OnShowHide, item)
            item = wx.MenuItem(self, wx.NewId(), "Assign Material")
            self.Append(item)
            self.Bind(wx.EVT_MENU, self.OnAssignMaterial, item)
            if obj_type=="CUBE_NOROT":
                item = wx.MenuItem(self, wx.NewId(), "Set as Simulation Box")
                self.Append(item)
                self.Bind(wx.EVT_MENU, self.OnSetSimulationBox, item)
                
        elif obj_type==None:
            item = wx.MenuItem(self, wx.NewId(), "Set Units")
            self.Append(item)
            self.Bind(wx.EVT_MENU, self.OnSetUnits, item)

    def OnSetUnits(self, event):
        self.request = self.request_types.none_set_units
        print("Set Units")
        event.Skip()

    def OnAssignMaterial(self, event):
        self.request = self.request_types.volobj_assign_material
        print("Assign material")
        event.Skip()
        
    def OnShowHide(self, event):
        self.request = self.request_types.volobj_show_hide
        print(r"Show/Hide")
        event.Skip()

    def OnSetSimulationBox(self, event):
        self.request = self.request_types.cube_set_as_simulation_box
        print(r"Set as simulation box.")
        event.Skip()



class PopupMenu_Source(wx.Menu):
    def __init__(self, obj_type=None):
        wx.Menu.__init__(self)
        
        self.request = None
        
        class reuqest_types(Enum): 
            dipole_change_parameters = 1
        
        self.request_types = reuqest_types

        if obj_type=="dipole":
            item = wx.MenuItem(self, wx.NewId(), "Change Parameters")
            self.Append(item)
            item.Enable(False)
            self.Bind(wx.EVT_MENU, self.OnChangeParameters, item)

    def OnChangeParameters(self, event):
        self.request = self.request_types.dipole_change_parameters
        print("Change Parameters")
        event.Skip()



##--------------------  Exception handling
import traceback
import wx.lib.agw.genericmessagedialog as GMD

class ExceptionDialog(GMD.GenericMessageDialog):
    def __init__(self, msg):
        GMD.GenericMessageDialog.__init__(self, None, msg, "Exception!", wx.OK|wx.ICON_ERROR)        
 
 

def MyExceptionHookDlg(etype, value, trace):
    """
    Handler for all unhandled exceptions.

    :param `etype`: the exception type (`SyntaxError`, `ZeroDivisionError`, etc...);
    :type `etype`: `Exception`
    :param string `value`: the exception error message;
    :param string `trace`: the traceback header, if any (otherwise, it prints the
     standard Python header: ``Traceback (most recent call last)``.
    """
    frame = wx.GetApp().GetTopWindow()
    tmp = traceback.format_exception(etype, value, trace)
    exception = "".join(tmp)

    dlg = ExceptionDialog(exception)
    dlg.ShowModal()
    dlg.Destroy()    
    

def MyExceptionHook(etype, value, trace):
    tmp = traceback.format_exception(etype, value, trace)
    exception = "".join(tmp)

    print(exception)
    sys.stdout.flush()
    
    
    
    




