# chat client



__all__ = ["GUIChatClient"]




#from multiprocessing import Process, Queue
from threading import Thread
import queue
from tkinter import *
import time
from random import randint

import socket

##TODO:
""" - assure proper closing of sockets
    - send send socket connectionOk confirmation to gui
    - handle send-receive sockets exceptions
"""

class GUIChatClient(Thread):
    def __init__(self, que=None):
        self.queue = que
        self.que_send = None
        self.que_send_callBack = None
        self.que_recv = None
        self.quit = False
        self.connStablished = False
        self.connRecvOk = False
        self.connSendOk = False
        Thread.__init__(self)
        self.start()
        return
        
    def SetCommChannel(self, que):
        self.queue = que
        return
                
    def run(self):
        self.root = Tk()

        label = Label(self.root, text="Chat Client")#, font=LARGE_FONT)
        label.pack(pady=10,padx=10)
        #self.label = label

        self.Tb0 = Text(self.root, height=30, width=150)
        self.Tb0.pack()
        self.Tb0.config(state=DISABLED)
        
        self.Tb1 = Text(self.root, height=3, width=150)
        self.Tb1.pack()
                
        self.But = Button(self.root, text ="Send", command=self.OnSendBut)
        self.But.pack()
        
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
        del self.Tb0
        del self.Tb1
        del self.But
        del self.root
        return

    def CommandLoop(self):
        if not self.quit:
            if self.queue != None:
                if not self.queue.empty():                
                    comm = self.queue.get()
                    self.ProcessCommand(comm)
                    self.queue.task_done()
            if self.que_recv != None:
                if not self.que_recv.empty():                
                    comm = self.que_recv.get()
                    if comm[0]=='msg':
                        msg_rec = comm[1]
                        self.PrintMsg('OO: '+msg_rec)
                    elif comm[0]=='err':
                        err = comm[1]
                        self.PrintMsg('Error: {} \n'.format(err))
                    elif comm[0]=='conn recv':
                        if comm[1]==True:
                            self.connRecvOk = True
                            if self.connSendOk==True:
                                self.connStablished = True
                                self.PrintMsg('Connection stablished\n')
                    self.que_recv.task_done()
            if self.que_send_callBack != None:
                if not self.que_send_callBack.empty():                
                    comm = self.que_send_callBack.get()
                    if comm[0]=='err':
                        err = comm[1]
                        self.PrintMsg('Error: {} \n'.format(err))
                    elif comm[0]=='conn send':
                        if comm[1]==True:
                            self.connSendOk = True
                            if self.connRecvOk==True:
                                self.connStablished = True
                                self.PrintMsg('Connection stablished\n')
                    self.que_send_callBack.task_done()
            self.root.after(1000, self.CommandLoop)
        return

    def ProcessCommand(self, comm):
        return
 
    def OnSendBut(self):
        msg = self.Tb1.get(1.0, END)
        self.ClearInput()        
        if msg.startswith('##'):
            msg = msg[2:]
            msg = msg.strip(' ').strip('\n').strip('\r')
            msg_arr = msg.split(' ')
            err = None
            if msg_arr[0]=='listen':
                if len(msg_arr)==4:
                    ip = msg_arr[1]
                    port_send = int(msg_arr[2])
                    port_recv = int(msg_arr[3])
                    
                    try:
                        #if ip=='localhost':
                        #    ip = self.GetLocalIP()
                        self.PrintMsg('Listening: ip:{}, send port:{}, receive port:{} ... \n'.format(ip, port_send, port_recv))

                        self.que_send = queue.Queue()
                        self.que_send_callBack = queue.Queue()
                        self.soc_send = ServerSocketSend(ip, port_send, self.que_send, self.que_send_callBack)

                        self.que_recv = queue.Queue()
                        self.soc_recv = ServerSocketRecv(ip, port_recv, self.que_recv)

                    except Exception as e:
                        self.PrintMsg(str(e))
                    except:
                        self.PrintMsg(sys.exc_info()[0])

                else:
                    err = 'bad arguments'
            elif msg_arr[0]=='connect':
                if len(msg_arr)==4:
                    ip = msg_arr[1]
                    port_send = int(msg_arr[2])
                    port_recv = int(msg_arr[3])
                    
                    try:
                        self.PrintMsg('Connecting: ip:{}, send port:{}, receive port:{} ... \n'.format(ip, port_send, port_recv))

                        self.que_send = queue.Queue()
                        self.que_send_callBack = queue.Queue()
                        self.soc_send = ClientSocketSend(ip, port_send, self.que_send, self.que_send_callBack)

                        self.que_recv = queue.Queue()
                        self.soc_recv = ClientSocketRecv(ip, port_recv, self.que_recv)

                    except RuntimeError as e:
                        self.PrintMsg(str(e))
                        self.PrintMsg(sys.exc_info()[0])
                    except Exception as e:
                        self.PrintMsg(str(e))
                    except:
                        self.PrintMsg(sys.exc_info()[0])

                else:
                    err = 'bad arguments'
            elif msg_arr[0]=='getip':
                ip = self.GetLocalIP()
                self.PrintMsg('ip = {} \n'.format(ip))
            if err!=None:
                for msg_i in msg_arr:
                    self.PrintMsg(msg_i+' - ')
        elif self.connStablished:
            self.PrintMsg('Me: '+msg)
            try:
                self.que_send.put(['msg', msg])
            except RuntimeError as e:
                self.PrintMsg(str(e))
            except:
                self.PrintMsg(sys.exc_info()[0])
        else:
            self.PrintMsg('Not connected :::: ' + msg)
    
        
    def PrintMsg(self, msg):
        self.Tb0.config(state=NORMAL)
        self.Tb0.insert(END, msg)
        self.Tb0.config(state=DISABLED)
    
    def ClearInput(self):
        self.Tb1.delete(1.0, END)
        
    def GetLocalIP(self):
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(('8.8.8.8', 0))  # connecting to a UDP address doesn't send packets
        local_ip_address = s.getsockname()[0]

        return local_ip_address
    

TIME_SLEEP = 1.0
class ServerSocketRecv(Thread):
    def __init__(self, ip, port, que):
        self.queue = que
        self.buff_size = 1024
        self.stop = False
        self.ip = ip
        self.port = port
        Thread.__init__(self)
        self.daemon = True
        self.start()
        return
                        
    def run(self):
        self.StartSocket()
        while not self.stop:
            try:
                msg_rec = self.conn.recv(self.buff_size)
                #print('msg_rec', msg_rec)
                self.queue.put(['msg', msg_rec.decode('utf-8')])
            except Exception as e:
                self.queue.put(['err', str(e)])
            time.sleep(TIME_SLEEP)

    def StartSocket(self):
        try:
            self.soc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.soc.bind((self.ip.encode('utf-8'), self.port))
            self.soc.listen(1)
            self.conn, addr = self.soc.accept()  
            #print('ServerSocketRecv connected: ', addr)
            self.queue.put(['conn recv', True])
        except Exception as e:
            self.queue.put(['err', str(e)])


class ServerSocketSend(Thread):
    def __init__(self, ip, port, que, que_callback):
        self.queue = que
        self.queue_callback = que_callback
        self.buff_size = 1024
        self.stop = False
        self.ip = ip
        self.port = port
        Thread.__init__(self)
        self.daemon = True
        self.start()
        return
                        
    def run(self):
        self.StartSocket()
        while not self.stop:
            if not self.queue.empty():                
                comm = self.queue.get()
                if comm[0]=='msg':
                    try:
                        msg_send = comm[1]
                        #print('msg_send: ', msg_send)
                        self.conn.send(msg_send.encode('utf-8'))
                    except Exception as e:
                        self.queue_callback.put(['err', str(e)])
                self.queue.task_done()
            time.sleep(TIME_SLEEP)

    def StartSocket(self):
        try:
            self.soc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.soc.bind((self.ip.encode('utf-8'), self.port))
            self.soc.listen(1)
            self.conn, addr = self.soc.accept()  
            #print('ServerSocketSend connected', addr)
            self.queue_callback.put(['conn send', True])
        except Exception as e:
            self.queue_callback.put(['err', str(e)])
        

class ClientSocketRecv(Thread):
    def __init__(self, ip, port, que):
        self.queue = que
        self.buff_size = 1024
        self.stop = False
        self.ip = ip
        self.port = port
        Thread.__init__(self)
        self.daemon = True
        self.start()
        return
                        
    def run(self):
        self.StartSocket()
        while not self.stop:
            try:
                msg_rec = self.soc.recv(self.buff_size)
                #print('msg_rec', msg_rec)
                self.queue.put(['msg', msg_rec.decode('utf-8')])
            except Exception as e:
                self.queue.put(['err', str(e)])
            time.sleep(TIME_SLEEP)

    def StartSocket(self):
        try:
            self.soc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.soc.connect((self.ip.encode('utf-8'), self.port))
            #print('ClientSocketRecv connected')
            self.queue.put(['conn recv', True])
        except Exception as e:
            self.queue.put(['err', str(e)])


class ClientSocketSend(Thread):
    def __init__(self, ip, port, que, que_callback):
        self.queue = que
        self.queue_callback = que_callback
        self.buff_size = 1024
        self.stop = False
        self.ip = ip
        self.port = port
        Thread.__init__(self)
        self.daemon = True
        self.start()
        return
                        
    def run(self):
        self.StartSocket()
        while not self.stop:
            if not self.queue.empty():                
                comm = self.queue.get()
                if comm[0]=='msg':
                    msg_send = comm[1]
                    #print('msg_send: ', msg_send)
                    try:
                        self.soc.send(msg_send.encode('utf-8'))
                    except Exception as e:
                        self.queue_callback.put(['err', str(e)])
                self.queue.task_done()
            time.sleep(TIME_SLEEP)

    def StartSocket(self):
        try:
            self.soc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.soc.connect((self.ip.encode('utf-8'), self.port))
            #print('ClientSocketSend connected')
            self.queue_callback.put(['conn send', True])
        except Exception as e:
            self.queue_callback.put(['err', str(e)])
        


        

