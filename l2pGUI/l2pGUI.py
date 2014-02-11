#!/usr/bin/env python

import sys, os
import numpy as np
import Tkinter as Tk
import socket, multiprocessing
import subprocess
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
import time
import datetime as dt
import noaasun
#import funplot as fp
#import plot_plist2 as pp2
#import satpar as sp


def dataFakeRead(f, init_pos=None):
    """Reads lines from file from specified position to end.
    
    For testing purposes: let f be a saved l2planes output; calling
    this function at regular intervals will simulate real time data,
    provided the f is also being written continuously. Launch function 
    dataFakeGen from another Python terminal to generate such a file.
    
    This function is a generator which behaves similarly to Unix tail,
    but the caller is responsible for keeping track of file position
    and passing it as an argument.
    """
    if init_pos is None:
        f.seek(0, 2)
        init_pos = f.tell()
    f.seek(init_pos)
    lines = []
    i = 0
    while True:
        line = f.readline()
        if not line:
            print 'EOF ' + '%d lines read' % i
            pos = f.tell()
            yield lines, pos
        l = line.split()
        if len(l) == 12:
            lines.append(line)
        i += 1
        if i > 5000:
            print line
            break
    pos = f.tell()
    yield lines, pos
    

def dataFakeGen():
    """Reads data chunks from a saved l2planes output file and writes
    them to a second file at periodic intervals. Useful to simulate 
    incoming real time data from l2planes in conjunction with dataFakeRead().
    """
    import time
    f_write = open('fakeoutput.out', 'w')
    f_read = open('/home/jose/data/l2planes/l2p18vlong.out', 'r')
    N = 60
    while True:
        for n in range(N):
            line = f_read.readline()    
            f_write.write(line)
            f_write.flush()    
        sys.stdout.flush()
        time.sleep(0.5)

    
def colourMaplimits(value_limits=(0, 80), colourmap_limits=(0.05, 0.85)):
    """Compute a and b coefficients that bring values from value_limits
    to colourmap_limits in a linear way.
    
    I.e., solve y = (x - a) / b
    
    where x is a value from value_limits and y the corresponding one
    in colourmap_limits.
    
    This function is not called during L2pRadar execution,
    just used to pre-compute a, b to choose line colours.
    """
    vlow, vhigh = value_limits
    clow, chigh = colourmap_limits
    a = vhigh * clow / (clow - chigh)
    b = -a / clow
    return a, b
    
        
class Plane():
    '''Makes planes.'''
    def __init__(self, line, minel=10):
        l = line.split()
        self.minel = minel
        self.id = l[2]
        self.epc = []
        self.last_epoch = float(l[1])
        #self.lat = []
        #self.lon = []
        #self.alt = []
        #self.ran = []
        self.az = []
        self.el = []
        self.maxel = 0
        self.gaps = 0
        self.addLine(l)

    # 56395 40400.326   4ca626 RYR8JT   50.97158 -0.61729 29525 68.6683
    # 280.17873692 7.16197030   -474.0 191.0 1088   0.00 0.00

    def addLine(self, l):
        '''Adds data lines'''
        # l is an already splitted line (list of line contents)
        if len(l) == 13:
            if float(l[9]) < self.minel:
                return
            self.epc.append(float(l[1]))
            if self.epc[-1] < self.last_epoch:
                self.epc[-1] += 86000
            if self.epc[-1] - self.last_epoch > 600:
                self.gaps = 1
                del(self.epc[-1])
                return
            self.last_epoch = self.epc[-1]
            #self.lat.append(float(l[4]))
            #self.lon.append(float(l[5]))
            #self.alt.append(float(l[6]) * 0.3048)
            #self.ran.append(float(l[7]))
            self.az.append(np.pi / 180 * float(l[8]))
            self.el.append(90 - float(l[9]))
            self.maxel = self.el[-1] if self.el[-1] > self.maxel else self.maxel


def addPlanes(planeLines, planes_dict, minel=15):
    """Processes list of plane data lines and updates planes dictionary
    accordingly.
    
    Parameters
    ----------
    planeLines: list of plane lines from l2planes to process
    planes_dict: dictionary storing planes (keys: id; values: Plane instances)
    minel: elevation cutoff
    """
    P = planes_dict
    for line in planeLines:
        l = line.split()
        plane_id = l[2]
        if not P.has_key(plane_id):
            P[plane_id] = Plane(line, minel)
        else:
            P[plane_id].addLine(l)

    initial_lenP = len(P)
    if len(P) == 0 or len(planeLines) == 0:
        return P
    # Remove planes with no epochs (below elevation cutoff) and those
    # for which no beacons have been received for more than 15 seconds
    keys1 = {k for k,v in P.iteritems() if len(v.epc) == 0}
    last_epoch = float(l[1])
    keys2 = {k for k,v in P.iteritems() if abs(last_epoch - v.last_epoch) > 15}
    keys_to_remove = keys1.union(keys2)
    for key in list(keys_to_remove):
        del P[key]
    print 'len(P): {} {} (-{})'.format(initial_lenP, len(P), len(keys_to_remove))
    return P


class ConnectionL2P():
    """Handles connection with l2planes server, creates queue subproces 
    and sends data from the former to the latter.
    """
    def __init__(self, l2p_HOST, l2p_PORT):
        self.planeQueue = multiprocessing.Queue()
        try:
            self.connSocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        except socket.error, msg:
            print 'Failed to create socket: ' + str(msg[0]) + ' ; ' + msg[1]
            raise SystemExit
        print '\nSocket Created\n'
        self.connSocket.settimeout(10)
        try:        
            self.connSocket.connect((l2p_HOST, l2p_PORT))
        except socket.error, msg:
            print 'Failed to connect: ' + str(msg)
            print '\nIs the listen2planes server running?\n'
            raise SystemExit
        
        self.procWorker = multiprocessing.Process(target=self.receive_proc,
                          args=(self.connSocket, self.planeQueue))
        self.procWorker.start()
        
        self.fdump = open('dump.out', 'w')

    def receive_proc(self, connSocket, planeQueue):
        """Requests one data line from l2planes and sends it to the queue.
        """      
        while True:
            # Send string expected by listen2planes from clients
            try:
                connSocket.send('reader')
            except socket.error, msg:
                print('Failed to receive from server: '.format(msg))
                time.sleep(0.5)                
            try:
                data = connSocket.recv(1024)
            except socket.error, msg:
                print('Failed to receive from sserver: {}'.format(msg))
                time.sleep(0.5)                
            else:
                planeQueue.put(data)
        return
        
    def dump_queue(self):
        """Retrieves all the data lines from the queue."""
        plines = []
        tlines = []
        #print('Queue size: {}'.format(self.planeQueue.qsize())),
        while True:
            this_line = self.planeQueue.get()
            self.fdump.write(this_line + ' \n')
            print this_line
            L = len(this_line.split())
            if L == 13:
                plines.append(this_line)
            elif L == 6:
                tlines.append(this_line)
            if self.planeQueue.qsize() == 0:
                break
        return plines, tlines


class L2pRadar(Tk.Tk):
    """Real-time polar plot of ADS-B planes data received from listen2planes
    """
    def __init__(self, parent):
        Tk.Tk.__init__(self, parent)
        self.root = Tk.Tk._root(self)
        self.root.configure(background='black')
        self.l2p_HOST = '193.61.194.29'
        self.tmpath = os.path.expanduser('~/.plotsched_tmp')
        self.appRunning = False
        self.visHEO = False
        self.P = {}
        self.frameCtrls = Tk.Frame()
        self.frameCtrls.pack(side='left')
        self.buttonLimitUp = Tk.Button(self.frameCtrls, text='Up',
                             command=self.plotLimitUp, bg='grey')
        self.buttonLimitDown = Tk.Button(self.frameCtrls, text='Dn',
                               command=self.plotLimitDown, bg='grey')
        self.buttonHEO = Tk.Button(self.frameCtrls, text='HEO',
                                   command=self.displayHEO, bg='grey')
        self.buttonLimitUp.pack(side='top', fill=Tk.X)
        self.buttonLimitDown.pack(side='top', fill=Tk.X)
        self.buttonHEO.pack(side='top', fill=Tk.X)
        self.protocol("WM_DELETE_WINDOW", self.close)
        self.framePlot = Tk.Frame()
        self.framePlot.pack(side='left', fill=Tk.BOTH, expand=1)
        # ----------------------------------------------------------------- #
        #   Uncomment the following lines for testing L2pRadar with fake    #
        # data input read from a file. See dataFakeRead() and dataFakeGen() #
        # ----------------------------------------------------------------- #
        #self.datafile = open('fakeoutput.out', 'r')
        #self.datafile.seek(0, 2)
        #self.pos = self.datafile.tell()        
        self.setFig()
        self.run(newcon=True)
        
    def setFig(self):
        '''Sets figure up'''
        self.fig1 = plt.figure(facecolor='black', figsize=(6, 6))
        self.canvas = FigureCanvasTkAgg(self.fig1, master=self.framePlot)
        self.canvas.get_tk_widget().pack(fill=Tk.BOTH, expand=1)
        self.ax = self.fig1.add_subplot(111, projection='polar')
        plt.subplots_adjust(bottom=0.03, top=0.97, left=0.03, right=0.97)
        self.ax.set_axis_bgcolor('black')
        self.ax.spines['polar'].set_color('white')
        self.ax.grid(color='white', lw=2)
        self.ax.set_theta_direction(-1)
        self.ax.set_theta_offset(-np.pi)
        self.ax.set_yticks(range(0, 90, 10))
        self.ax.set_yticklabels([''] +  map(str, range(80, 0, -10)))
        for label in self.ax.get_xticklabels() + self.ax.get_yticklabels():
            label.set_color('white')
        self.lines = sum((self.ax.plot([], [], lw=4, markeredgewidth=0)
                          for n in range(15)), [])
        self.points = sum((self.ax.plot([], [], 'o', markeredgewidth=0, ms=6,
                           color='w') for n in range(15)), [])
        self.tel_line = self.ax.plot([], [], 'o', color='#00ff00', ms=10)
        self.sun_line = self.ax.plot([], [], 'o', color='#ffff00', ms=25, alpha=0.9)
        #self.sunav_line = self.ax.plot([], [], color='#ffff00')
        self.yhigh = 80
        self.ax.set_ylim(0, self.yhigh)
        
        self.heos = self.ax.plot([], [], 'o', color='red')
    
    def plotLimitUp(self):
        self.yhigh = self.yhigh - 10
        self.ax.set_yticks(range(0, self.yhigh, 10))
        self.ax.set_ylim(0, self.yhigh)
        self.anim._stop()
        self.appRunning = False
        self.run(newcon=False)
        
    def plotLimitDown(self):
        self.yhigh = self.yhigh + 10
        self.ax.set_yticks(range(0, self.yhigh, 10))
        self.ax.set_ylim(0, self.yhigh)    
        self.anim._stop()
        self.appRunning = False
        self.run(newcon=False)
        
    def displayHEO(self):
        """Toggle HEO visibility variable"""
        self.visHEO = not self.visHEO
        print('HEO visibility set to {}'.format(self.visHEO))
        if self.visHEO is True:
            self.buttonHEO.configure(bg='green', activebackground='green')
            self.plotHEO()
            self.anim._stop()
            self.appRunning = False
            self.run(newcon=False)    
        else:
            self.buttonHEO.configure(bg='grey', activebackground='grey')
    
    def plotHEO(self):
        pp2.getPlist(tmpath=self.tmpath)            
        with open(os.path.join(self.tmpath, 'plist.out')) as f:
            heo_list = [line for line in f if line.split()[1][:2].upper() in ['GL', 'CO', 'GI', 'ET', 'GA', 'GP']]
        utcnow = dt.datetime.utcnow()
        epoch = utcnow.second + 60 * (utcnow.minute + 60 * utcnow.hour)
        gazs, gels = [], []
        for line in heo_list:
            npass = line[1:4]
            function = os.path.join(self.tmpath, 'function.' + npass)
            if not os.path.exists(function):
                sp.ftpit('function.' + npass, path='pred')
                os.renames('function.' + npass, function)
            gazelr = fp.azelsat(npass, epoch, fun_path=self.tmpath)
            gazs.append(gazelr[0, 0])
            gels.append(90 - gazelr[0, 1] * 180 / np.pi)
        self.heos[0].set_data(gazs, gels)        
        
    
    def close(self):
        """Closes GUI application and worker subprocess"""
        self.l2p_conn.procWorker.terminate()
        self.destroy()
        print '\nExiting...\n'
        sys.exit()
        
    def anim_init(self):
        """Initial plot state"""
        #self.ax.set_ylim(0, self.yhigh)
        for line, point in zip(self.lines, self.points):
            line.set_data([], [])
            point.set_data([], [])
        self.tel_line[0].set_data([], [])
        self.sun_line[0].set_data([], [])
        self.heos[0].set_data([], [])
        return (self.lines + self.tel_line + self.sun_line + 
                self.points + self.heos)

        
    def animate(self, i):
        """Matplotlib animation function
        """
        # ----------------------------------------------------------------- #
        #       Uncomment the next two lines for off-line testing           #
        # ----------------------------------------------------------------- #
        #planeLines, self.pos = dataFakeRead(self.datafile, self.pos).next()
        #self.P = addPlanes(planeLines, self.P, minel=1)
        
        planeLines, telLines = self.l2p_conn.dump_queue()
        self.P = addPlanes(planeLines, self.P, minel=10)
        
        x = [p.az[-60:] for p in self.P.values()[:15]]
        y = [p.el[-60:] for p in self.P.values()[:15]]
        colours = ((p.el[-1] + 5) / 100 for p in self.P.values())

        try:
            # 56692  41847.094 telscp  75.00  65.00 1
            telPos = telLines[0].split()[3:5]
        except IndexError:
            self.telAz = self.telEl = 0
            #print '\n' + '#' * 40
            #print sys.exc_info()[0]
            #print '#' * 40 + '\n'
        else:
            self.telAz, self.telEl = float(telPos[0]), float(telPos[1][:4])
            
        telAz = float(self.telAz * np.pi / 180)
        telEl = float(90 - self.telEl)
        if i % 10 == 0:
            sunAz, sunEl = noaasun.sunpos(JD='now')
            sunAz = sunAz * np.pi / 180
            self.sun_line[0].set_data(sunAz, 90 - sunEl)
            #t = np.arange(0, 2 * np.pi, 0.2)
            #x = sunAz
            #self.sunav_line[0].set_data(t, sunEl + 15 * np.cos(t))
        nplanes = len(x)
        if nplanes > 0:
            self.tel_line[0].set_data(telAz, telEl)
            for j, line, point in (
                        zip(range(len(self.lines)), self.lines, self.points)):
                line.set_data(x[j], y[j])
                line.set_color(cm.jet_r(colours.next()))
                point.set_data(x[j][-1], y[j][-1])
                if j == nplanes - 1:
                    for point, line in (
                            zip(self.points[nplanes:], self.lines[nplanes:])):
                        line.set_data([], [])
                        point.set_data([], [])
                    break
        else:
            self.lines[0].set_data([], [])
            self.points[0].set_data([], [])
        
        if (self.visHEO is True) and (i % 15 == 0):
            self.plotHEO()
        elif self.visHEO is False:
            self.heos[0].set_data([], [])
            
        sys.stdout.flush()
        
        return (self.lines + self.tel_line + self.sun_line + 
                self.points +  self.heos)
            
    def run(self, newcon=False):
        """Matplotlib animation loop"""
        if self.appRunning is True:
            return
        self.appRunning = True
        if newcon is True:
            self.l2p_conn = ConnectionL2P(self.l2p_HOST, 2020)
        self.anim = animation.FuncAnimation(self.fig1, self.animate,
               init_func=self.anim_init, blit=True, interval=1000, repeat=True)
        self.canvas.draw()
        

def main():
    # First, some house keeping with a hammer
    if (sys.platform == 'linux') or (sys.platform == 'linux2'):
        p = subprocess.Popen(['ps', 'aux'], stdout=subprocess.PIPE)
        ps, err = p.communicate()
        ps = ps.splitlines()
        lines = [l for l in ps if l.find('l2pGUI.py') != -1]
        pids = [l.split()[1] for l in lines]
        pids = pids[:-3]
        for pid in pids:
            subprocess.call(['kill', pid])
    
    # Start now
    app = L2pRadar(None)
    app.mainloop()
    
if __name__ == "__main__":
    main()




