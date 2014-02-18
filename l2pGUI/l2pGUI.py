#!/usr/bin/env python
'''Real-time display for ADS/B data. This application implements a client 
for listen2planes, an application that receives data from an ADS/B box 
and sends it via TCP/IP upon request. Data is displayed on a polar plot 
as traces for each individual plane. It can also be used offline by reading 
a file of collected data, which can be re-played at many times the regular 
collection speed.

The application needs to know the IP address of the l2planes server, and 
the coordinates of the observing station for correct determination of the 
Sun position. These can be set below.
'''
import sys, os
import numpy as np
import Tkinter as Tk
import socket
import multiprocessing
#import subprocess
import signal
import argparse
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import time
import datetime as dt
import noaasun
#import funplot as fp
#import plot_plist2 as pp2
#import satpar as sp


L2P_HOST = ('193.61.194.29', 2020)
LON = 0.3361
LAT = 50.8674


def dataFakeRead(f, init_pos=None, N_lines=140, print_lines=False):
    """Reads lines from file from specified position to end.
    
    Let f be a saved l2planes output; calling this function at 
    regular intervals will simulate real time data.
    
    At the moment this function has no concept of time, it simply reads
    a number of lines specified by parameter N_lines. This means that
    the display speed will only be similar to reality if the data packets 
    were detected at similar rates. Thus, one can speed up the replayed 
    animation by increasing this number, at the cost of increased choppiness.
    
    This function is a generator which behaves similarly to Unix tail,
    but the caller is responsible for keeping track of file position
    and passing it as an argument.
    """
    if init_pos is None:
        f.seek(0, 2)
        init_pos = f.tell()
    f.seek(init_pos)
    plines, tlines = [], []
    i = 0
    while True:
        line = f.readline()
        if print_lines:
            print('{}\n'.format(line))
        if not line:
            print 'EOF ' + '%d lines read' % i
            pos = f.tell()
            yield plines, tlines, pos
        l = line.split()
        if len(l) == 13:
            plines.append(line)
        elif len(l) == 6:
            tlines.append(line)
        i += 1
        if i > N_lines:
            break
    pos = f.tell()
    yield plines, tlines, pos


def loadPlanesFile(fname, minel=-5):
    """Loads planes data from file. Useful for analysis.
    
    Parameters
    ----------
    fname: l2planes dump file
    minel: elevation cutoff
    
    Returns
    -------
    P: dictionary containing Plane instances
    """
    t0 = time.time()
    with open(fname, 'r') as f:
        P = {}
        for line in f:
            if len(line.split()) < 13:
                continue
            P = addPlanes((line,), P, minel=-5)
    t = time.time() - t0
    print('{} planes loaded in {:<4.2f} seconds'.format(len(P), t))
    return P

    
def colourMaplimits(value_limits=(0, 80), colourmap_limits=(0.05, 0.85)):
    """Compute a and b coefficients that bring values from value_limits
    to colourmap_limits in a linear way.
    
    I.e., simply solve y = (x - a) / b
    
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
    """Makes planes.
    
    This class takes plane lines from l2planes and stores the relevant
    information for later use
    """
    def __init__(self, line, minel=10):
        l = line.split()
        self.minel = minel
        self.mjd = []
        self.id = l[2]
        self.code = l[3]
        self.epc = []
        self.last_epoch = float(l[1])
        self.lat = []
        self.lon = []
        self.alt = []
        self.ran = []
        self.az = []
        self.el = []
        self.maxel = -10    # maximum observed plane elevation (starting value)
        self.gaps = 0       # times the same plane id has been observed - 1
        self.addLine(l)

    # 56395 40400.326   4ca626 RYR8JT   50.97158 -0.61729 29525 68.6683
    # 280.17873692 7.16197030   -474.0 191.0 1088   0.00 0.00

    def addLine(self, l):
        """Adds one line of data"""
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
            self.mjd.append(float(l[0]))
            self.last_epoch = self.epc[-1]
            self.lat.append(float(l[4]))
            self.lon.append(float(l[5]))
            self.alt.append(float(l[6]) * 0.3048)
            self.ran.append(float(l[7]))
            self.az.append(np.pi / 180 * float(l[8]))
            self.el.append(float(l[9]))
            self.maxel = self.el[-1] if self.el[-1] > self.maxel else self.maxel


def addPlanes(planeLines, planes_dict, minel=-5, time_alive=-1):
    """Processes plane data lines and updates planes dictionary accordingly.
    
    Parameters
    ----------
    planeLines: list of plane lines from l2planes
    planes_dict: dictionary storing planes
                 keys: plane id; values: Plane instances
    minel: elevation cutoff
    time_alive: seconds to wait before discarding planes for which
                no beacons have been received. No limit if set to negative
    """
    P = planes_dict
    for line in planeLines:
        l = line.split()
        plane_id = l[2]
        el = float(l[9])
        if el > minel:
            if not P.has_key(plane_id):
                P[plane_id] = Plane(line, minel)
            else:
                P[plane_id].addLine(l)

    if len(P) == 0 or len(planeLines) == 0:
        return P
    # Remove planes for which no beacons have been 
    # received for more than given time
    if time_alive >= 0:
        last_epoch = float(l[1])
        keys = [k for k,v in P.iteritems() if 
                                 abs(last_epoch - v.last_epoch) > time_alive]
        for key in keys:
            del P[key]
    return P


class L2pRadar(Tk.Tk):
    """Real-time polar plot of ADS-B planes data received from listen2planes
    
    Parameters
    ----------
    replay: if specified, listen2planes data previously written 
            to this file will be displayed
    dump2file: if True, collected data will be written to a file
    """
    def __init__(self, replay=None, dump2file=None, print_lines=None, Tstep=1000):
        Tk.Tk.__init__(self)
        self.replay = replay
        self.dump2file = dump2file
        self.print_lines = print_lines
        self.Tstep = Tstep
        self.MaxPlanes = 25
        self.tmpath = os.path.expanduser('~/.plotsched_tmp')
        self.visHEO = False
        self.P = {}
        self.last_mjd = 0
        
        self.root = Tk.Tk._root(self)
        self.root.configure(background='black')
        self.root.title('l2pGUI')
        self.frameCtrls = Tk.Frame()
        self.frameCtrls.pack(side='left')
        self.buttonLimitUp = Tk.Button(self.frameCtrls, text='Up',
                                       command=self.plotLimitUp, bg='grey')
        self.buttonLimitDown = Tk.Button(self.frameCtrls, text='Down',
                                         command=self.plotLimitDown, bg='grey')
        self.buttonRotate = Tk.Button(self.frameCtrls, text='Rot',
                                      command=self.plotRotate, bg='grey')
        self.buttonHEO = Tk.Button(self.frameCtrls, text='HEO',
                                   command=self.displayHEO, bg='grey')
        self.buttonQuit = Tk.Button(self.frameCtrls, text='Quit',
                                    command=self.close, bg='grey')
        self.buttonLimitUp.pack(side='top', fill=Tk.X, pady=2)
        self.buttonLimitDown.pack(side='top', fill=Tk.X, pady=2)
        self.buttonRotate.pack(side='top', fill=Tk.X, pady=2)
        self.buttonHEO.pack(side='top', fill=Tk.X, pady=2)
        self.buttonQuit.pack(side='top', fill=Tk.X, pady=2)
        self.protocol("WM_DELETE_WINDOW", self.close)
        self.framePlot = Tk.Frame(self.root)
        self.framePlot.pack(side='left', fill=Tk.BOTH, expand=1)
        
        # Read data from file if requested...
        if self.replay:
            self.datafile = open(self.replay, 'r')
            self.pos = self.datafile.tell()
            self.setFig()
            self.run(newcon=False)
        # otherwise proceed normally
        else:
            self.setFig()
            self.run(newcon=True)
            if self.dump2file:
                self.outFile = open(self.dump2file, 'w')

    def setFig(self):
        """Sets figure up"""
        self.fig1 = Figure(facecolor='black', figsize=(6, 6))
        self.canvas = FigureCanvasTkAgg(self.fig1, master=self.framePlot)
        self.canvas.get_tk_widget().pack(fill=Tk.BOTH, expand=1)
        self.ax = self.fig1.add_subplot(111, projection='polar')
        self.fig1.subplots_adjust(bottom=0.03, top=0.97, left=0.03, right=0.97)
        self.ax.set_axis_bgcolor('black')
        self.ax.spines['polar'].set_color('white')
        self.ax.grid(color='white', lw=2)
        self.ax.set_theta_direction(-1)
        self.theta_offset = 2
        self.ax.set_theta_offset(self.theta_offset * np.pi/2)
        self.ax.set_yticks(range(0, 90, 10))
        self.ax.set_yticklabels([''] +  map(str, range(80, 0, -10)))
        for label in self.ax.get_xticklabels() + self.ax.get_yticklabels():
            label.set_color('white')
        self.lines = sum((self.ax.plot([], [], lw=4, markeredgewidth=0)
                          for n in range(self.MaxPlanes)), [])
        self.points = self.ax.plot([], [], 'o', markeredgewidth=0, 
                                                   ms=6, color='w')
        self.tel_line = self.ax.plot([], [], 'o', color='#00ff00', ms=10)
        self.sun_line = self.ax.plot([], [], 'o', color='gold',
                                     ms=25, alpha=0.8)
        self.sunav_line = self.ax.plot([], [], color='gold', lw=2, alpha=0.8)
        self.yhigh = 90
        self.ax.set_ylim(0, self.yhigh)
        self.heos = self.ax.plot([], [], 'o', color='red')
    
        self.time = time.time()
    
    def plotLimitUp(self):
        """Decrease plot elevation range"""
        self.yhigh = self.yhigh - 10
        self.ax.set_yticks(range(0, self.yhigh, 10))
        self.ax.set_ylim(0, self.yhigh)
        # Since the animation is blitted we need to restart it or
        # the changes above will only last one animation cycle
        self.anim._stop()
        self.run(newcon=False)
        
    def plotLimitDown(self):
        """Increse plot elevation range"""
        self.yhigh = self.yhigh + 10
        self.ax.set_yticks(range(0, self.yhigh, 10))
        self.ax.set_ylim(0, self.yhigh)
        self.anim._stop()
        self.run(newcon=False)
        
    def plotRotate(self):
        """Rotate plot"""
        self.theta_offset = np.mod(self.theta_offset + 1, 4)
        self.ax.set_theta_offset(self.theta_offset * np.pi / 2)
        dirs = {0:'RIGHT', 1:'TOP', 2:'LEFT', 3:'BOTTOM'}
        print('\nNORTH set to {}\n'.format(dirs[self.theta_offset]))
        self.anim._stop()
        self.run(newcon=False)
        
    def displayHEO(self):
        """Toggle HEO visibility variable"""
        self.visHEO = not self.visHEO
        print('HEO visibility set to {}'.format(self.visHEO))
        if self.visHEO is True:
            self.buttonHEO.configure(bg='green', activebackground='green')
            self.plotHEO()
            self.anim._stop()
            self.run(newcon=False)    
        else:
            self.buttonHEO.configure(bg='grey', activebackground='grey')
    
    def plotHEO(self):
        """Plot predicted HEO satellites"""
        pp2.getPlist(tmpath=self.tmpath)            
        with open(os.path.join(self.tmpath, 'plist.out')) as f:
            heo_list = [line for line in f if line.split()[1][:2].upper() in
                        ['GL', 'CO', 'GI', 'ET', 'GA', 'GP']]
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

    def anim_init(self):
        """Initial plot state"""
        #for line, point in zip(self.lines, self.points):
            #line.set_data([], [])
            #point.set_data([], [])
        for line in self.lines:
            line.set_data([], [])
        for line in [self.tel_line, self.sun_line, self.heos, self.sunav_line, self.points]:
            line[0].set_data([], [])
        return (self.lines + self.tel_line + self.sun_line + 
                self.points + self.heos  + self.sunav_line)
    
    def el2zdist(self, x):
        """Elevation to zenith distance (degrees)"""
        return 90 - x
    
    def formattedOutput(self):
        """Prints planes present in the queue"""
        # hex      id       Az   El     Lon      Lat       Alt   Dist
        #-------------------------------------------------------------        
        #a2b728  UPS203    291 15.4    -0.199   50.994    11278  41.7
        #aa7974  SOO275    343  5.1    -0.101   51.754    10058 103.8
        print(' hex      id       Az  El      Lon     Lat        Alt   Dist')
        print '-' * 62
        for plane in self.P.values():
            strf = '{:8s}{:8s} {:4.0f}{:5.1f} {:>9.3f}{:>9.3f} {:>8.0f}{:>6.1f}' 
            print(strf.format(plane.id, plane.code, 
                              plane.az[-1] * 180 / np.pi, plane.el[-1],
                              plane.lon[-1], plane.lat[-1],
                              plane.alt[-1], plane.ran[-1]))
        print('')
        sys.stdout.flush()
            
    def process_lines(self, data_lines, print_lines=False, dump2file=False):
        """Processes data lines according to length
        
        Parameters
        ----------
        data_lines: list of data lines
        print_output: boolean flag to request printed output
        dump: boolean flag to request written output
        
        Returns
        -------
        plines, tlines: lists containing plane and telescope lines
        """
        tlines, plines = [], []
        for line in data_lines:
            if print_lines is True:
                print('{}\n'.format(line))
            if dump2file:
                self.outFile.write(line + ' \n')
            L = len(line.split())
            if L == 13:
                plines.append(line)
            elif L == 6:
                tlines.append(line)
            elif L == 2:
                time.sleep(2)
                self.reconnect()
        return plines, tlines
    
    def updateData(self):
        """Update planes dictionary with data from queue or from dump file"""
        # Grab data via TCP/IP normally...
        if not self.replay:
            data_lines = dump_queue(self.planeQueue)
            planeLines, telLines = self.process_lines(data_lines, 
                                                  print_lines=self.print_lines,
                                                  dump2file=self.dump2file)
            self.P = addPlanes(planeLines, self.P, minel=0, time_alive=15)
        # or read data from dump.out if requested
        elif self.replay:
            planeLines, telLines, self.pos = (
                             dataFakeRead(self.datafile, self.pos,
                                          print_lines=self.print_lines).next())
            self.P = addPlanes(planeLines, self.P, minel=0, time_alive=15)
        
        if len(telLines) > 0:
            self.telLines = telLines[-1]
        else:
            self.telLines = '0 0 0 00.00 00.00 0'
        
    def animate(self, i):
        """Matplotlib animation function
        
        NB here 'lines' refers to plot lines, not data lines
        """
        self.updateData()
        if not self.print_lines:
            self.formattedOutput()

        # Az/El from planes present in the dictionary, grabbing only 
        # the last 80 positions available in steps of 5
        Azs = [p.az[-80::5] for p in self.P.values()]
        Els = [p.el[-80::5] for p in self.P.values()]
        
        colours = ((p.el[-1] + 5) / 100 for p in self.P.values() if len(p.el) > 0)
        Nplanes = len(Azs)
        if Nplanes > 0:
            # Update plot lines with Az/El data
            for j, line in enumerate(self.lines):
                line.set_data(Azs[j], map(self.el2zdist, (Els[j])))
                line.set_color(cm.jet(colours.next()))
                if j == Nplanes - 1:
                    # Reset plot lines from previous planes
                    for line in self.lines[Nplanes:]:
                        line.set_data([], [])
                    break
        else:
            self.lines[0].set_data([], [])
        
        x = [n[-1] for n in Azs]
        y = [n[-1] for n in Els]
        self.points[0].set_data(x, map(self.el2zdist, y))

        
        if (self.visHEO is True) and (i % 15 == 0):
            self.plotHEO()
        elif self.visHEO is False:
            self.heos[0].set_data([], [])
            
        # Telescope position
        # 56692  41847.094 telscp  75.00  65.00 1
        telPos = self.telLines.split()[3:5]
        telAz = float(telPos[0]) * np.pi / 180
        telEl = float(telPos[1][:4])
        self.tel_line[0].set_data(telAz, 90 - telEl)
        
        # Update Sun position every 20 animation steps
        if i % 20 == 0:
            # Print FPS
            #newtime = time.time()
            #print('\nFPS: {:4.1f}\n'.format(20 / (newtime - self.time)))
            #self.time = newtime
            if not self.replay:
                d, JD = noaasun.parseDates(mode='now')
                sunAz, sunEl = noaasun.sunpos(JD, lon=LON, lat=LAT)
            else:
                # Read MJD from data so that the Sun is in the right place
                if len(self.P) > 0:
                    # Update last date if planes found
                    self.last_mjd = (self.P.values()[0].mjd[-1] + 
                                 self.P.values()[0].epc[-1] / 86400)
                sunAz, sunEl = noaasun.sunpos(JD=2400000.5 + self.last_mjd,
                                              lon=LON, lat=LAT)
                
            sunAz = sunAz * np.pi / 180
            self.sun_line[0].set_data(sunAz, 90 - sunEl)
            
            # Draw Sun avoidance region         
            if sunEl > -20:
                radius = 15
                theta = np.linspace(0, 2 * np.pi, 40)
                X = radius * np.cos(theta)
                Y = radius * np.sin(theta)
                A = (sunAz + X * np.pi / 180)            
                B = sunEl + Y
                az_corr = 1 / np.cos(B * np.pi / 180)
                A = sunAz + (A - sunAz) * az_corr
                self.sunav_line[0].set_data(A, 90 - B)
            else:
                self.sunav_line[0].set_data(0, 0)

        return (self.lines + self.tel_line + self.sun_line + 
                self.points +  self.heos + self.sunav_line)
    
    def run(self, newcon=False):
        """Start subprocesses and Matplotlib animation loop"""
        if newcon is True:
            self.planeQueue = multiprocessing.Queue()
            self.procWorker = multiprocessing.Process(target=receive_proc,
                                                      args=[self.planeQueue])
            self.procWorker.start()
            
        self.anim = animation.FuncAnimation(self.fig1, self.animate,
               init_func=self.anim_init, blit=True, interval=self.Tstep)
        
        signal.signal(signal.SIGINT, self.signal_handler)

    def signal_handler(self, signal, frame):
        """Catch SIGINT and close everything properly"""
        print('CTRL-C')
        self.close()
        
    def reconnect(self):
        """Try to re-establish connections"""
        print('\nAttempting to reconnect to l2p server...\n')
        self.planeQueue.close()
        self.planeQueue = multiprocessing.Queue()
        self.procWorker.terminate()
        self.procWorker = multiprocessing.Process(target=receive_proc,
                                                      args=[self.planeQueue])
        self.procWorker.start()
        
    def close(self):
        """Closes application and worker subprocess as appropriate"""
        if not self.replay:
            self.procWorker.terminate()
            self.planeQueue.close()
        if self.dump2file:
            self.outFile.close()
        self.root.destroy()
        print '\nExiting...\n'
        sys.exit()


def receive_proc(planeQueue):
    """Requests data lines from l2planes and sends it to the queue"""
    connSocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        connSocket.connect(L2P_HOST)
    except socket.error, msg:
        print('Error connecting: {}'.format(msg))

    while True:
        # Send string expected by listen2planes from clients
        try:
            connSocket.send('reader\0')
        except socket.error, msg:
            print('Failed to send to server: '.format(msg))
            time.sleep(1.5)
            # We will look for the following error string
            # and attempt to reconnect if found
            planeQueue.put('CONN ERROR')
            break
            
        try:
            data = connSocket.recv(256)
        except socket.error, msg:
            print('Failed to receive from server: {}'.format(msg))
            time.sleep(1.5)
            planeQueue.put('CONN ERROR')
            break
        else:
            planeQueue.put(data)
    return


def dump_queue(planeQueue):
    """Retrieves all the data lines from the queue"""
    data_lines = []
    while True:
        data_lines.append(planeQueue.get())
        if planeQueue.qsize() == 0:
            break
    return data_lines


def main(argv=None):
    """Deal with command line arguments and launch the program"""
    parser = argparse.ArgumentParser(description='listen2planes display client')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-r', '--replay', help='Replay plane data from file')
    group.add_argument('-d', '--dump2file', help='Write plane data to file')
    parser.add_argument('-pl', '--print-lines', action='store_true',
                        help='Print data lines')
    parser.add_argument('-t', '--time-step', type=int, default=1000,
                        help='Time in milliseconds between animation steps')
    args = parser.parse_args()
    
    # Some house keeping with a hammer to clean processes from previous 
    # runs in Linux systems. Not needed anymore since we're properly 
    # terminating child processes. Left here just in case.
    #if (sys.platform == 'linux') or (sys.platform == 'linux2'):
        #p = subprocess.Popen(['ps', 'aux'], stdout=subprocess.PIPE)
        #ps, err = p.communicate()
        #ps = ps.splitlines()
        #lines = [l for l in ps if l.find('l2pGUI.py') != -1]
        #pids = [l.split()[1] for l in lines]
        #pids = pids[:-3]
        #for pid in pids:
            #subprocess.call(['kill', pid])
    
    app = L2pRadar(replay=args.replay, 
                   dump2file=args.dump2file, 
                   print_lines=args.print_lines,
                   Tstep=args.time_step)
    app.mainloop()
    
    
if __name__ == "__main__":
    sys.exit(main())
