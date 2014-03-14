===========
  l2pGUI
===========

l2pGUI displays ADS/B aircraft beacons [1]. For real-time use it needs 
to connect to l2pserver [2], another program that collects the data 
from the hardware receiver. For offline use it simply reads previously 
collected data from a file.

l2pGUI is essentially a Matplotlib polar plot animation embedded in 
Tkinter with some extra controls (e.g. zoom level). The animation 
performs sufficiently well and is low on resources at the refresh rates 
that make practical sense for this application (i.e. 1-2 FPS). It is 
possible to increase the frame rate when working offline to obtain nice 
animations, achieving tens of frames per second at the expense of high 
CPU load.

When called with no arguments, l2pGUI will attempt to connect to a 
running instance of l2pserver. The IP address and port of the server are 
specified in the configuration file l2pGUI.cfg, whose location depends on 
the platform (probably /usr/l2pGUI/ on Linux). The coordinates 
of the observing station are needed to display the position of the Sun 
correctly, and can be entered in the same configuration file.

Requirements: Python (tested with 2.7), Matplotlib

Install with pip: 

pip install l2pGUI-master.zip

Run with the command l2pgui_run.py

Or extract contents, modify and copy conf/l2pGUI.cfg to the working
directory and run from there:

python l2pGUI.py



listen2planes display client

optional arguments:

  -h, --help            show this help message and exit
  -r REPLAY, --replay REPLAY
                        Replay plane data from file
  -d DUMP2FILE, --dump2file DUMP2FILE
                        Write plane data to file
  -pl, --print-lines    Print data lines
  -t TIME_STEP, --time-step TIME_STEP
                        Time in milliseconds between animation steps
                       

License: GPLv2
Disclaimer:  http://www.bgs.ac.uk/downloads/softdisc.html

[1] http://en.wikipedia.org/wiki/Automatic_dependent_surveillance-broadcast
[2] https://github.com/matwiNERC/listen2planes

