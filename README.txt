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
that make practical sense for this application, i.e. 1-2 FPS. It is 
possible to increase the frame rate when working offline, achieving tens 
of frames per second at the expense of high CPU load.

When called with no arguments, l2pGUI will attempt to connect to a 
running instance of l2pserver.


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
                       

[1] http://en.wikipedia.org/wiki/Automatic_dependent_surveillance-broadcast
[2] 

Urls are http://like.this and links can be
written `like this <http://www.example.com/foo/bar>`_.

