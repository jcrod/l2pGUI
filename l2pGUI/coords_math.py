#!/usr/bin/env python
'''Functions to rotate from geodetic coordinates to Az/El.

These functions are identical to those in module coords.py, 
except that here the calculations are performed at the vector 
component level so that we only use math operations (instead of 
Numpy)
'''
import math
from math import sin
from math import cos


D2R = math.pi / 180


def geod2geo(Lat=50.867387222, Lon=0.33612916666, H=75.357):
    '''Compute geocentric cartesian coordinates from geodetic
    coordinates in WGS84 ellipsoid.
        
    From the Explanatory Supplement to the Astronomical Almanac 1992
    
    Parameters
    ----------
    Lat, Lon, H: latitude and longitude (degrees) and height (m)
    
    Returns
    -------
    x, y, z: cartesian geocentric coordinates vector components (m)
    '''
    Lat = Lat * D2R
    Lon = Lon * D2R
    a = 6378137.0
    f = 1 / 298.257223563
    b = a * (1 - f)
    e2 = (a**2 - b**2) / a**2
    N = a / math.sqrt(1 - e2 * sin(Lat)**2)

    x = (N + H) * cos(Lat) * cos(Lon)
    y = (N + H) * cos(Lat) * sin(Lon)
    z = (N * (1 - e2) + H) * sin(Lat)
    return x, y, z


def geo2top(x, y, z, Lat=50.867387222, Lon=0.33612916666):
    '''Rotate geocentric vector xyz to alt-az frame.
   
    Two steps: l degrees rotation about Z
               90 - p degrees rotation about Y'
                   
    Two equivalent rotation matrices for the second rotation can
    be used (rotating 90 - p or p). The local frame is defined with
    the Z'' axis pointing to the zenith, and the X'' axis to the
    North, therefore a reflection of the X'' axis over the Y''Z''
    plane is needed at the end (this is to conform with the
    astronomical convention of setting the azimuth origin in the
    North.
    '''
    Lat = Lat * D2R
    Lon = Lon * D2R
          
    clon = cos(Lon)
    slon = sin(Lon)
    clat = cos(Lat)
    slat = sin(Lat)
    
    v1 = slat*clon*x + slat*slon*y -clat*z
    v2 = -slon*x +clon*y
    v3 = clat*clon*x + clat*slon*y + slat*z
    return -v1, v2, v3


def top2azel(x, y, z):
    '''Topocentric rectangular coordinates to azimuth and elevation
    '''
    r = math.sqrt(x**2 + y**2 + z**2)
    e = math.atan2(z, math.sqrt(x**2 + y**2))
    _, a = divmod(math.atan2(y, x), 2 * math.pi)
    return a, e, r


def geo2azel(x, y, z, Lat=50.867387222, Lon=0.33612916666, H=75.357):
    '''Geocentric body-fixed (m) to azimuth and elevation
    '''
    X, Y, Z = geod2geo(Lat, Lon, H)
    # Translate to topocentric origin
    x = x - X
    y = y - Y
    z = z - Z
    # Rotate to topocentric frame
    x, y, z = geo2top(x, y, z, Lat, Lon)
    a, e, r = top2azel(x, y, z)
    return a, e, r
    
