#!/usr/bin/env python

import math
from math import sin
from math import cos
import coords_math as cds


D2R = math.pi / 180
R2D = 180 / math.pi
S2R = 1. / 3600 * D2R


def moonpos(JD):
    '''Moon position in geocentric, body-fixed coordinates.
    
    Algorithm from Montenbruck and Gill
    '''
    d = JD - 2451545.0
    T = d / 36525.0

    # These are in degrees
    L0 = 218.31617 + 481267.88088 * T - 1.3972 * T
    l = 134.96292 + 477198.86753 * T
    lp = 357.52543 + 35999.04944 * T
    F = 93.27283 + 483202.01873 * T
    D = 297.85027 + 445267.11135 * T

    l = l * D2R
    lp = lp * D2R
    F = F * D2R
    D = D * D2R
    L0 = L0 * 3600

    # Moon's longitude with respect to equinox and ecliptic J2000 (seconds)
    lambda_m = (L0 + 22640 * sin(l) + 769 * sin(2*l) -
                4586 * sin(l - 2*D) + 2370 * sin(2*D) -
                668 * sin(lp) - 412 * sin(2*F) - 
                212 * sin(2*l - 2*D) - 206 * sin(l + lp - 2*D) +
                192 * sin(l + 2*D) - 165 * sin(lp - 2*D) +
                148 * sin(l - lp) - 125 * sin(D) -
                110 * sin(l + lp) - 55 * sin(2*F - 2*D) )

    _, lambda_m = divmod(lambda_m / 3600, 360)
    lm = lambda_m * D2R
    L0 = L0 / 3600 * D2R

    # Moon's latitude (seconds)
    beta_m = (
        18520 * sin(F + lm - L0 + 412*S2R * sin(2*F) + 541*S2R * sin(lp)) -
        526 * sin(F - 2*D) + 44 * sin(l + F - 2*D) -
        31 * sin(-l + F - 2*D) - 25 * sin(-2*l + F) -
        23 * sin(lp + F - 2*D) + 21 * sin(-l + F) +
        11 * sin(-lp + F - 2 *D) )
            
    # Moon's distance (km)
    rM = (385000 - 20905 * cos(l) - 3699 * cos(2*D - l) - 
        2956 * cos(2*D) - 570 * cos(2*l) + 246 * cos(2*l - 2*D) -
        205 * cos(lp - 2*D) - 171 * cos(l + 2*D) -
        152 * cos(l + lp - 2*D)) * 1000
    
    # Convert ecliptic to equatorial Cartesian coordinates
    lm = lm + (1.3972 * T) * D2R

    x = rM * cos(lm) * cos(beta_m * S2R)
    y = rM * sin(lm) * cos(beta_m * S2R)
    z = rM * sin(beta_m * S2R)

    ecl = 23.43929111 * D2R
    s = sin(-ecl)
    c = cos(-ecl)
    y, z = (y * c + z * s, -y * s + z * c)

    # Rotate to body-fixed
    theta = math.radians(280.46061837 + 360.98564736629 * d)
    s = math.sin(theta)
    c = math.cos(theta)
    x, y  = (x * c + y * s, -x * s + y * c)
    return x, y, z

    
def sunpos_lacc(JD):
    '''Sun position (km). Low-precision algorithm from Montenbruck and Gill
    '''
    d = (JD - 2451545.0)
    T = d / 36525.0
    
    # Omega + omega = 282.94 degrees
    M = 357.5256 + 35999.049 * T
    l = 282.94 + M + 6892./3600 * sin(M * D2R) + 72./3600 * sin(2 * M * D2R)
    r = (149.619 - 2.499 * cos(M * D2R) - 0.021 * cos(2 * M * D2R)) * 1e6
    ecl = 23.4392911
    
    x = r * cos(l * D2R)
    y = r * sin(l * D2R) * cos(ecl * D2R)
    z = r * sin(l * D2R) * sin(ecl * D2R)

    theta = (280.46061837 + 360.98564736629 * (JD - 2451545.0) +
             0.000387933 * T**2 - T**3 / 38710000)
    sint = math.sin(theta * D2R)
    cost = math.cos(theta * D2R)
    x, y = (x*cost + y*sint, y*cost - x*sint)
    return x, y, z


def sunpos(JD=2448908.5):
    '''Sun position. Jean Meeus, Astronomical Algorithms, 1998, chs 25 & 26.
    
    Stated accuracy: 0.01 degrees in RA/DEC
    '''
    # Time is Julian centuries from J2000.0
    T = (JD - 2451545.0) / 36525
    
    # Geometric mean longitude (mean equinox of date)
    L0 = 280.46646 + 36000.76983 * T + 0.0003032 * T**2
    # Mean anomaly
    M = 357.52911 + 35999.05029 * T - 0.0001537 * T**2
    # Eccentricity of Earth's orbit
    e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T**2
    
    # Sun's equation of centre
    C = ((1.914602 - 0.004817 * T - 0.000014 * T**2) * sin(M * D2R) +
         (0.019993 - 0.000101 * T) * sin(2 * M * D2R) + 
         0.000289 * sin(3 * M * D2R) )
    
    # Sun's true geometric longitude (mean equinox of date)
    Theta = L0 + C
    # True anomaly
    v = M + C
    
    
    # Sun's radius vector (m) (eq. 25.5 times 1 A.U)
    R = 1.000001018 * (1 - e**2) / (1 + e * cos(v * D2R)) * 1.49597870e11

    # Apparent (true of date) longitude (nutation and aberration corrections)
    Omg = 125.04 - 1934.136 * T
    lon = Theta - 0.00569 - 0.00478 * sin(Omg * D2R)
    beta = 0    # never greater than 1''.2
    
    # Mean obliquity of the ecliptic (eq. 22.2)
    eps = (23 + (26 + (21.448 / 60)) / 60 - 
           46.8150 / 3600 * T -
           0.00059 / 3600 * T**2 +
           0.001813/ 3600 * T**3)
    
    # Right ascension and declination
    #RA = math.atan2(cos(eps * D2R) * sin(Theta * D2R), cos(Theta * D2R))
    #DEC = math.asin(sin(eps * D2R) * sin(Theta * D2R))
    # If apparent RA/DEC is required, use lon instead of Theta,
    # and correct eps by the quantity 0.00256 * cos(Omg) (degrees)
    
    # Longitude and latitude to rectangular coordinates
    x = R * cos(lon * D2R)
    y = R * sin(lon * D2R) * cos(eps * D2R)
    z = R * sin(lon * D2R) * sin(eps * D2R)
    
    # Mean sidereal time (degrees) (eq. 12.4)
    the = (280.46061837 + 360.98564736629 * (JD - 2451545.0) +
           0.000387933 * T**2 - T**3 / 38710000)
    
    # Rotate to body-fixed coordinates
    sint = sin(the * D2R)
    cost = cos(the * D2R)
    x, y = (x*cost + y*sint, y*cost - x*sint)
    return x, y, z


def sunazel_lacc(JD, Lat=50.867387222, Lon=0.33612916666, H=75.357):
    '''Sun's azimuth and elevation at given JD and geodetic location.
    
    Low-precision algorithm'''
    x, y, z = sunpos_lacc(JD)
    az, el, r = cds.geo2azel(x, y, z, Lat, Lon, H)
    return az, el, r    
    
    
def sunazel(JD, Lat=50.867387222, Lon=0.33612916666, H=75.357):
    '''Sun's azimuth and elevation at given JD and geodetic location'''
    x, y, z = sunpos(JD)
    az, el, r = cds.geo2azel(x, y, z, Lat, Lon, H)
    return az, el, r        


def test_sunazel(JD=2448908.5):
    a1, e1, r1 = sunazel_lacc(JD)
    a2, e2, r2 = sunazel(JD)
    
    print a1 * R2D, e1 * R2D
    print a2 * R2D, e2 * R2D
    

def moonazel(JD=2448908.5, Lat=50.867387222, Lon=0.33612916666, H=75.357):
    x, y, z = moonpos(JD)
    az, el, r = cds.geo2azel(x, y, z, Lat, Lon, H)
    return az, el, r








