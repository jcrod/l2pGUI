#!/usr/bin/env python
'''This module provides some commonly used space and time transformations.

It also bundles some convenient vector and matrix operations to perform 
calculations on arrays efficiently.
'''
from __future__ import division
import numpy as np
import datetime as dt

class Transf(): 
    '''Several coordinate transformations packed in a class.
    
    The geodetic coordinates of the place of interest are specified 
    at instantiation. Else, Herstmonceux coordinates are chosen 
    as a default.
    
    Initial Class Parameters
    ------------------------
    Lon, Lat: geodetic longitude and latitude (degrees)
    h: height(m)
    ellipsoid: reference ellipsoid, e.g. WGS84 (default)
    '''
    def __init__(self, Lat=50.867387222, Lon=0.33612916666, H=75.357,
                 ellipsoid='WGS84'):

        if ellipsoid == 'WGS84':             
            a = 6378137.0
            f = 1 / 298.257223563
        if ellipsoid == 'IAU76':
            a = 6378140.0
            f = 1 / 298.257
        if ellipsoid == 'made up':
            a = 6378137.0
            f = 1 / 298.257
        self.Lon = np.double(Lon * np.pi / 180)
        self.Lat = np.double(Lat * np.pi / 180)
        self.H = np.double(H)
        self.f = np.double(f)
        self.a = a
        self.b = a * (1 - f)
        self.e2 = (a**2 - self.b**2) / a**2
        self.N = a / np.sqrt(1 - self.e2 * np.sin(self.Lat)**2)

        #print Lat, Lon, H
        
    def geod2geo(self):
        '''Compute geocentric cartesian coordinates from the geodetic
        coordinates specified during class initialization.
        
        From the Explanatory Supplement to the Astronomical Almanac 1992
        
        Returns
        -------
        [x, y, z]: cartesian geocentric coordinates vector (m)
        '''
        x = (self.N + self.H) * np.cos(self.Lat) * np.cos(self.Lon)
        y = (self.N + self.H) * np.cos(self.Lat) * np.sin(self.Lon)
        z = (self.N * (1 - self.e2) + self.H) * np.sin(self.Lat)
        #print x
        #print y
        #print z
        return np.array([x, y, z])


    def vgeo2top(self, v):
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
        Lon = self.Lon
        Lat = self.Lat

        Rlz = np.array(( [ np.cos(Lon), np.sin(Lon),            0] ,
                         [-np.sin(Lon), np.cos(Lon),            0] ,
                         [           0,           0,            1] ))
                          
        Rpy = np.array(( [ np.sin(Lat),           0, -np.cos(Lat)] ,
                         [           0,           1,            0] ,
                         [ np.cos(Lat),           0,  np.sin(Lat)] ))
        
        Rpylz = np.dot(Rpy, Rlz)
        result = np.dot(Rpylz, v.T).T
        
        if sum(np.shape(v)) == 3:
            result[0] = -result[0]
        else:
            result[:, 0] = -result[:, 0]
        return result

        
    def vtop2geo(self, vector):
        '''Rotate from local rectangular frame to geocentric coordinates.
        
        First, a reflection of the X axis over the YZ plane, followed by
        a rotation about the Y and another about Z. In practical terms, 
        it just uses the inverse (i.e. transpose in this case) of the
        rotation matrix employed in transf.geo2top()
        '''
        Lon = self.Lon
        Lat = self.Lat
        v = np.copy(vector)
        if sum(np.shape(v)) == 3:
            v[0] = -v[0]
        else:
            v[:, 0] = -v[:, 0]

        Rlz = np.array(( [ np.cos(Lon), np.sin(Lon),            0] ,
                         [-np.sin(Lon), np.cos(Lon),            0] ,
                         [           0,           0,            1] ))
                          
        Rpy = np.array(( [ np.sin(Lat),           0, -np.cos(Lat)] ,
                         [           0,           1,            0] ,
                         [ np.cos(Lat),           0,  np.sin(Lat)] ))
        
        RpylzT = np.dot(Rpy, Rlz).T
        result = np.dot(RpylzT, v.T).T
        return result
        
        
    def vgeo2azel(self, v):
        '''Body-fixed (XYZ) to azimuth and elevation
        '''
        vherl = self.geod2geo() * 1e-6
        vtop = self.vgeo2top(v - vherl)
        return vtop2azel(vtop)


def vtop2azel(V):
    '''Topocentric rectangular coordinates to azimuth and elevation.
    '''
    if isinstance(V, np.ndarray) is False:
        print '\n(vxyz2azel) ERROR: Input value not a vector.\n'
        return
    if np.shape(V) == (3,):
        V = V[None, :]
    elif np.shape(V) == (3, 1):
        V = V.T
        
    r = np.sqrt(V[:, 0]**2 + V[:, 1]**2 + V[:, 2]**2)
    e = np.arctan2(V[:, 2], np.sqrt(V[:, 0]**2 + V[:, 1]**2))
    a = np.mod(np.arctan2(V[:, 1], V[:, 0]), 2*np.pi)
    return np.array([a, e, r]).T


def vazel2top(Vaer):
    '''Azimuth and elevation to topocentric rectangular coordinates.
    
    ***  Positive X axis pointing North  ***
    '''
    if isinstance(Vaer, np.ndarray) is False:
        print '\n(coords.vazel2top) ERROR: Input value not a vector.\n'
        return
    if np.shape(Vaer) == (3,):
        Vaer = Vaer[None, :]
    elif np.shape(Vaer) == (3, 1):
        Vaer = Vaer.T
        
    x = -Vaer[:, 2] * np.cos(Vaer[:, 1]) * np.cos(Vaer[:, 0])
    y =  Vaer[:, 2] * np.cos(Vaer[:, 1]) * np.sin(Vaer[:, 0])
    z =  Vaer[:, 2] * np.sin(Vaer[:, 1])
    return np.array([x, y, z]).T
    

def vazel2top2(Vaer):
    '''Azimuth and elevation to rectangular coordinates.
    
    *** Positive X axis pointing South ***'''
    if isinstance(Vaer, np.ndarray) is False:
        print '\n(coords.vazel2top2) ERROR: Input value not a vector.\n'
        return
    if np.shape(Vaer) == (3,):
        Vaer = Vaer[None, :]
    elif np.shape(Vaer) == (3, 1):
        Vaer = Vaer.T
    
    x = Vaer[:, 2] * np.cos(Vaer[:, 1]) * np.cos(Vaer[:, 0])
    y = Vaer[:, 2] * np.cos(Vaer[:, 1]) * np.sin(Vaer[:, 0])
    z = Vaer[:, 2] * np.sin(Vaer[:, 1])
    return np.array([x, y, z]).T


def geod2azel(Lat, Lon, Alt, alt_units='m'):
    '''Geodetic coordintes to azimuth and elevation.
    
    Parameters
    ----------
    Lon = geodetic longitude (degrees)
    Lat = geodetic latitude (degrees)
    Alt = altitude (feet)
    
    Returns
    -------
    aer = vector containing azimuth, elevation and range 
          in degrees and meters.
    '''
    if alt_units == 'ft':
        Alt = Alt * 0.3048
    pT = Transf(Lat, Lon, Alt)
    p = pT.geod2geo()

    hT = Transf()
    h = hT.geod2geo()

    v = p - h
    w = hT.vgeo2top(v)
    aer = vtop2azel(w)[0]
    aer[:2] = aer[:2] * 180 / np.pi
    return aer

    
def s2hms(epoch):
    h, m = divmod(epoch, 3600)
    m, s = divmod(m, 60)
    return int(h), int(m), float(s)


def h2hms(h):
    m, h = np.modf(h)
    m = m * 60
    s, m = np.modf(m)
    s = s * 60
    return h, m, s


def hms2s(h, m=0, s=0):
    return ((h * 60) + m) * 60 + s
    

def gcal2jd(Y=1978, M=7, D=27):
    '''From Explanatory Supplement 12.92
    '''
    A = -((14 - M) // 12)
    JD = (1461 * (Y + 4800 + A)) // 4
    JD += (367 * (M - 2 - 12 * A)) // 12
    JD -= 3 * ((Y + 4900 + A) // 100) // 4
    JD += D - 32075.5       # subtract 0.5 to start at 00:00
    return JD

         
def jd2gcal(JD=2443716):
    '''From Jean Meeus' Astronomical Algorithms. It starts at 0000.
    '''
    F, Z = np.modf(JD + 0.5)
    
    if Z >= 2299161:
        alpha = (Z - 1867216.25) // 36524.25
        A = Z + 1 + alpha - alpha // 4
    else:
        A = Z

    B = A + 1524
    C = (B - 122.1) // 365.25
    D = np.modf(365.25 * C)[1]
    E = (B - D) // 30.6001
    D = B - D - np.modf(30.6001 * E)[1] + F
    fD, D = np.modf(D)    
    if E < 14:
        M = E - 1
    else:
        M = E - 13
    if M > 2:
        Y = C - 4716
    else:
        Y = C - 4715
    return int(Y), int(M), int(D), fD
    

def jdToday():
    """Returns JD for current date
    """
    today = dt.date.today()
    jd = gcal2jd(today.year, today.month, today.day)
    print jd
    return
    

def jd2gcal2(JD=2443716):
    '''From Explanatory Supplement. It starts at 12:00.
    '''
    L = JD + 68569
    N = (4 * L) // 146097
    L = L - (146097 * N + 3) //4
    I = (4000 * (L + 1)) // 1461001
    L = L - (1461 * I) // 4 + 31
    J = (80 * L) // 2447
    D = L - (2447 * J) // 80
    L = J // 11
    M = J + 2 - 12 * L
    Y = 100 * (N - 49) + I + L
    return Y, M, D

    
def jcal2jd(Y=1978, M=7, D=27):
    '''From Explanatory Supplement.
    '''
    JD = (367 * Y - (7 *(Y + 5001 - ((9 - M) // 7))) // 4 + 
         (275 * M) // 9 + D + 1729777)
    return JD    

    
def jd2jcal(JD=2443716):
    '''From Explanatory Supplement.
    '''
    J = JD + 1402
    K = (J - 1) // 1461
    L = J - 1461 * K
    N = (L - 1) // 365 - L // 1461
    I = L - 365 * N + 30
    J = (80 * I) // 2447
    D = I - (2447 * J) // 80
    I = J // 11
    M = J + 2 - 12 * I
    Y = 4 * K + N + I - 4716
    return Y, M, D
    
    
def epsilon(JD=2443716):
    '''Mean obliquity
    '''
    a = 23 * 3600 + 26 * 60 + 21.448
    T = (JD - 2451545.0) / 36525
    e0 = a - 46.8150 * T - 0.00059 * T**2 + 0.001813 * T**3
    return e0
    

def gmst(JD=2443716):
    '''Explanatory Supplement, for apparent sidereal time...
    I'll be back to this at some point...
    '''
    Tu = (JD - 2451545) / 36525.
    gmst = (67310.54841 + (87660*3600 + 8640184.812866) * Tu +
            0.093104 * Tu**2 - 6.2e-6 * Tu**3)
    print gmst
    return gmst

    
def jd2gmst(JD=2443716, verbose=False):
    '''Greenwich Mean Sideral Time at JD.
    
    From Explanatory Supplement (2.24), implemented as in :
    http://celestrak.com/columns/v02n02/
    
    Earth's rotation rate (rad/s):
    w = 7.2921151467e-5
    Ratio of UT-day to period of rotation:
    w * 86400 * 180 / np.pi * 1 / 360 = 1.0027378
    '''
    UT = np.modf(JD + 0.5)[0]
    JD = JD - UT
    Tu = (JD - 2451545.0) / 36525.  # centuries of 36525s since J2000
    # GMST at 0h:
    GMST0 = (24110.54841 + 8640184.812866 * Tu + 0.093104 * Tu**2 -
            6.2e-6 * Tu**3)

    # GMST at UT (s):
    GMST = np.mod(GMST0 + 86400 * 1.00273781191135448 * UT, 86400)
    thetag = 2 * np.pi * GMST / 86400  # (rad)
    H, M, S = s2hms(GMST)

    if verbose == True:
        print('GMST: {:<10.4f}s ({:<6.4f} rad)').format(GMST, thetag)
        print('      {:0=2.0f}h {:<2.0f}m {:<6f}s').format( H, M, S)
    return thetag


def zrot(v, t):
    '''Positive rotation about the Z axis.
    '''
    Rg = np.array([ [ np.cos(t), -np.sin(t), 0],
                    [ np.sin(t),  np.cos(t), 0],
                    [         0,          0, 1] ])
    return np.dot(Rg, v.T).T
    
    
def gei2geo(v, JD=2443716):
    '''Inertial coordinates to body-fixed rectangular geocentric.
    '''
    t = jd2gast(JD)
    return zrot(v, -t)

def geo2gei(v, JD):
    '''Body-fixed rectangular geocentric coordintates to inertial.
    '''
    t = jd2gast(JD)
    return zrot(v, t)

    
def DpsiDeps(JD=2443716):
    '''Computes angular corrections to psi and eps
    '''
    d2r = np.pi / 180
    D = JD - 2451545.0
    Dpsi = -0.0048*d2r * np.sin(125.0*d2r - 0.05295*d2r * D) - \
            0.0004*d2r * np.sin(200.9*d2r + 1.97129*d2r * D)
    Deps =  0.0026*d2r * np.cos(125.0*d2r - 0.05295*d2r * D) + \
            0.0002*d2r * np.cos(200.9*d2r + 1.97129*d2r * D)
    return Dpsi, Deps    


def epsilon_0(JD=2443716):
    '''Mean obliquity of date (radians)
    '''
    T = (JD - 2451545.0) / 36525.
    a = 23 * 3600 + 26 * 60 + 21.448
    eps_0 = a - 46.8150 * T - 0.00059 * T**2 + 0.001813 * T**3
    eps_0 = eps_0 / 3600. * np.pi / 180
    return eps_0
    

def jd2gast(JD=2443716, verbose=False):
    '''Julian date to Greenwich apparent sidereal time
    '''
    thetag = jd2gmst(JD, verbose)
    GMST = thetag / (2 * np.pi) * 86400
    Dpsi, Deps = DpsiDeps(JD)
    eps = epsilon_0(JD) + Deps
    GAST = thetag + Dpsi * np.cos(eps)
    theta = GAST
    GAST = GAST * 86400 / (2 * np.pi)
    H, M, S = s2hms(GAST)
    
    if verbose == True:
        print('GAST: {:<10.4f}s ({:<6.4f} rad)').format(GAST, theta)
        print('      {:0=2.0f}h {:<2.0f}m {:<6f}s').format( H, M, S)
    return theta


def mean2true_equinox(V, JD=2443716):
    '''Mean to true equinox of date
    '''
    Dpsi, Deps = DpsiDeps(JD)
    eps = epsilon_0(JD) + Deps
    
    if np.shape(V) == (3,):
        V = V[None, :]
    v = np.zeros_like(V)
    
    Dx = -(V[:,1] * np.cos(eps) + V[:,2] * np.sin(eps)) * Dpsi
    Dy =   V[:,0] * np.cos(eps) * Dpsi - V[:,2] * Deps
    Dz =   V[:,0] * np.sin(eps) * Dpsi + V[:,1] * Deps    
    v[:,0] = V[:,0] + Dx
    v[:,1] = V[:,1] + Dy
    v[:,2] = V[:,2] + Dz
    return v


def vnorms(v):
    '''Returns norms of vectors in array v. Faster than numpy.linalg
    '''
    if isinstance(v, np.ndarray) is False:
        print '\n(coords.vnorms) ERROR: Input value not a vector.\n'
        return
    if np.shape(v) == (3,):
        v = v[None, :]
    elif np.shape(v) == (3, 1):
        v = v.T
    N =  np.sqrt(v[:,0]**2 + v[:,1]**2 + v[:,2]**2)
    return N


def vnormd(v):
    '''Returns vector(s) x normalised. Faster than numpy.linalg
    '''
    if isinstance(v, np.ndarray) is False:
        print '\n(coords.vnormd) ERROR: Input value not a vector.\n'
    if np.shape(v) == (3,):
        v = v[None, :]
    elif np.shape(v) == (3, 1):
        v = v.T
    N =  np.sqrt(v[:,0]**2 + v[:,1]**2 + v[:,2]**2)
    N = v / N[:, np.newaxis]
    return N
    

def ttut1(jd=2454195.5, dut1=-0.072073685, dat=33):
    '''TT time to UT1
    '''
    utc = jd
    ut1 = utc + dut1 / 86400.    
    tai = utc + dat / 86400.
    tt = tai + 32.184 / 86400.
    
    #print('utc: {:<18.10f}').format(utc)
    #print('tai: {:<18.10f}').format(tai)
    #print('tt:  {:<18.10f}').format(tt)
    #print('ut1: {:<18.10f}').format(ut1)
    return tt, ut1
    

if __name__ == '__main__':
    # Some tests to compare coordinates.
    import sys
    import matplotlib.pyplot as plt
    import funplot as fp
    #import pysofa as ps
    import novas.compat as nc
    
    if len(sys.argv) == 2:
        npass = sys.argv[-1]
    tjd, epoch, dur, scode = fp.read_function_header(npass)

    T = Transf(ellipsoid='bollocks')
    tfn = np.arange(648 * 60, 940 * 60 + 1, 60)
    JD = 2440000.5 + tjd + tfn / 86400.
    
    Afn, Tfn = fp.azelsat(npass, tfn, rect=True, debug=0, geompos=1)
    Gfn = T.vtop2geo(Tfn) + T.geod2geo() * 1e-6
    Ifn = np.array([geo2gei(v, jd) for v, jd in zip(Gfn, JD)])
    #Ifn_tes = [np.dot(ps.nutm80(jd, 0), n) for jd, n in zip(JD, Ifn)]
    #Ifn_tes = np.array(Ifn_tes).reshape(len(Ifn_tes), 3)

    orbit = np.loadtxt('/home/jose/data/orbits/orbit.tod', skiprows=1)
    orbit = orbit[(orbit[:,0] >= 648) & (orbit[:,0] <= 940)]
    torb = orbit[:, 0]
    Iorb = orbit[:, 1:4] * 1e-6
    Gorb = np.array([gei2geo(v, jd) for v, jd in zip(Iorb, JD)])
    Torb = T.vgeo2top(Gorb - T.geod2geo() * 1e-6)
    Aorb = vtop2azel(Torb)

    cpf = np.loadtxt('/home/jose/code/fortran/beta/interp/fort.2',
                     skiprows=4, usecols=[2,3,4,5,6])
    cpf = cpf[(cpf[:,0] == 55983) & (cpf[:,1] >= 648 * 60) & 
              (cpf[:,1] <= 940 * 60)]
    tcpf = cpf[:,0] + 2400000.5 + cpf[:,1] / 86400.
    Gcpf = cpf[:,2:] * 1e-6
    Tcpf = T.vgeo2top(Gcpf - T.geod2geo()*1e-6)
    Acpf = vtop2azel(Tcpf)    
    Icpf = np.array([geo2gei(v, jd) for v, jd in zip(Gcpf, tcpf)])
    
    bf = np.loadtxt('/home/jose/data/orbits/ETA_55983_bf', skiprows=1,
                    usecols=[0, 1, 2, 3, 4])
    bf = bf[(bf[:,0] == 55983) & (bf[:,1] >= 648) & (bf[:,1] <= 940)]
    Gbf = bf[:,2:] * 1e-6
    Tbf = T.vgeo2top(Gbf - T.geod2geo() * 1e-6)
    Abf = vtop2azel(Tbf)    
    Ibf = np.array([geo2gei(v, jd) for v, jd in zip(Gbf, JD)])
    
    j2k = np.loadtxt('/home/jose/data/orbits/ETA_55983_J2K', skiprows=1,
                     usecols=[0, 1, 2, 3, 4])
    j2k = j2k[(j2k[:,0] == 55983) & (j2k[:,1] >= 648) &
              (j2k[:,1] <= 940)]
    Jj2k = j2k[:,2:]
    Jj2k = np.array([nc.nutation(jd, v) for jd, v in zip(JD, Jj2k)])
    
    utc = JD
    dut1 = -0.4712715
    dat = 34
    
    TT = [ttut1(jd, dut1=dut1, dat=dat)[0] for jd in utc]
    
    #Ij2k = np.array([np.dot(ps.pnm80(tt, 0), v) for tt, v in zip(TT, Jj2k)])
    #Ij2k = Ij2k.reshape(293,3)
    Gj2k= np.array([gei2geo(v, jd) for v, jd in zip(Ij2k, JD)])

    
    f1 = plt.figure(1); f1.clear()
    ax = f1.add_subplot(211)
    ax.plot(1e6 * vnorms(Ifn - Iorb), 'b', label='fn-orb')
    ax.plot(1e6 * vnorms(Ifn - Icpf), 'g', label='fn-cpf')
    ax.plot(1e6 * vnorms(Ifn - Ibf), 'm', label='fn-bf')
    ax.set_ylabel('$V_1-V_2$ inertial (m)'); ax.grid('on')
    ax.legend(loc=1, frameon=False)

    ax2 = f1.add_subplot(212)
    ax2.plot(3600 * 180 / np.pi * (np.mod(Afn[:,0], 2 * np.pi) - Aorb[:,0]),
             label='fn-orb $\Delta$az (arcs)')
    ax2.plot(3600 * 180 / np.pi * (Afn[:,1] - Aorb[:,1]),
             label='fn-orb $\Delta$el (arcs)')
    ax2.plot(1e6 * (Afn[:,2] - Aorb[:,2]), label='$\Delta$r (m)')
    ax2.set_ylabel('fn-orb $V_1-V_2$ topocentric')
    ax2.plot(3600 * 180 / np.pi * (np.mod(Afn[:,0], 2 * np.pi) - Acpf[:,0]),
             label='fn-cpf $\Delta$az')
    ax2.plot(3600 * 180 / np.pi * (np.mod(Afn[:,1], 2 * np.pi) - Acpf[:,1]),
             label='fn-cpf $\Delta$el')
    ax2.plot(1e6 * (Afn[:,2] - Acpf[:,2]), label='$\Delta$r')
    ax2.legend(loc=0, frameon=False)
    ax2.set_xlabel('time (min)')
    ax2.grid('on')
    
    f3 = plt.figure(2); f3.clear()
    ax = f3.add_subplot(111, projection='3d')
    ax.plot(Icpf[:,0], Icpf[:,1], Icpf[:,2], color='red')
    
    #f2 = plt.figure(2); f1.clear()
    #ax = f2.add_subplot(111, projection='3d')
    #ax.scatter(Xgei[:,0], Xgei[:,1], Xgei[:,2], s=1, color='r')
    #ax.scatter(Xgei[-1,0], Xgei[-1,1], Xgei[-1,2], s=20, color='r')
    #ax.scatter(Xgei_te[:,0], Xgei_te[:,1], Xgei_te[:,2], s=1, color='g')
    #ax.scatter(ogei[:,0], ogei[:,1], ogei[:,2], color='b',s=1)
    #ax.scatter(ogei[-1,0], ogei[-1,1], ogei[-1,2], color='b', s=20)
    #fp.plot_limits()
    plt.show()
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
