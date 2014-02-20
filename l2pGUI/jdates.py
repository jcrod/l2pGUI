#!/usr/bin/env python
'''Various functions to deal with Julian Dates'''

import numpy as np
import datetime as dt

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
    

def jdNow():
    '''JD right now'''
    d = dt.datetime.utcnow()
    fd = ((d.hour * 60 + d.minute) * 60 + d.second) / 86400.
    return gcal2jd(d.year, d.month, d.day + fd)


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