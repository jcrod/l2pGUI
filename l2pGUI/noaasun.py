#!/usr/bin/env python

'''Provides functions to compute Sun position and times of sunrise, 
sunset, twilight and sunlight duration. Based on Astronomical Algorithms,
by Jean Meeus. Implementation from the online NOAA solar calculator:

    http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
'''

import sys
import string
import datetime as dt
import numpy as np
import argparse
import coords as cds

D2R = np.pi / 180
R2D = 180 / np.pi


def sunradeceqtime(JD):
    '''Computes Sun's right ascension and declination and
    the equation of time.
    
    Method from Astronomical Algorithms, by Jean Meeus. Implementation 
    based on NOAA Solar calculator:
    
    http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
    '''
    # Julian century
    JC = (JD - 2451545) / 36525.
    
    # Geom mean longitude and mean anomaly of the Sun (deg)
    gmls = np.mod(280.46646 + JC * (36000.76983 + JC * 0.0003032), 360)
    gmas = 357.52911 + JC * (35999.05029 - 0.0001537 * JC)
    
    # Eccentriciy of Earth orbit
    ecceo = 0.016708634 - JC * (0.000042037 + 0.0000001267 * JC)
    
    # Sun equation of the centre
    seqctr = (np.sin(gmas * D2R) * (1.914602 - JC * (0.004817 +
             0.000014 * JC)) + np.sin(2 * gmas * D2R) * (0.019993 -
             0.000101 * JC) + np.sin(3 * gmas * D2R) * 0.000289)
             
    # Sun true longitude and true anomaly (deg)
    struel = gmls + seqctr
    struea = gmas + seqctr
    
    # Sun radius vector (AU)
    sradv = (1.000001018 * (1 - ecceo**2)) / (1 + ecceo * 
            np.cos(struea * np.pi / 180))
            
    # Sun apparent longitude (deg)
    sappl = struel - 0.00569 - 0.00478 * np.sin(D2R * 
            (125.04 - 1934.136 * JC))
            
    # Mean and corrected obliquity of the ecliptic (deg)
    moblecl = 23 + (26 + ((21.448 - JC * (46.815 + JC * (0.00059 - JC *
              0.001813)))) / 60) / 60
    oblcorr = moblecl + 0.00256 * np.cos(D2R * (125.04 - 1934.136 * JC))

    # Sun right ascension and declination (deg)
    sra = R2D * np.arctan2(np.cos(D2R * oblcorr) * np.sin(D2R * sappl), 
          np.cos(D2R * sappl))      
    sde = R2D * (np.arcsin(np.sin(D2R * oblcorr) * np.sin(D2R * sappl)))

    # Equation of time (min)
    vary = np.tan(D2R * oblcorr/2)**2
    eqtime = 4*R2D * (vary * np.sin(2*D2R * gmls) - 2*ecceo * 
             np.sin(D2R * gmas) + 4*ecceo * vary * np.sin(D2R * gmas) *
             np.cos(2*D2R * gmls) - 0.5*vary**2 * np.sin(4*D2R * gmls) -
             1.25*ecceo**2 * np.sin(2*D2R * gmas))

    return (sra, sde), eqtime
    
    
def sunradec2azel(sde, eqtime, JD, TZ=0, Lon=0.3361, Lat=50.8674):
    '''Computes Sun's azimuth and elevation, sunrise, sunset and
    solar noon times and sun light duration.
    
    Method from Astronomical Algorithms, by Jean Meeus. Implementation 
    based on NOAA Solar calculator:
    
    http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
    '''
    # Julian century
    JC = (JD - 2451545) / 36525.
    
    # Time past local midnight
    t = JD - np.floor(JD) - 0.5
    
    # Hour angle sunrise (deg).
    # Sun true elevation: -50' = 0.833 deg. -16' semidiam. -34' diffraction
    HAsunr = R2D * (np.arccos(np.cos(D2R * 90.833) / 
             (np.cos(D2R * Lat) * np.cos(D2R * sde)) -
             np.tan(D2R * Lat) * np.tan(D2R * sde)))
             
    # Civil twilight
    HAcivt = R2D * (np.arccos(np.cos(D2R * 96) / 
             (np.cos(D2R * Lat) * np.cos(D2R * sde)) -
             np.tan(D2R * Lat) * np.tan(D2R * sde)))
    
    # Nautical twilight
    HAnaut = R2D * (np.arccos(np.cos(D2R * 102) / 
             (np.cos(D2R * Lat) * np.cos(D2R * sde)) -
             np.tan(D2R * Lat) * np.tan(D2R * sde)))    
    
    # Astronomical twilight
    HAastr = R2D * (np.arccos(np.cos(D2R * 108) / 
             (np.cos(D2R * Lat) * np.cos(D2R * sde)) -
             np.tan(D2R * Lat) * np.tan(D2R * sde)))    
    
    # Solar noon, sunrise and sunset (LST, local sidereal time)
    Solnoon = (720 - 4*Lon - eqtime + TZ*60) / 1440
    Sunrise = Solnoon - HAsunr * 4 / 1440
    Sunset = Solnoon + HAsunr * 4 / 1440
    Civtw0 = Solnoon - HAcivt * 4 / 1440
    Civtw1 = Solnoon + HAcivt * 4 / 1440
    Nautw0 = Solnoon - HAnaut * 4 / 1440
    Nautw1 = Solnoon + HAnaut * 4 / 1440
    Asttw0 = Solnoon - HAastr * 4 / 1440
    Asttw1 = Solnoon + HAastr * 4 / 1440
    
    
    # Sunlight duration (minutes)
    Sunldur = 8*HAsunr
    
    # True solar time (min)
    trueSolt = np.mod(1440*t + eqtime + 4*Lon - 60*TZ, 1440)

    # Hour angle (deg)
    if trueSolt / 4 < 0:
        HA = trueSolt / 4 + 180
    else:
        HA = trueSolt / 4 - 180

    # Zenith angle, elevation and azimuth (deg)
    za = R2D * (np.arccos(np.sin(D2R * Lat) * 
         np.sin(D2R * sde) + np.cos(D2R * Lat) *
         np.cos(D2R * sde) * np.cos(D2R * HA)))
    elsun = 90 - za

    if HA > 0:
        azsun = np.mod(R2D* (np.arccos(((np.sin(D2R * Lat) *
                np.cos(D2R * za)) - np.sin(D2R * sde)) /
                (np.cos(D2R * Lat) * np.sin(D2R * za)))) + 180, 360)
    else: 
        azsun = np.mod(540 - R2D * (np.arccos(((np.sin(D2R * Lat) *
                np.cos(D2R * za)) - np.sin(D2R * sde)) /
                (np.cos(D2R * Lat) * np.sin(D2R * za)))), 360)

    return (azsun, elsun), (Sunrise, Solnoon, Sunset, Sunldur,
            Civtw0, Civtw1, Nautw0, Nautw1, Asttw0, Asttw1)


def refr_corr(elsun):
    '''Refraction correction
    '''
    h = elsun * D2R
    hh = elsun
    if elsun < -0.575:
        c = -20.774  / (3600 * np.tan(h))
    elif (elsun < 5):
        c =(1735 + hh*(-518.2 + hh*(103.4 + hh*(-12.79 + hh*0.711)))) / 3600
    elif (elsun < 85):
        c = (58.1/np.tan(h) - 0.07/np.tan(h)**3 + 0.000086/np.tan(h)**5) / 3600
    else:
        c = 0

    return c    
    

def sunpos(JD, Lon=0.3361, Lat=50.8674):
    '''Returns Sun's azimuth and elevation at specified epoch and location
    '''
    sun_radec, eq_time = sunradeceqtime(JD)
    (azsun, elsun), suntimes = sunradec2azel(sun_radec[1], eq_time, JD,
                                             Lon=Lon, Lat=Lat)
    elsun = elsun + refr_corr(elsun)
    return azsun, elsun


def sunpos_and_times(JD, Lon=0.3361, Lat=50.8674):
    '''Returns Sun's azimuth and elevation at specified epoch and location 
    and sunrise, sunset and twilight times.
    '''
    sun_radec, eq_time = sunradeceqtime(JD)
    (azsun, elsun), suntimes = sunradec2azel(sun_radec[1], eq_time, JD,
                                             Lon=Lon, Lat=Lat)
    elsun = elsun + refr_corr(elsun)
    return (azsun, elsun), suntimes


def pretty_print(name, ih, im, iss):
    line_str = '   {:22s}: {:0=2.0f}:{:0=2.0f}:{:0=2.0f}'
    print(line_str).format(name, ih, im, iss)
    return

    
def parseDates(date=None, mode='mjd'):
    if mode == 'greg':
        try:
            d = dt.datetime.strptime(date, '%Y-%m-%d %H:%M:%S')
        except ValueError, error:
            print('Could not parse date: {}'.format(error))
            sys.exit()
        fd = cds.hms2s(d.hour, d.minute, d.second) / 86400.
        JD = cds.gcal2jd(d.year, d.month, d.day) + fd
    elif mode == 'mjd':
        iy, im, iday, fd = cds.jd2gcal(2400000.5 + date)
        ih, imin, iss = cds.s2hms(fd * 86400)
        d = dt.datetime(*[int(n) for n in (iy, im, iday, ih, imin, iss)])
        JD = 2400000.5 + date
    elif mode == 'now':
        d = dt.datetime.utcnow()
        fd = cds.hms2s(d.hour, d.minute, d.second) / 86400.    
        JD = cds.gcal2jd(d.year, d.month, d.day) + fd
    return d, JD


def main(date=None):
    parser = argparse.ArgumentParser(description='Sun position and rise times')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-n', '--now', action='store_true',
                        help='Current date and time')
    group.add_argument('-g', '--gregorian', 
                        help='Gregorian calendar input yyyy-mm-dd hh:mm:ss',
                        nargs='+')
    group.add_argument('-mjd', '--modifiedJulianDate', 
                        help='Modified Julian date input', type=float)
    group.add_argument('-jd', '--julianDate', 
                        help='Julian date input', type=float)
    parser.add_argument('-lon', '--longitude', 
                        help='Longitude of observing station', type=float)
    parser.add_argument('-lat', '--latitude', 
                        help='Latitude of observing station', type=float)
    args = parser.parse_args()
    
    if args.gregorian:
        if len(args.gregorian) == 1:
            args.gregorian.append('00:00:00')
        d, JD = parseDates(' '.join([s for s in args.gregorian]), mode='greg')
    elif args.modifiedJulianDate:
        d, JD = parseDates(args.modifiedJulianDate, mode='mjd')
    elif args.julianDate:
        d, JD = parseDates(args.modifiedJulianDate - 2400000.5, mode='mjd')
    else:
        d, JD = parseDates(mode='now')
    
    iy, im, iday, ih, imin, iss, _, _, _ = d.timetuple()
    
    Lon = 0.3361 if not args.longitude else args.longitude
    Lat = 50.8674 if not args.latitude else args.latitude
    
    (azsun, elsun), suntimes = sunpos_and_times(JD, Lon=Lon, Lat=Lat)
    (sunrise, solnoon, sunset, sunldur, 
     ct0, ct1, nt0, nt1, at0, at1) = suntimes
    print '\n' * 2 + '   ' + '=' * 5,
    date_str = '{:4.0f}-{:0=2.0f}-{:0=2.0f} {:0=2.0f}:{:0=2.0f}:{:0=2.0f}'
    print(date_str).format(iy, im, iday, ih, imin, iss),
    print '=' * 6 + '\n'
    print('   Azimuth/Elevation: {:6.2f} {:6.2f}').format(azsun, elsun)
    pretty_print('Sunlight duration', *cds.s2hms(sunldur * 60))
    print '   ' + '-' * 33
    if not np.isnan(at0):
        pretty_print('Astronomical twilight', *cds.s2hms(86400 * at0))
    else:
        print '   No astronomical twilight'
    pretty_print('Nautical twilight', *cds.s2hms(86400 * nt0))
    pretty_print('Civil twilight', *cds.s2hms(86400 * ct0))
    print '   ' + '-' * 33
    pretty_print('Sunrise', *cds.s2hms(86400 * sunrise))
    pretty_print('Solar noon', *cds.s2hms(86400 * solnoon))
    pretty_print('Sunset', *cds.s2hms(86400 * sunset))
    print '   ' + '-' * 33
    pretty_print('Civil twilight', *cds.s2hms(86400 * ct1))
    pretty_print('Nautical twilight', *cds.s2hms(86400 * nt1))    
    if not np.isnan(at1):
        pretty_print('Astronomical twilight', *cds.s2hms(86400 * at1))
    else:
        print '   No astronomical twilight'
    print '\n'

if __name__ == '__main__':
    main()
 
    
