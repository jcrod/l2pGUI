#!/usr/bin/env python

'''Provides functions to compute Sun position and times of sunrise, 
sunset, twilight and sunlight duration. Based on Astronomical Algorithms,
by Jean Meeus. Implementation from the online NOAA solar calculator:

    http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
'''

import sys
import string
import datetime as dt
import math
import argparse
import jdates as jd

D2R = math.pi / 180
R2D = 180 / math.pi


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
    #gmls = math.mod(280.46646 + JC * (36000.76983 + JC * 0.0003032), 360)
    gmls = 280.46646 + JC * (36000.76983 + JC * 0.0003032) % 360
    gmas = 357.52911 + JC * (35999.05029 - 0.0001537 * JC)
    
    # Eccentriciy of Earth orbit
    ecceo = 0.016708634 - JC * (0.000042037 + 0.0000001267 * JC)
    
    # Sun equation of the centre
    seqctr = (math.sin(gmas * D2R) * (1.914602 - JC * (0.004817 +
             0.000014 * JC)) + math.sin(2 * gmas * D2R) * (0.019993 -
             0.000101 * JC) + math.sin(3 * gmas * D2R) * 0.000289)
             
    # Sun true longitude and true anomaly (deg)
    struel = gmls + seqctr
    struea = gmas + seqctr
    
    # Sun radius vector (AU)
    sradv = (1.000001018 * (1 - ecceo**2)) / (1 + ecceo * 
            math.cos(struea * math.pi / 180))
            
    # Sun apparent longitude (deg)
    sappl = struel - 0.00569 - 0.00478 * math.sin(D2R * 
            (125.04 - 1934.136 * JC))
            
    # Mean and corrected obliquity of the ecliptic (deg)
    moblecl = 23 + (26 + ((21.448 - JC * (46.815 + JC * (0.00059 - JC *
              0.001813)))) / 60) / 60
    oblcorr = moblecl + 0.00256 * math.cos(D2R * (125.04 - 1934.136 * JC))

    # Sun right ascension and declination (deg)
    sra = R2D * math.atan2(math.cos(D2R * oblcorr) * math.sin(D2R * sappl), 
          math.cos(D2R * sappl))      
    sde = R2D * (math.asin(math.sin(D2R * oblcorr) * math.sin(D2R * sappl)))

    # Equation of time (min)
    vary = math.tan(D2R * oblcorr/2)**2
    eqtime = 4*R2D * (vary * math.sin(2*D2R * gmls) - 2*ecceo * 
             math.sin(D2R * gmas) + 4*ecceo * vary * math.sin(D2R * gmas) *
             math.cos(2*D2R * gmls) - 0.5*vary**2 * math.sin(4*D2R * gmls) -
             1.25*ecceo**2 * math.sin(2*D2R * gmas))

    return (sra, sde), eqtime
    
    
def sunradec2azel(sde, eqtime, JD, TZ=0, lon=0.3361, lat=50.8674):
    '''Computes Sun's azimuth and elevation, sunrise, sunset and
    solar noon times and sun light duration.
    
    Method from Astronomical Algorithms, by Jean Meeus. Implementation 
    based on NOAA Solar calculator:
    
    http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
    '''
    # Julian century
    JC = (JD - 2451545) / 36525.
    
    # Time past local midnight
    t = JD - math.floor(JD) - 0.5
    
    # Hour angle sunrise (deg).
    # Sun true elevation: -50' = 0.833 deg. -16' semidiam. -34' diffraction
    HAsunr = R2D * (math.acos(math.cos(D2R * 90.833) / 
             (math.cos(D2R * lat) * math.cos(D2R * sde)) -
             math.tan(D2R * lat) * math.tan(D2R * sde)))
             
    # Civil twilight
    try:
        HAcivt = R2D * (math.acos(math.cos(D2R * 96) / 
             (math.cos(D2R * lat) * math.cos(D2R * sde)) -
             math.tan(D2R * lat) * math.tan(D2R * sde)))
    except ValueError:
        HAcivt = -1
    
    # Nautical twilight
    try:
        HAnaut = R2D * (math.acos(math.cos(D2R * 102) / 
             (math.cos(D2R * lat) * math.cos(D2R * sde)) -
             math.tan(D2R * lat) * math.tan(D2R * sde)))
    except ValueError:
        HAnaut = -1
    
    # Astronomical twilight
    try:
        HAastr = R2D * (math.acos(math.cos(D2R * 108) / 
             (math.cos(D2R * lat) * math.cos(D2R * sde)) -
             math.tan(D2R * lat) * math.tan(D2R * sde)))
    except ValueError:
        HAastr = -1
    
    # Solar noon, sunrise and sunset (LST, local sidereal time)
    Solnoon = (720 - 4*lon - eqtime + TZ*60) / 1440

    keys = ['Srise', 'Sset', 'Ctw0', 'Ctw1', 'Ntw0', 'Ntw1', 'Atw0', 'Atw1']
    Stimes = {k: 0 for k in keys}
    Stimes['Solnoon'] = Solnoon
    
    Stimes['Srise'] = Solnoon - HAsunr * 4 / 1440
    Stimes['Sset'] = Solnoon + HAsunr * 4 / 1440
    Stimes['Ctw0'] = Solnoon - HAcivt * 4 / 1440
    Stimes['Ctw1'] = Solnoon + HAcivt * 4 / 1440
    Stimes['Ntw0'] = Solnoon - HAnaut * 4 / 1440
    Stimes['Ntw1'] = Solnoon + HAnaut * 4 / 1440
    Stimes['Atw0'] = Solnoon - HAastr * 4 / 1440
    Stimes['Atw1'] = Solnoon + HAastr * 4 / 1440
    
    for k in ['Srise', 'Ctw0', 'Ntw0', 'Atw0']:
        if Stimes[k] < 0: Stimes[k] += 1

    # Sunlight duration (minutes)
    Stimes['Sunldur'] = 8 * HAsunr
    
    # True solar time (min)
    trueSolt = 1440*t + eqtime + 4*lon - 60*TZ % 1440

    # Hour angle (deg)
    if trueSolt / 4 < 0:
        HA = trueSolt / 4 + 180
    else:
        HA = trueSolt / 4 - 180

    # Zenith angle, elevation and azimuth (deg)
    za = R2D * (math.acos(math.sin(D2R * lat) * 
         math.sin(D2R * sde) + math.cos(D2R * lat) *
         math.cos(D2R * sde) * math.cos(D2R * HA)))
    elsun = 90 - za

    if HA > 0:
        azsun = (R2D* (math.acos(((math.sin(D2R * lat) *
                math.cos(D2R * za)) - math.sin(D2R * sde)) /
                (math.cos(D2R * lat) * math.sin(D2R * za)))) + 180) % 360
    else: 
        azsun = (540 - R2D * (math.acos(((math.sin(D2R * lat) *
                math.cos(D2R * za)) - math.sin(D2R * sde)) /
                (math.cos(D2R * lat) * math.sin(D2R * za))))) % 360

    return (azsun, elsun), Stimes


def refr_corr(elsun):
    '''Refraction correction
    '''
    h = elsun * D2R
    hh = elsun
    if elsun < -0.575:
        c = -20.774  / (3600 * math.tan(h))
    elif (elsun < 5):
        c =(1735 + hh*(-518.2 + hh*(103.4 + hh*(-12.79 + hh*0.711)))) / 3600
    elif (elsun < 85):
        c = (58.1 / math.tan(h) - 0.07 / math.tan(h)**3 + 
             0.000086 / math.tan(h)**5) / 3600
    else:
        c = 0

    return c    
    

def sunpos(JD, lon=0.3361, lat=50.8674, TZ=0):
    '''Returns Sun's azimuth and elevation at specified epoch and location
    '''
    sun_radec, eq_time = sunradeceqtime(JD)
    (azsun, elsun), suntimes = sunradec2azel(sun_radec[1], eq_time, JD,
                                             lon=lon, lat=lat, TZ=TZ)
    elsun = elsun + refr_corr(elsun)
    return azsun, elsun


def sunpos_and_times(JD, lon=0.3361, lat=50.8674, TZ=0):
    '''Returns Sun's azimuth and elevation at specified epoch and location 
    and sunrise, sunset and twilight times.
    '''
    sun_radec, eq_time = sunradeceqtime(JD)
    (azsun, elsun), suntimes = sunradec2azel(sun_radec[1], eq_time, JD,
                                             lon=lon, lat=lat, TZ=TZ)
    elsun = elsun + refr_corr(elsun)
    return (azsun, elsun), suntimes


def pretty_print(name, ih, im, iss):
    line_str = '   {:22s}: {:0=2.0f}:{:0=2.0f}:{:0=2.0f}'
    print(line_str).format(name, ih, im, iss)
    return


def print_time(t, string):
    if (t > 0) and (t < 1):
        pretty_print('{} twilight'.format(string), *jd.s2hms(86400 * t))
    else:
        print('   No {} twilight'.format(string))
        

def parseDates(date=None, mode='mjd'):
    if mode == 'greg':
        try:
            d = dt.datetime.strptime(date, '%Y-%m-%d %H:%M:%S')
        except ValueError, error:
            print('Could not parse date: {}'.format(error))
            sys.exit()
        fd = jd.hms2s(d.hour, d.minute, d.second) / 86400.
        JD = jd.gcal2jd(d.year, d.month, d.day) + fd
    elif mode == 'mjd':
        iy, im, iday, fd = jd.jd2gcal(2400000.5 + date)
        ih, imin, iss = jd.s2hms(fd * 86400)
        d = dt.datetime(*[int(n) for n in (iy, im, iday, ih, imin, iss)])
        JD = 2400000.5 + date
    elif mode == 'now':
        d = dt.datetime.utcnow()
        fd = jd.hms2s(d.hour, d.minute, d.second) / 86400.    
        JD = jd.gcal2jd(d.year, d.month, d.day) + fd
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
                        help='latitude of observing station', type=float)
    parser.add_argument('-tz', '--time-zone', type=int, default=0,
                        help='Time zone')
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
    TZ = args.time_zone
    
    (azsun, elsun), stimes = sunpos_and_times(JD, lon=Lon, lat=Lat, TZ=TZ)
    
    print(d.strftime('\n   ===== %Y-%m-%d %H:%M:%S =====\n'))
    print('   Azimuth/Elevation: {:6.2f} {:6.2f}').format(azsun, elsun)
    pretty_print('Sunlight duration', *jd.s2hms(stimes['Sunldur'] * 60))
    print '   ' + '-' * 33
    print_time(stimes['Atw0'], 'Astronomical')
    print_time(stimes['Ntw0'], 'Nautical')
    print_time(stimes['Ctw0'], 'Civil')
    print '   ' + '-' * 33
    pretty_print('Sunrise', *jd.s2hms(86400 * stimes['Srise']))
    pretty_print('Solar noon', *jd.s2hms(86400 * stimes['Solnoon']))
    pretty_print('Sunset', *jd.s2hms(86400 * stimes['Sset']))
    print '   ' + '-' * 33
    print_time(stimes['Ctw1'], 'Civil')
    print_time(stimes['Ntw1'], 'Nautical')
    print_time(stimes['Atw1'], 'Astronomical')
    print '\n'


if __name__ == '__main__':
    main()
 
    
