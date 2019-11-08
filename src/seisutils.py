import numpy as np
from scipy.signal import butter, lfilter

# This file contains various functions which can be used for routine seismic data processing.
# ---------------------------------------------------------------------------------------------
# Last updated 10/22/2019 by aburky@princeton.edu
# ---------------------------------------------------------------------------------------------


def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance (degrees) between a pair of latitude/longitude points, as well as the
    azimuth from the first coordinate to the second.

    :param lat1: Latitude of first coordinate (degrees)
    :param lon1: Longitude of first coordinate (degrees)
    :param lat2: Latitude of second coordinate (degrees)
    :param lon2: Longitude of second coordinate (degrees)

    :return: Two floats: the great circle distance (degrees) and azimuth (degrees)
    """
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
    if lon1 < 0:
        lon1 = np.deg2rad(lon1+360)
    else:
        lon1 = np.deg2rad(lon1)
    if lon2 < 0:
        lon2 = np.deg2rad(lon2+360)
    else:
        lon2 = np.deg2rad(lon2)
    # Radius of Earth (km)
    r = 6371
    # Kilometers in one degree of lat/lon distance
    d2km = 111.1949
    # Calculate great circle distance
    gcarc = 2*r*np.arcsin(np.sqrt(np.sin((lat2-lat1)/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2)**2))
    delta = gcarc/d2km
    # Calculate azimuth from point 1 to point 2
    azimuth = np.arctan2(np.sin(lon2-lon1)*np.cos(lat2), np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)
                         * np.cos(lon2-lon1))
    azimuth = np.rad2deg(azimuth)
    if azimuth < 0:
        azimuth = azimuth + 360
    return delta, azimuth


def seisrt(ndat, edat, baz):
    """
    Rotate seismic data from a North-East coordinate system into the ray-based, radial-transverse coordinate system.

    :param ndat: North component of seismic data
    :param edat: East component of seismic data
    :param baz: Back-azimuth (degrees) from the seismic station to the seismic source, measured clockwise from North

    :return: rdat, tdat: Arrays containing the rotated radial and transverse seismic data
    """
    # Add 180 degrees to back-azimuth to properly orient to raypath
    if baz < 180:
        baz = baz + 180
    elif baz >= 180:
        baz = baz - 180
    # Construct rotation matrix
    m2d = np.array([[np.cos(baz*(np.pi/180)), np.sin(baz*(np.pi/180))],
                   [-np.sin(baz*(np.pi/180)), np.cos(baz*(np.pi/180))]])
    dvc = np.array([ndat, edat])
    rot = np.matmul(m2d, dvc)
    rdat = rot[0]
    tdat = rot[1]
    return rdat, tdat


def seisne(dat1, dat2, azimuth):
    """
    Rotate horizontal channels of seismic data (BH1, BH2) to North and East orientations (BHN, BHE) using
    the azimuth of channel 'BH1'

    :param dat1: Channel 1 of seismic data (BH1, HH1, etc.)
    :param dat2: Channel 2 of seismic data (BH2, HH2, etc.)
    :param azimuth: Azimuth (degrees) of first channel of seismic data

    :return: ndat, edat: Arrays containing properly oriented North and East components of seismic data
    """
    # Construct the rotation matrix
    m2d = np.array([[np.cos(azimuth*(np.pi/180)), -np.sin(azimuth*(np.pi/180))],
                    [np.sin(azimuth*(np.pi/180)), np.cos(azimuth*(np.pi/180))]])
    dvc = np.array([dat1, dat2])
    rot = np.matmul(m2d, dvc)
    ndat = rot[0]
    edat = rot[1]
    return ndat, edat


def bp_butter(data, lowcut, highcut, fs, order=5):
    """
    Bandpass filter an array using a Butterworth filter, default order is 5.

    :param data: Array data to be filtered
    :param lowcut: Low corner frequency (Hz)
    :param highcut: High corner frequency (Hz)
    :param fs: Sample rate of data
    :param order: Order of filter. Default is 5

    :return: y: Array containing filtered input data
    """
    nyq = 0.5*fs
    low = lowcut/nyq
    high = highcut/nyq
    b, a = butter(order, [low, high], btype='band')
    y = lfilter(b, a, data, axis=0)
    return y