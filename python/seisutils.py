import numpy as np
from scipy.signal import butter, lfilter
from scipy import signal

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


def transfer(data, delta, freq_limits, units, pz_file):
    """
    Remove the instrument response from a seismogram using a SAC Pole Zero file.

    :param data: Vector containing seismic data
    :param delta: sample rate of the seismic data (s)
    :param freq_limits: Vector containing corner frequencies of cosine filter applied to seismic data
                        before deconvolution. This step is necessary to stabilize the division.
    :param units: Desired output units. Currently supported values are 'displacement', 'velocity', or 'acceleration'
    :param pz_file: Full path to the SAC Pole Zero file

    :return: data: Vector containing the seismic data with the instrument response removed
    """

    # Need to add some input validation checks to this function!

    # Maybe make two sub-functions - one to parse the pole zero file, and one to apply the cosine filter?

    # Parse the Pole Zero file
    z_string = 'ZEROS'
    p_string = 'POLES'
    c_string = 'CONSTANT'

    with open(pz_file) as pz:
        for num, line in enumerate(pz, 1):
            if z_string in line:
                z_line = num
            if p_string in line:
                p_line = num
            if c_string in line:
                c_line = num

    z = np.loadtxt(pz_file, skiprows=z_line, max_rows=(p_line - z_line - 1))
    z = z.view(np.complex)
    p = np.loadtxt(pz_file, skiprows=p_line, max_rows=(c_line - p_line - 1))
    p = p.view(np.complex)
    k = np.loadtxt(pz_file, skiprows=(c_line - 1), usecols=1)

    # FFT parameters
    npts = len(data)
    nextpow2 = np.ceil(np.log2(npts))
    nfft = int(2 ** nextpow2)
    dfreq = 1 / (nfft * delta)
    nfreq = int((nfft / 2) + 1)

    # Construct frequency grid
    f = np.linspace(0, 20, nfreq)
    f = 2 * np.pi * f

    # Convert arrays to lists
    z = list(z.flatten())
    p = list(p.flatten())

    # Construct the transfer function from the poles, zeros, and constant
    t_function = signal.ZerosPolesGain(z, p, k)
    w, h = signal.freqresp(t_function, f)

    # Invert the transfer function spectrum
    iresp = np.zeros(nfreq, dtype=complex)
    for i in np.arange(1, nfreq):
        denr = (h[i].real ** 2 + h[i].imag ** 2)
        denr = 1.0 / denr
        if denr < 1e-37:
            iresp[i] = complex(0, 0)
        else:
            iresp[i] = complex(h[i].real * denr, -h[i].imag * denr)

    # Filter the data (maybe put this in a separate function..?)
    f1 = freq_limits[0]
    f2 = freq_limits[1]
    f3 = freq_limits[2]
    f4 = freq_limits[3]
    for i in np.arange(1, nfreq):
        freq = i * dfreq
        if freq < f1:
            fac = 0.0
        elif f1 <= freq <= f2:
            fac = 0.5 * (1 - np.cos(np.pi * (freq - f1) / (f2 - f1)))
        elif f3 <= freq <= f4:
            fac = 0.5 * (1 + np.cos(np.pi * (freq - f3) / (f4 - f3)))
        elif freq > f4:
            fac = 0.0
        else:
            fac = 1.0

        iresp[i] = iresp[i] * fac

    # Take the Fourier transform of the data
    data_fft = np.fft.fft(data, nfft)

    # Multiply transformed data by the transfer function
    for i in np.arange(1, nfreq):
        tempR = iresp[i].real * data_fft[i].real - iresp[i].imag * data_fft[i].imag
        tempI = iresp[i].real * data_fft[i].imag + iresp[i].imag * data_fft[i].real
        data_fft[i] = complex(tempR, tempI)

        if i < (nfreq - 1):
            j = nfft - i
            data_fft[j] = complex(tempR, -tempI)

    data_fft[0] = complex(0, 0)
    data_fft[-1] = complex(np.sqrt(data_fft[-1].real * data_fft[-1].real + data_fft[-1].imag * data_fft[-1].imag), 0)

    data = np.fft.ifft(data_fft, nfft)[0:npts].real

    return data
