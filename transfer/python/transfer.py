import numpy as np
from scipy import signal

# This file contains functions relevant to the removal of the instrument response from seismic data

# ---------------------------------------------------------
# Last updated 8/25/2021 by aburky@princeton.edu
# ---------------------------------------------------------


def transfer(data, delta, freq_limits, units, pz_file):
    """
    Remove the instrument response from a seismogram using a SAC Pole Zero file.

    :param data: Vector containing seismic data
    :param delta: sample rate of the seismic data (s)
    :param freq_limits: Vector containing corner frequencies of cosine filter applied to seismic data
                        before deconvolution. This step is necessary to stabilize the division.
    :param units: Desired output units. Currently supported values are 'displacement', 'velocity', or
                  'acceleration'
    :param pz_file: Full path to the SAC Pole Zero file

    :return data: Vector containing the seismic data with the instrument response removed
    """

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

    zeros = np.loadtxt(pz_file, skiprows=z_line, max_rows=(p_line - z_line - 1))
    zeros = zeros.view(np.complex)
    p = np.loadtxt(pz_file, skiprows=p_line, max_rows=(c_line - p_line - 1))
    p = p.view(np.complex)
    k = np.loadtxt(pz_file, skiprows=(c_line - 1), usecols=1)

    # Construct zeros array based on desired output units
    nz = np.nonzero(zeros)
    if units == 'displacement':
        z = np.array([[0j], [0j], [0j]])
    elif units == 'velocity':
        z = np.array([[0j], [0j]])
    elif units == 'acceleration':
        z = np.array([[0j]])
    else:
        print('Error: Currently accepted output units are displacement, velocity, or acceleration')

    if len(nz[0]) != 0:
        for i in range(0, len(nz[0])):
            z = np.append(z, zeros[nz[0][i]])

    # FFT parameters
    npts = len(data)
    nextpow2 = np.ceil(np.log2(npts))
    nfft = int(2 ** nextpow2)
    dfreq = 1 / (nfft * delta)
    nfreq = int((nfft / 2) + 1)
    nyq = (1/delta)*0.5

    # Pre-initialize vectors?
    iresp = np.zeros(nfreq, dtype=complex)
    denr = np.zeros(nfreq)

    # Construct frequency grid
    f = np.linspace(0, nyq, nfreq)
    omega = 2 * np.pi * f

    # Convert arrays to lists
    z = list(z.flatten())
    p = list(p.flatten())

    # Construct the transfer function from the poles, zeros, and constant
    t_function = signal.ZerosPolesGain(z, p, k)
    w, h = signal.freqresp(t_function, omega)

    # Invert the transfer function spectrum
    denr[1:] = 1.0 / (pow(h[1:].real, 2) + pow(h[1:].imag, 2))
    iresp[1:] = np.vectorize(complex)(np.multiply(h[1:].real, denr[1:]), np.multiply(-h[1:].imag, denr[1:]))

    # Get frequency limits of cosine filter
    f1 = freq_limits[0]
    f2 = freq_limits[1]
    f3 = freq_limits[2]
    f4 = freq_limits[3]

    # Move the construction of the cosine filter to a different function!
    cosine_filter = np.zeros(nfreq, dtype=float)
    idx = np.where(f < f1)
    cosine_filter[idx] = 0.0
    idx = np.where(np.logical_and(f >= f1, f <= f2))
    cosine_filter[idx] = 0.5 * (1 - np.cos(np.pi * (f[idx] - f1) / (f2 - f1)))
    idx = np.where(np.logical_and(f > f2, f < f3))
    cosine_filter[idx] = 1.0
    idx = np.where(np.logical_and(f >= f3, f <= f4))
    cosine_filter[idx] = 0.5 * (1 + np.cos(np.pi * (f[idx] - f3) / (f4 - f3)))
    idx = np.where(f > f4)
    cosine_filter[idx] = 0.0

    # Multiply the inverted response by the cosine filter
    iresp = iresp * cosine_filter

    # Take the Fourier transform of the data
    data_fft = np.fft.fft(data, nfft)

    # Multiply the data by the inverted, filtered instrument response
    data_fft[0:nfreq] = data_fft[0:nfreq] * iresp
    data_fft[nfreq:] = np.conj(np.flipud(data_fft[1:(nfreq - 1)]))

    data_fft[0] = complex(0, 0)
    data_fft[-1] = complex(np.sqrt(data_fft[-1].real * data_fft[-1].real + data_fft[-1].imag * data_fft[-1].imag), 0)

    # Inverse transform back to the time domain
    data = np.fft.ifft(data_fft, nfft)[0:npts].real

    return data

    return data
