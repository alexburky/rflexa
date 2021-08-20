import numpy as np

# This file contains functions relevant to the removal of the instrument response from seismic data

# ---------------------------------------------------------
# Last updated 8/20/2021 by aburky@princeton.edu
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

    # FFT parameters
    npts = len(data)
    nextpow2 = np.ceil(np.log2(npts))
    nfft = int(2 ** nextpow2)
    dfreq = 1 / (nfft * delta)
    nfreq = int((nfft / 2) + 1)

    # Initialize vectors
    iresp = np.zeros(nfreq, dtype=complex)
    denr = np.zeros(nfreq)

    return data
