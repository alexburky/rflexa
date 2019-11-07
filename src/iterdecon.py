import numpy as np
from scipy import fftpack

# This file contains all of the functions necessary to perform iterative time-domain deconvolution
# as outlined in Ligorria & Ammon 1999.
# --------------------------------------------------------------------------------------------------
# Last updated 10/22/2019 by aburky@princeton.edu
# --------------------------------------------------------------------------------------------------


def gauss_filter(dt, nft, f0):
    """
    Construct a gaussian filter of a prescribed width in the frequency-domain.

    :param dt: sampling interval (seconds) -> typically the sampling interval of the seismic data
    :param nft: length of gaussian filter (samples)
    :param f0: gaussian width factor (the constant in the denominator of the exponential function)

    :return: gauss: Array containing the desired gaussian filter in the frequency-domain
    """
    df = 1.0/(nft * dt)
    nft21 = int(0.5*nft + 1)
    # Construct x-axis (frequency) vectors
    f = df*np.arange(0, nft21, 1)
    w = 2*np.pi*f
    # Construct gaussian filter
    gauss = np.zeros(nft)
    gauss[0:nft21] = np.exp(-0.25*(w/f0)**2)/dt
    gauss[nft21:] = np.flip(gauss[1:nft21-1])
    return gauss


def gauss_convolve(a, b, nf, ndt):
    """
    Convolve time-series data with a gaussian filter by Fourier transforming into the frequency-domain, multiplying,
    and inverse-Fourier transforming back into the time-domain.

    :param a: Array of time series data
    :param b: Gaussian filter (output of gauss_filter() function)
    :param nf: Number of points (samples) in Fourier transform
    :param ndt: Sampling interval of data (seconds)

    :return: ab: Array containing the convolution of the data and the gaussian filter
    """
    afr = fftpack.fft(a, nf)
    bfr = (afr * b) * ndt
    ab = np.real(np.fft.ifft(bfr, nf))
    return ab


def correlate(R, W, nfft):
    """
    Calculate the cross-correlation between two time-series by Fourier transforming and multiplying by the complex
    conjugate of one of the time-series.

    :param R: Array of time series data (radial component for Z -> R receiver functions)
    :param W: Array of time series data (vertical component for Z -> R receiver functions)
    :param nfft: Number of points (samples) in Fourier transform

    :return: x: Array containing the resulting cross-correlation
    """
    x = fftpack.ifft(fftpack.fft(R, nfft) * np.conj(fftpack.fft(W, nfft)), nfft).real
    return x


def phase_shift(x, nfft, dt, tshift):
    """
    Shift a vector of time-series data by a given number of seconds.

    :param x: Input vector of time-series data
    :param nfft: Number of points (samples) in Fourier transform
    :param dt: Sampling interval of time-series data (seconds)
    :param tshift: Desired time shift (seconds)

    :return: x: Time shifted version of input vector
    """
    # Go into the frequency domain
    xf = fftpack.fft(x, nfft)
    # Phase shift in radians
    shift_i = round(tshift / dt)
    p = 2 * np.pi * np.arange(1, nfft + 1) * shift_i / (nfft)
    # Apply shift
    xf = xf * (np.cos(p) - 1j * np.sin(p))
    # Back into time
    x = fftpack.ifft(xf, nfft).real / np.cos(2 * np.pi * shift_i / nfft)
    return x


def next_power_2(x):
    """
    Determine the next integer that is 2 raised to some power.

    :param x: Number which you would like to find the next power of 2 for

    :return: x: Number which is 2 raised to some power
    """
    # Function which finds the nearest number that is 2 raised to some power
    return 1 if x == 0 else 2**(x-1).bit_length()


def iterdecon(num, den, dt, nt, tshift, f0, itmax, errtol):
    """
    Calculate a receiver function using the iterative time-domain deconvolution algorithm outlined by
    Ligorria & Ammon 1999.

    :param num: Numerator in deconvolution (radial component data for Z -> R receiver functions)
    :param den: Denominator in deconvolution (vertical component data for Z -> R receiver functions)
    :param dt: Sampling interval of data (seconds)
    :param nt: Length of input data vectors (samples)
    :param tshift: Desired time shift of resulting receiver function (seconds)
    :param f0: Gaussian width factor determining width of gaussian filter used in deconvolution
    :param itmax: Maximum allowed number of iterations before giving up and outputting receiver function
    :param errtol: Minimum change in error between iterations before giving up and outputting receiver function

    :return: RFI, RMS: An array containing the resulting receiver function (time-domain) and an array containing the
                       RMS error from each iteration.
    """
    # Initiate iterative deconvolution
    rms = np.zeros(itmax)
    nfft = next_power_2(nt)
    p0 = np.zeros(nfft)

    # Resize and rename numerator and denominator
    u0 = np.zeros(nfft)
    w0 = u0.copy()
    u0[0:nt] = num
    w0[0:nt] = den

    # Construct Gaussian Filter
    gF = gauss_filter(dt, nfft, f0)

    # Apply Gaussian Filter to signals
    u = gauss_convolve(u0, gF, nfft, dt)
    w = gauss_convolve(w0, gF, nfft, dt)
    # Get Vertical Component in Frequency Domain
    wf = fftpack.fft(w0, nfft)
    r = u.copy()

    # Get power in numerator for error scaling
    powerU = np.sum(u**2)

    # Loop through the iterations
    it = -1
    sumsq_i = 1
    d_error = 100*powerU + errtol
    maxlag = int(0.5*nfft)

    while abs(d_error) > errtol and it < itmax:
        # Advance iteration
        it = it + 1

        # Cross - correlate signals
        rw = correlate(r, w, nfft)
        rw = rw/np.sum(w**2)
        # Get index of maximum of cross correlation
        i1 = np.argmax(abs(rw[0:maxlag]))
        amp = rw[i1]/dt

        # Compute predicted deconvolution
        p0[i1] += amp
        p = gauss_convolve(p0, gF, nfft, dt)
        p = gauss_convolve(p, wf, nfft, dt)

        # Compute residual with filtered numerator
        r = u - p.copy()
        sumsq = np.sum(r**2)/powerU
        rms[it] = sumsq
        d_error = 100*(sumsq_i - sumsq)
        sumsq_i = sumsq

    # Compute final receiver function
    p = gauss_convolve(p0, gF, nfft, dt)

    # Apply optional time shift to receiver function
    p = phase_shift(p, nfft, dt, tshift)
    # Output first nt samples of final receiver function
    rfi = p[0:nt]
    # Output RMS values
    rms = rms[0:it]

    return rfi, rms
