import obspy
import seisutils as su
import iterdecon as rfunc
import numpy as np
import os
import shutil
from scipy import signal

# ------------------------------------------------------------------------------------------
# computeRFs.py
# This function calculates receiver functions and saves the results to a directory upon
# completion. The user can specify a choice of filter and Gaussian width factors for use
# in the data processing. The resulting SAC files are saved with additional metadata in the
# header:
#
# USER0: Vertical component SNR
# USER1: Radial component SNR
# USER2: RMS error calculated from iterative time domain deconvolution
# USER3: Quality metric, nu
#
# ------------------------------------------------------------------------------------------
# Last updated 8/18/2020 by aburky@princeton.edu
# ------------------------------------------------------------------------------------------


def computeRFs(network, station, location, data_directory, gaussian_width=1.0, high_cut=0, low_cut=0,
               response_removed=True):
    # Define network, station and location
    ntwk = network
    stat = station
    loc = location
    gw = gaussian_width

    show_traces = 0
    # Determine directory containing earthquakes and output directory
    if response_removed:
        sac_dir = data_directory + ntwk + "/" + stat + "/" + loc + "/RFQUAKES/"
        if low_cut != 0 and high_cut != 0:
            rf_dir = data_directory + ntwk + "/" + stat + "/" + loc + "/RFUNCS/FILTERED_" + str(low_cut) + "_" + \
                     str(high_cut) + "/GW" + ''.join(str(gw).split('.')) + "/"
        else:
            rf_dir = data_directory + ntwk + "/" + stat + "/" + loc + "/RFUNCS/UNFILTERED/GW" + \
                     ''.join(str(gw).split('.')) + "/"
    else:
        sac_dir = data_directory + ntwk + "/" + stat + "/" + loc + "/RFQUAKES_COUNTS/"
        if low_cut != 0 and high_cut != 0:
            rf_dir = data_directory + ntwk + "/" + stat + "/" + loc + "/RFUNCS_COUNTS/FILTERED_" + str(low_cut) + \
                     "_" + str(high_cut) + "/GW" + ''.join(str(gw).split('.')) + "/"
        else:
            rf_dir = data_directory + ntwk + "/" + stat + "/" + loc + "/RFUNCS_COUNTS/UNFILTERED/GW" + \
                     ''.join(str(gw).split('.')) + "/"

    if os.path.exists(rf_dir):
        overwrite = input('Receiver function directory exists! Would you like to overwrite it? [y/n]')
        if overwrite == 'y':
            shutil.rmtree(rf_dir)
        elif overwrite == 'n':
            print('Terminating process...')
            quit()
        else:
            print('ERROR: Invalid response. Terminating process...')
            quit()
    if not os.path.exists(rf_dir):
        os.makedirs(rf_dir)

    st = obspy.read(sac_dir + "*.sac")

    # Filter data
    if low_cut != 0 and high_cut != 0:
        for i in range(0, len(st)):
            # For data with no instrument response removed, demean the data first
            if not response_removed:
                st[i].detrend('demean')
            # Test with rounded vs. un-rounded data
            fs = 1/st[i].meta.sac.delta
            lo = low_cut
            hi = high_cut
            st[i].data = su.bp_butter(st[i].data, lo, hi, fs, 3)

    # Rotate data
    for i in range(0, len(st)):
        if (st[i].stats.sac.kcmpnm == u'BH1') or (st[i].stats.sac.kcmpnm == u'BHN'):
            tb = st[i].meta.starttime
            te = st[i].meta.endtime
            for j in range(-2, 3):
                if (st[(i+j)].stats.sac.kcmpnm == u'BH2') or (st[(i+j)].stats.sac.kcmpnm == u'BHE') and \
                   (tb - 0.5 <= st[(i+j)].meta.starttime <= tb + 0.5 or te - 0.5 <= st[(i+j)].meta.endtime <= te + 0.5):
                    # print('Rotating horizontals...', st[i].id, st[(i+j)].id)
                    ch1 = i
                    ch2 = int(i + j)
                    # Check if length of data vectors is equal
                    if st[ch1].meta.npts != st[ch2].meta.npts:
                        npmin = np.min([st[ch1].meta.npts, st[ch2].meta.npts])
                        st[ch1].data = st[ch1].data[0:npmin]
                        st[ch2].data = st[ch2].data[0:npmin]
                    st[ch1].data, st[ch2].data = su.seisne(st[ch1].data.astype(float), st[ch2].data.astype(float),
                                                           st[ch1].meta.sac.cmpaz)
                    st[ch1].data, st[ch2].data = su.seisrt(st[ch1].data.astype(float), st[ch2].data.astype(float),
                                                           st[ch1].meta.sac.baz)
                    break
                elif j == 2:
                    print('Couldn''t find matching horizontal data for event:', tb)

    # Cut and taper data
    for i in range(0, len(st)):
        print(i)
        # Get indices of window around P arrival
        window_start = 30
        window_end = 90
        bidx = int(round((st[i].meta.sac.t0 - window_start) * st[i].meta.sampling_rate) - 1)
        eidx = int(round((st[i].meta.sac.t0 + window_end) * st[i].meta.sampling_rate))
        if eidx > st[i].data.size:
            print('Array size is problematic, moving on to the next trace...')
            continue
        else:
            st[i].data = st[i].data[bidx:eidx]
            # Taper data
            window = signal.tukey(eidx - bidx, alpha=0.25)
            st[i].data = st[i].data*window
            print('Taper successful')

        if show_traces == 1:
            plt.ion()
            plt.plot(st[i].data, 'k', linewidth=0.25)
            plt.title(st[i].meta.starttime)
            plt.show(block=False)
            plt.pause(0.05)
            plt.cla()

    # Set deconvolution parameters
    tshift = 10
    itmax = 1000
    tol = 0.001
    # Deconvolve!
    rf = []
    rms = []
    rad_idx = []
    snr_z = []
    snr_r = []
    k = 0
    for i in range(0, len(st)):
        if (st[i].stats.sac.kcmpnm == u'BH1') or (st[i].stats.sac.kcmpnm == u'BHE'):
            tb = st[i].meta.starttime
            te = st[i].meta.endtime
            for j in range(-2, 3):
                if (st[(i+j)].stats.sac.kcmpnm == u'BHZ') and \
                   (tb - 0.5 <= st[(i+j)].meta.starttime <= tb + 0.5 or te - 0.5 <= st[(i+j)].meta.endtime <= te + 0.5) and\
                   (st[i].data.size == st[(i+j)].data.size):
                    chz = int(i+j)
                    rf.append(k)
                    rms.append(k)
                    rad_idx.append(k)
                    snr_z.append(k)
                    snr_r.append(k)
                    # Get indices for calculating signal to noise ratio (method of Gao and Liu, 2014)
                    noise_begin = int(round((window_start - 20) * st[i].meta.sampling_rate) - 1)
                    noise_end = int(round((window_start - 10) * st[i].meta.sampling_rate))
                    signal_begin = int(round((window_start - 8) * st[i].meta.sampling_rate) - 1)
                    signal_end = int(round((window_start + 12) * st[i].meta.sampling_rate))
                    # Calculate vertical component SNR
                    vertical_noise = np.abs(np.mean(st[chz].data[noise_begin:noise_end]))
                    vertical_signal = np.max(np.abs(st[chz].data[signal_begin:signal_end]))
                    vertical_snr = vertical_signal / vertical_noise
                    snr_z[k] = vertical_snr
                    # Calculate radial component SNR
                    radial_noise = np.abs(np.mean(st[i].data[noise_begin:noise_end]))
                    radial_signal = np.max(np.abs(st[i].data[signal_begin:signal_end]))
                    radial_snr = radial_signal / radial_noise
                    snr_r[k] = radial_snr
                    # Calculate receiver function
                    [rf[k], rms[k]] = rfunc.iterdecon(st[i].data, st[chz].data,
                                                      st[i].meta.delta, len(st[i].data),
                                                      tshift, gw, itmax, tol)
                    rad_idx[k] = i
                    k += 1

    # Write receiver functions to obspy stream/trace objects
    rfstream = obspy.Stream(traces=rf[:])
    qc = []
    fit = []
    for i in range(0, len(rfstream)):
        ch1 = rad_idx[i]
        rfstream[i] = obspy.Trace(rf[i])
        rfstream[i].stats = st[ch1].stats
        rfstream[i].stats.sac.b = 0
        # Save vertical SNR to 'USER0' header
        rfstream[i].stats.sac.user0 = snr_z[i]
        # Save radial SNR to 'USER1' header
        rfstream[i].stats.sac.user1 = snr_r[i]
        # Save RMS error to 'USER2' header
        fit.append(i)
        fit[i] = 100*(1 - rms[i][-1])
        rfstream[i].stats.sac.user2 = fit[i]
        # Format filename
        evid = st[ch1].stats.starttime.isoformat().replace('-', '.').replace('T', '.').replace(':', '.').split('.')
        evid.extend([ntwk, stat, loc, 'RF', 'sac'])
        evid = ".".join(evid)
        # Calculate quality ratio and save to 'USER3' header
        qc.append(i)
        qc[i] = rfunc.rf_quality(rf[i], st[ch1].meta.delta, gw, tshift=tshift)
        rfstream[i].stats.sac.user3 = qc[i]
        # Write to SAC files
        rfstream[i].write(rf_dir + evid, format='SAC')

    print('Receiver function computation complete!')
