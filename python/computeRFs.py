import obspy
import seisutils as su
import iterdecon as rfunc
import numpy as np
import os
import shutil
import glob
from scipy import signal
import time
import fnmatch

# ------------------------------------------------------------------------------------------
# computeRFs.py
#
# This function makes a receiver function and saves it to a directory upon completion. The
# user can specify a choice of filter and Gaussian width factors for use in the data
# processing. The resulting SAC files are saved with additional metadata in the header.
#
# USER0: Vertical component SNR
# USER1: Radial component SNR
# USER2: RMS error calculated during iterative time domain deconvolution
# USER3: Quality metric, nu
#
# To do: Add more options for instrument response removal, this is currently hardcoded
#
# ------------------------------------------------------------------------------------------
# Last updated 11/12/2021 by aburky@princeton.edu
# ------------------------------------------------------------------------------------------


def computeRFs(network, location, data_directory, input_units, gaussian_width, high_cut, low_cut, station):

    # Log the execution time for statistical purposes
    t1 = time.time()

    # Define network, station, and location
    ntwk = network
    stat = station
    loc = location
    gw = gaussian_width

    print('Calculating receiver functions for station: ' + stat)

    # If location is "ALL" - need to loop over all location folders
    dirs = [location]
    if loc == "ALL":
        root, dirs, files = os.walk(data_directory + ntwk + "/" + stat + "/").__next__()

    # Loop over all possible locations
    for locDir in np.arange(0, len(dirs)):

        loc = dirs[locDir]

        # Determine directory containing the earthquake data
        sac_dir = data_directory + ntwk + "/" + stat + "/" + loc + "/RFQUAKES_COUNTS/"

        # Construct path to output directory (this needs to be dependent on desired output units!)
        rf_dir = data_directory + ntwk + "/" + stat + "/" + loc + "/"
        if input_units == 'displacement':
            rf_dir = rf_dir + "RFUNCS_DISP/"
        elif input_units == 'velocity':
            rf_dir = rf_dir + "RFUNCS_VEL/"
        elif input_units == 'acceleration':
            rf_dir = rf_dir + "RFUNCS_ACC/"
        else:
            print('Error: Invalid output units selected, acceptable options are displacement, velocity, or acceleration')
            return

        if low_cut != 0 and high_cut != 0:
            rf_dir = rf_dir + "FILTERED_" + str(low_cut) + "_" + str(high_cut) + "/GW" + ''.join(str(gw).split('.')) + "/"
        else:
            rf_dir = rf_dir + "UNFILTERED/GW" + ''.join(str(gw).split('.')) + "/"

        # Check if the directory already exists
        if os.path.exists(rf_dir):
            # if glob.glob(rf_dir + '*.log'):
            #     print('Directory with receiver functions already exists, terminating process...')
            #     return
            # else:
            shutil.rmtree(rf_dir)
        if not os.path.exists(rf_dir):
            os.makedirs(rf_dir)

        # Log potential errors to a .log file
        logFile = rf_dir + ntwk + '.' + stat + '.RF.log'

        # Read in all of the seismic data in the directory
        st = obspy.read(sac_dir + "*.SAC")

        # Loop over data and perform processing steps
        for i in range(0, len(st)):
            # First, remove the instrument response
            pz_id = str(st[i].meta.sac.user0).split('.')[0]
            pz_file = glob.glob(sac_dir + 'SAC_PZs*' + pz_id)[0]

            # Remove the mean, trend, and taper!
            st[i].data = st[i].data - np.mean(st[i].data)
            st[i].data = signal.detrend(st[i].data, type='linear')
            taper = signal.tukey(len(st[i].data), alpha=0.1)
            st[i].data = st[i].data*taper

            # Set filter parameters
            f = [0.002, 0.004, 5.0, 10.0]

            st[i].data = su.transfer(st[i].data, st[i].meta.sac.delta, f, input_units, pz_file, file_type='sacpz')

        # Filter data
        if low_cut != 0 and high_cut != 0:
            for i in range(0, len(st)):
                fs = 1/st[i].meta.sac.delta
                st[i].data = su.bp_butter(st[i].data, low_cut, high_cut, fs, 3)

        # Rotate data
        for i in range(0, len(st)):
            if (fnmatch.fnmatch(st[i].stats.sac.kcmpnm, '?H1')) or (fnmatch.fnmatch(st[i].stats.sac.kcmpnm, '?HN')):
                tb = st[i].meta.starttime
                te = st[i].meta.endtime
                for j in range(-2, 3):
                    if (fnmatch.fnmatch(st[(i+j)].stats.sac.kcmpnm, '?H2')) or (fnmatch.fnmatch(st[(i+j)].stats.sac.kcmpnm, '?HE')) and \
                            (tb - 0.5 <= st[(i+j)].meta.starttime <= tb + 0.5 or te - 0.5 <= st[(i+j)].meta.endtime <= te + 0.5):
                        ch1 = i
                        ch2 = int(i + j)
                        # Check to make sure length of data vectors is equal
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
                        with open(logFile, "a") as log:
                            log.write('Couldn''t find matching horizontal data for event: ' + str(tb) + '\n')

        # Cut and taper data
        for i in range(0, len(st)):
            window_start = 30
            window_end = 90
            bidx = int(round((st[i].meta.sac.t0 - window_start) * st[i].meta.sampling_rate) - 1)
            eidx = int(round((st[i].meta.sac.t0 + window_end) * st[i].meta.sampling_rate))
            # print('Data Size: %i' % st[i].data.size)
            # print('bidx: %i' % bidx)
            # print('eidx: %i' % eidx)
            if eidx > st[i].data.size:
                # print('Array size is problematic, moving on to the next trace...')
                with open(logFile, "a") as log:
                    log.write('Array size is problematic for event: ' + str(st[i].meta.starttime) + '\n')
                continue
            else:
                st[i].data = st[i].data[bidx:eidx]
                # Taper data
                window = signal.tukey(eidx - bidx, alpha=0.25)
                st[i].data = st[i].data * window

        # Set deconvolution parameters
        tshift = 10
        itmax = 1000
        tol = 0.001
        rf = []
        rms = []
        rad_idx = []
        snr_z = []
        snr_r = []
        k = 0

        for i in range(0, len(st)):
            if (fnmatch.fnmatch(st[i].stats.sac.kcmpnm, '?H1')) or (fnmatch.fnmatch(st[i].stats.sac.kcmpnm, '?HN')):
                if (fnmatch.fnmatch(st[i].stats.sac.kcmpnm, '?H1')):
                    bidx = -2
                    eidx = 3
                else:
                    bidx = -3
                    eidx = 2
                tb = st[i].meta.starttime
                te = st[i].meta.endtime
                for j in range(bidx, eidx):
                    if (fnmatch.fnmatch(st[(i+j)].stats.sac.kcmpnm, '?HZ')) and \
                            (tb - 0.5 <= st[(i+j)].meta.starttime <= tb + 0.5 or
                             te - 0.5 <= st[(i+j)].meta.endtime <= te + 0.5) and (st[i].data.size == st[(i+j)].data.size):

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
                        [rf[k], rms[k]] = rfunc.iterdecon(st[i].data, st[chz].data, st[i].meta.delta, len(st[i].data),
                                                          tshift, gw, itmax, tol)
                        rad_idx[k] = i
                        k += 1

        # Write receiver functions to obspy stream/trace objects
        rfstream = obspy.Stream(traces=rf[:])
        for i in range(0, len(rfstream)):
            ch1 = rad_idx[i]
            rfstream[i] = obspy.Trace(rf[i])
            rfstream[i].stats = st[ch1].stats
            rfstream[i].stats.sac.b = 0
            # Save vertical SNR to 'USER0' Header
            rfstream[i].stats.sac.user0 = snr_z[i]
            # Save radial SNR to 'USER1' Header
            rfstream[i].stats.sac.user1 = snr_r[i]
            # Save RMS Error to 'USER2' Header
            fit = 100*(1 - rms[i][-1])
            rfstream[i].stats.sac.user2 = fit
            # Format filename
            evid = st[ch1].stats.starttime.isoformat().replace('-', '.').replace('T', '.').replace(':', '.').split('.')[:-1]
            # loc = st[ch1].stats.sac.khole
            evid.extend([ntwk, stat, loc, 'RF', 'SAC'])
            evid = ".".join(evid)
            # Calculate quality control ratio and save to 'USER3' Header
            qc = rfunc.rf_quality(rf[i], st[ch1].meta.delta, gw, tshift=tshift)
            rfstream[i].stats.sac.user3 = qc
            # Write to SAC files
            rfstream[i].write(rf_dir + evid, format='SAC')

        elapsed = time.time() - t1
        with open(logFile, "a") as log:
            log.write('Time required to generate receiver functions: ' + str(elapsed) + 's')

        # print('Receiver function computation complete!')

