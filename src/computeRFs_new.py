import obspy
import seisutils as su
import iterdecon as rfunc
import numpy as np
import matplotlib
import os
import shutil
from scipy import signal
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

# Define network, station and location
ntwk = "II"
stat = "OBN"
loc = "00"
gw = 1.0

show_traces = 0

# sac_dir = "/Users/aburky/PycharmProjects/bermudaRFs/data/rfQuakes/"
sac_dir = "/mnt/usb/aburky/IFILES/STATIONS/" + ntwk + "_" + stat + "/RFQUAKES/"
# rf_dir = "/Users/aburky/PycharmProjects/bermudaRFs/data/rfuncs/"
rf_dir = "/mnt/usb/aburky/IFILES/STATIONS/" + ntwk + "_" + stat + "/RFUNCS/GW" + ''.join(str(gw).split('.')) + "/"
if os.path.exists(rf_dir):
    # response = raw_input("Receiver function directory exists. Overwrite? [y/n]:")
    # if response == "y":
        shutil.rmtree(rf_dir)
        os.makedirs(rf_dir)
    # elif response == "n":
    #     print('Terminating process...')
    #     quit()
    # else:
    #     print('Invalid response! Terminating process...')
    #     quit()
if not os.path.exists(rf_dir):
    os.makedirs(rf_dir)

st = obspy.read(sac_dir + "*.sac")

# Filter data
for i in range(0, len(st)):
    # Test with rounded vs. un-rounded data
    fs = 1/st[i].meta.sac.delta
    lo = 0.034
    hi = 0.1
    st[i].data = su.bp_butter(st[i].data, lo, hi, fs, 3)
    if show_traces == 1:
        plt.ion()
        plt.plot(st[i].data, 'k', linewidth=0.25)
        plt.title(st[i].meta.starttime)
        plt.show(block=False)
        plt.pause(0.05)
        plt.cla()

# Rotate data
for i in range(0, len(st)):
    if (st[i].stats.sac.kcmpnm == u'BH1') or (st[i].stats.sac.kcmpnm == u'BHN'):
        tb = st[i].meta.starttime
        te = st[i].meta.endtime
        for j in range(-2, 3):
            if (st[(i+j)].stats.sac.kcmpnm == u'BH2') or (st[(i+j)].stats.sac.kcmpnm == u'BHE') and \
               (tb - 0.5 <= st[(i+j)].meta.starttime <= tb + 0.5 or te - 0.5 <= st[(i+j)].meta.endtime <= te +0.5):
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
    # Check to make sure data exists
    #if st[i].data.size != st[i].meta.sac.npts:
    #    print('Array size is problematic, moving on to the next trace...')
    #else:
        # Get indices of window around P arrival
    bidx = int(round((st[i].meta.sac.t0 - 30) * st[i].meta.sampling_rate) - 1)
    eidx = int(round((st[i].meta.sac.t0 + 90) * st[i].meta.sampling_rate))
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
k = 0
for i in range(0, len(st)):
    if (st[i].stats.sac.kcmpnm == u'BH1') or (st[i].stats.sac.kcmpnm == u'BHN'):
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
    # Save RMS error to 'USER0' header
    fit.append(i)
    fit[i] = 100*(1 - rms[i][-1])
    rfstream[i].stats.sac.user0 = fit[i]
    # Format filename
    evid = st[ch1].stats.starttime.isoformat().replace('-', '.').replace('T', '.').replace(':', '.').split('.')
    evid.extend([ntwk, stat, loc, 'RF', 'sac'])
    evid = ".".join(evid)
    # Calculate quality ratio and save to 'USER1' header
    qc.append(i)
    qc[i] = rfunc.rf_quality(rf[i], st[ch1].meta.delta, gw, tshift=tshift)
    rfstream[i].stats.sac.user1 = qc[i]
    # Write to SAC files
    rfstream[i].write(rf_dir + evid, format='SAC')

print('Receiver function computation complete!')
# Make scatter plot showing statistics
# plt.scatter(fit, qc)
# plt.show()

# for i in range(0, len(rf)):
#     if rfstream[i].stats.sac.user1 > 0.35:
#         plt.ion()
#         plt.plot(rf[i], 'k', linewidth=0.5)
#         plt.title(evid + ' Quality: ' + str(rfstream[i].stats.sac.user0))
#         plt.ylim(-0.5, 1.0)
#         plt.show(block=False)
#         plt.pause(0.5)
#         plt.cla()
