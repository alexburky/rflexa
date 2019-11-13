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
stat = "ARU"
loc = "00"
gw = 1.0

# sac_dir = "/Users/aburky/PycharmProjects/bermudaRFs/data/rfQuakes/"
sac_dir = "/mnt/usb/aburky/IFILES/STATIONS/" + ntwk + "_" + stat + "/RFQUAKES/"
# rf_dir = "/Users/aburky/PycharmProjects/bermudaRFs/data/rfuncs/"
rf_dir = "/mnt/usb/aburky/IFILES/STATIONS/" + ntwk + "_" + stat + "/RFUNCS/GW" + ''.join(str(gw).split('.')) + "/"
if os.path.exists(rf_dir):
    response = raw_input("Receiver function directory exists. Overwrite? [y/n]:")
    if response == "y":
        shutil.rmtree(rf_dir)
        os.makedirs(rf_dir)
    elif response == "n":
        print('Terminating process...')
        quit()
    else:
        print('Invalid response! Terminating process...')
        quit()
if not os.path.exists(rf_dir):
    os.makedirs(rf_dir)

st = obspy.read(sac_dir + "*.sac")
show_traces = 0

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

show_traces = 0
# Rotate data
ch1 = np.arange(0, len(st), 3)
for i in range(0, len(ch1)):
    ch2 = ch1[i] + 1
    # Check if length of traces is equal!
    if st[ch1[i]].meta.npts != st[ch2].meta.npts:
        npmin = np.min([st[ch1[i]].meta.npts, st[ch2].meta.npts])
        st[ch1[i]].data = st[ch1[i]].data[0:npmin]
        st[ch2].data = st[ch2].data[0:npmin]
    # st[ch1[i]].data, st[ch2].data = su.seisne(st[ch1[i]].data, st[ch2].data, st[ch1[i]].meta.sac.cmpaz)
    st[ch1[i]].data, st[ch2].data = su.seisne(st[ch1[i]].data.astype(float), st[ch2].data.astype(float),
                                              st[ch1[i]].meta.sac.cmpaz)
    # st[ch1[i]].data, st[ch2].data = su.seisrt(st[ch1[i]].data, st[ch2].data, st[ch1[i]].meta.sac.baz)
    st[ch1[i]].data, st[ch2].data = su.seisrt(st[ch1[i]].data.astype(float), st[ch2].data.astype(float),
                                              st[ch1[i]].meta.sac.baz)
    if show_traces == 1:
        plt.ion()
        plt.plot(st[ch1[i]].data, 'k', linewidth=0.25)
        plt.title(st[ch1[i]].meta.starttime)
        plt.show(block=False)
        plt.pause(0.05)
        plt.cla()

show_traces = 0
# Cut and taper data
for i in range(0, len(st)):
    # Get indices of window around P arrival
    bidx = int(round((st[i].meta.sac.t0 - 30) * st[i].meta.sampling_rate) - 1)
    eidx = int(round((st[i].meta.sac.t0 + 90) * st[i].meta.sampling_rate))
    st[i].data = st[i].data[bidx:eidx]
    # Taper data
    window = signal.tukey(eidx - bidx, alpha=0.25)
    st[i].data = st[i].data*window
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
rad_idx = np.arange(0, len(st), 3)
rf = []
rms = []
for i in range(0, len(rad_idx)):
    rf.append(i)
    rms.append(i)
    [rf[i], rms[i]] = rfunc.iterdecon(st[rad_idx[i]].data, st[(rad_idx[i]+2)].data, st[rad_idx[i]].meta.delta,
                                      len(st[rad_idx[i]].data), tshift, gw, itmax, tol)

# Write receiver functions to obspy stream/trace objects
rfstream = obspy.Stream(traces=rf[:])
qc = []
fit = []
for i in range(0, len(rfstream)):
    ch1 = i*3
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
