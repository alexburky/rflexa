import obspy
import seisutils as su
import iterdecon as rfunc
import numpy as np
import matplotlib
from scipy import signal
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

sac_dir = "/Users/aburky/PycharmProjects/bermudaRFs/data/rfQuakes/"
# rf_dir = "/Users/aburky/PycharmProjects/bermudaRFs/data/rfuncs/"
# if os.path.exists(rf_dir):
    # Maybe add an interface/dialogue that checks with user if they would like to overwrite folder?
#    shutil.rmtree(rf_dir)
#    os.makedirs(rf_dir)
# if not os.path.exists(rf_dir):
#    os.makedirs(rf_dir)

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
    st[ch1[i]].data, st[ch2].data = su.seisne(st[ch1[i]].data, st[ch2].data, st[ch1[i]].meta.sac.cmpaz)
    st[ch1[i]].data, st[ch2].data = su.seisrt(st[ch1[i]].data, st[ch2].data, st[ch1[i]].meta.sac.baz)
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
gw = 0.7
itmax = 1000
tol = 0.001
# Deconvolve!
rad_idx = np.arange(0, len(st), 3)
rf = []
rms = []
# rfstream = obspy.Stream()
for i in range(0, len(rad_idx)):
    rf.append(i)
    rms.append(i)
    [rf[i], rms[i]] = rfunc.iterdecon(st[rad_idx[i]].data, st[(rad_idx[i]+2)].data, st[rad_idx[i]].meta.delta,
                                      len(st[rad_idx[i]].data), tshift, gw, itmax, tol)
    # rfstream[i] = obspy.Stream(obspy.Trace(rf[i]))

# Write receiver functions to obspy stream/trace objects
rfstream = obspy.Stream(traces=rf[:])
for i in range(0, len(rfstream)):
    rfstream[i] = obspy.Trace(rf[i])
#    rfstream[i].stats.sac = {}

if show_traces == 1:
    for i in range(0, len(rf)):
        plt.ion()
        plt.plot(rf[i], 'k', linewidth=0.5)
        plt.ylim(-0.5, 1.0)
        plt.show(block=False)
        plt.pause(0.5)
        plt.cla()
