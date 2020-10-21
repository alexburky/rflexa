import numpy as np
import obspy
from obspy.taup import TauPyModel
from scipy import signal
import os
import matplotlib
from matplotlib import rc
# matplotlib.use('Qt4Agg')
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# This script makes a contour plot of receiver function data in the time domain as a function of epicentral
# distance
# ------------------------------------------------------------------------------------------
# Last updated 11/12/2019 by aburky@princeton.edu
# ------------------------------------------------------------------------------------------

# Read in receiver function data
# data_directory = "/mnt/usb/aburky/IFILES/NETWORKS/"
data_directory = "/Users/aburky/IFILES/NETWORKS/"
ntwk = "IU"
stat = "BBSR"
loc = "00"
gw = 1.0
qc = 0.05
fit = 90
# Undo time-shift if it existed
tshift = 10

# Construct path to receiver funcitons
rf_dir = data_directory + ntwk + "/" + stat + "/" + loc + "/RFUNCS/FILTERED_0.04_5/GW" + ''.join(str(gw).split('.')) + "/"
rfs = obspy.read(rf_dir + '*.sac')
# Construct path to figure
fig_dir = data_directory + ntwk + "/" + stat + "/" + loc + "/GRAPHICS/"
fig_name = ntwk + "." + stat + ".GW" + ''.join(str(gw).split('.')) + ".QC." + str(qc) + ".FIT." + str(fit) + ".eps"

# Make sure figure directory exists
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

# Get some necessary header information
npts = []
delta = []
bad = []
k = 0
for i in range(0, len(rfs)):
    npts.append(i)
    npts[i] = rfs[i].stats.sac.npts
    delta.append(i)
    delta[i] = rfs[i].stats.sac.delta
    # Check if receiver function is anomalously long (make this automatic!)
    if rfs[i].stats.sac.e >= 150:
        bad.append(k)
        bad[k] = i
        k += 1

# Delete bad receiver functions
bad.sort(reverse=True)
for i in range(0, len(bad)):
    del npts[bad[i]], delta[bad[i]], rfs[bad[i]]

# Resample data for consistent plotting
npmax = np.max(npts)
dt = np.min(delta)
norm_rfs = []
for i in range(0, len(rfs)):
    if rfs[i].stats.sac.npts != npmax:
        rfs[i].data = signal.resample(rfs[i].data, npmax)
    # Normalize receiver functions (maximum amplitude of 1)
    norm_rfs.append(i)
    norm_rfs[i] = np.array(rfs[i].data/np.max(rfs[i].data))
norm_rfs = np.squeeze(norm_rfs)

# Read in Gao travel time table
# pdstime = np.loadtxt(r"/mnt/usb/aburky/pds.ttime")
pdstime = np.loadtxt(r"/Users/aburky/pds.ttime")
# Construct list of time to depth conversion values
ttzm = []
ags = np.zeros((1, 801))
j = 0
for i in range(0, len(rfs)):
    if rfs[i].stats.sac.user0 > fit and rfs[i].stats.sac.user1 > qc:
        ttzm.append(j)
        idx = np.where(np.around(pdstime[:, 0], 2) == round(rfs[i].stats.sac.user9, 2))
        ttzm[j] = pdstime[idx, 2]
        # ai = rfs[i].data[np.round(ttzm[j]/rfs[i].stats.sac.delta).astype(int)]
        ai = norm_rfs[i][(np.round(ttzm[j]/dt).astype(int)+np.round(10/dt).astype(int))]
        ags += ai
        j += 1

# Set matplotlib to use LaTeX fonts
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
# Plot the result
x = np.arange(0, 801, 1)
plt.plot(x, ags[0]/len(rfs), 'k', linewidth=0.25)
plt.fill_between(x, 0, ags[0]/len(rfs), where=ags[0]/len(rfs) > 0, facecolor='red')
plt.xlim(200, 800)
plt.ylim(-0.01, 0.01)
plt.xlabel('Depth (km)')
plt.ylabel('Amplitude (relative to P)')
plt.title(str(j))
plt.show()
