from rfDepcon import rf_depth_conversion
import obspy
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
# timeToDepth.py
#
# This script gives an example of how to convert a single time-domain receiever function into a depth-domain receiver
# function using the function found in 'rfDepcon.py'
#
# ----------------------------------------------------------------------------------------------------------------------
# Last updated 11/22/2021 by aburky@princeton.edu
# ----------------------------------------------------------------------------------------------------------------------

# Path to the time domain receiver functions
data_dir = '/Users/aburky/IFILES/NETWORKS/TA/F36A/NULL/RFUNCS_VEL/FILTERED_0.02_0.2/GW10/'

# Read in all of the seismic data in the directory
st = obspy.read(data_dir + "*.SAC")

# Necessary metadata for depth conversion (sample rate, ray parameter)
dt = st[0].stats.delta
deg2km = 111.1949
p = st[0].stats.sac.user9/deg2km

# Depth resolution (km)
dz = 0.1

# Remove 10 second time shift (default option for computeRFs.py)
time_shift = 10

# Depth convert the first receiver function
depth_rf = rf_depth_conversion(st[0].data, dt, p, time_shift, dz, model='iasp91')
z = np.arange(0, len(depth_rf)/10, dz)

# Plot the original and depth converted receiver functions
plt.subplot(211)
plt.plot(st[0].times(), st[0].data, 'k')
plt.xlim(0, 120)
plt.xlabel('Time (s)')

plt.subplot(212)
plt.plot(z, depth_rf, 'r')
plt.xlim(0, 800)
plt.xlabel('Depth (km)')
