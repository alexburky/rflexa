import obspy
import numpy as np
from rfDepcon import rf_depth_conversion

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------------------------------------------------
# stackRFs.py
#
# This script depth converts and stacks receiver functions, and saves the result to a .mat file for subsequent
# plotting with a MATLAB script.
#
# ---------------------------------------------------------------------------------------------------------------------
# Last updated 8/18/2020 by aburky@princeton.edu
# ---------------------------------------------------------------------------------------------------------------------

# Specify the location of the receiver function data
rf_dir = '/Users/aburky/IFILES/NETWORKS/IU/BBSR/00/RFUNCS/UNFILTERED/GW10/*.sac'
st = obspy.read(rf_dir)

# Specify QC constraints
snr_z = 4
snr_r = 4
nu = 0.2
fit = 80
gcarc = 29

j = 0
rf_depth = []
stack = []
for i in range(0, len(st)):
    if st[i].meta.sac.user0 > snr_z and st[i].meta.sac.user1 > snr_r and st[i].meta.sac.user2 > fit  \
     and st[i].meta.sac.user3 > nu and st[i].meta.sac.gcarc >= gcarc:
        rf_depth.append(j)
        stack.append(0)
        P = st[i].meta.sac.user9 / 111.1949
        rf_depth[j] = rf_depth_conversion(st[i].data, st[i].meta.delta, P, 0.1, 'iasp91')
        stack[0] = stack[0] + rf_depth[j]
        j = j + 1

stack[0] = stack[0] / j

depth = np.arange(0, 759.6, 0.1)
plt.plot(depth, stack[0])
plt.show()
