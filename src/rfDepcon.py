import obspy
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------------------------------------------------
# rfDepcon.py
#
# This function takes a receiver function in the time domain as input and produces a receiver function converted
# to depth as output. It supports a variety of one-dimensional Earth models for depth conversion.
# ---------------------------------------------------------------------------------------------------------------------
# Last updated 8/18/2020 by aburky@princeton.edu
# ---------------------------------------------------------------------------------------------------------------------

# Directory containing data
sac_dir = '/Users/aburky/IFILES/NETWORKS/IU/BBSR/00/RFUNCS/UNFILTERED/GW10/'

# File to read
event = '2020.08.16.20.05*.sac'
st = obspy.read(sac_dir + event)

deg2km = 111.1949

dt = st[0].meta.delta
p = st[0].meta.sac.user9/deg2km
dz = 0.1
model = 'iasp91'

if model == 'iasp91':
    iasp91_file = open("/Users/aburky/IFILES/MODELS/IASP91/IASP91.tvel")
    iasp91 = np.loadtxt(iasp91_file, delimiter=" ", skiprows=2)
    z = iasp91[:, 0][0:124]
    vp = iasp91[:, 1][0:124]
    vs = iasp91[:, 2][0:124]
elif model == 'prem':
    print('Model PREM is currently unsupported... :(')
else:
    print('ERROR: Not a valid velocity model. Acceptable values are \'iasp91\' or \'prem\'')

# Verify that the user input a valid ray parameter
if p > 1:
    print('ERROR: Invalid ray parameter, units should be (s/km)')
    quit()

# Find discontinuities in the Earth model
idx = np.where(np.diff(z) == 0)[0] + 1

zp = []
vpp = []
vsp = []
# Handle discontinuities in the Earth model
for i in range(0, len(idx) + 1):
    zp.append(i)
    vpp.append(i)
    vsp.append(i)
    if i == 0:
        zp[i] = z[i:idx[i]]
        vpp[i] = vp[i:idx[i]]
        vsp[i] = vs[i:idx[i]]
    elif 0 < i < len(idx):
        zp[i] = z[idx[i - 1]:idx[i]]
        vpp[i] = vp[idx[i - 1]:idx[i]]
        vsp[i] = vs[idx[i - 1]:idx[i]]
    elif i == len(idx):
        zp[i] = z[idx[i - 1]:]
        vpp[i] = vp[idx[i - 1]:]
        vsp[i] = vs[idx[i - 1]:]

vp_int = []
vs_int = []
zzp = []
# Construct interpolants
for i in range(0, len(idx)):
    vp_int.append(i)
    vs_int.append(i)
    zzp.append(i)
    vp_int[i] = np.interp(np.arange(zp[i][0], zp[i][-1] + dz, dz), zp[i], vpp[i])
    vs_int[i] = np.interp(np.arange(zp[i][0], zp[i][-1] + dz, dz), zp[i], vsp[i])
    zzp[i] = np.arange(zp[i][0], zp[i][-1] + dz, dz)

# Construct depth conversion integrand terms (from Chevrot et al.)
a = []
b = []
v = []
for i in range(0, len(idx)):
    a.append(i)
    b.append(i)
    v.append(i)
    a[i] = np.sqrt(pow(vs_int[i], -2) - pow(p, 2) * pow(((6371 - zzp[i]) / 6371), -2))
    b[i] = np.sqrt(pow(vp_int[i], -2) - pow(p, 2) * pow(((6371 - zzp[i]) / 6371), -2))
    v[i] = a[i] - b[i]

# Solve the integral piecewise
tz = []
for i in range(0, len(idx)):
    tz.append(i)
    tz[i] = []
    for j in range(0, len(zzp[i])):
        tz[i].append(j)
        tz[i][j] = np.trapz(v[i][0:j], zzp[i][0:j])

# Concatenate into a time to depth conversion vector
ttd = tz
for i in range(0, len(idx) - 1):
    ttd[i + 1] = ttd[i + 1] + ttd[i][-1]

ttzm = np.concatenate(ttd)

# Delete duplicate entries in ttzm (corresponding to discontinuities)
dup_idx = np.where(np.diff(ttzm) == 0)
ttzm = np.delete(ttzm, dup_idx)

# Keep real components of ttzm
if not np.isrealobj(ttzm):
    ttzm = np.real(ttzm)
    print('Warning: Depth conversion matrix has imaginary values')

# Normalize receiver function amplitude
st[0].data = st[0].data / np.max(st[0].data)

# Align maximum to zero time
shidx = np.argmax(st[0].data)
st[0].data = st[0].data[shidx:-1]

# Depth convert the data!
ttzm_r = np.round(ttzm / dt)

depthRF = st[0].data[ttzm_r.astype(int)]

# Plot and verify?
depth = np.arange(0, 799.6, dz)

plt.plot(depth, depthRF)
plt.show()
