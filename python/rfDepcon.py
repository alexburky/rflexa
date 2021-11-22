import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
# rfDepcon.py
#
# This function takes a receiver function in the time domain as input and produces a receiver function converted
# to depth as output. It supports a variety of one-dimensional Earth models for depth conversion.
#
# Note: The file 'IASP91.tvel' is required to run this function. It is currently configured to assume that this file is
# located in the current working directory
#
# ----------------------------------------------------------------------------------------------------------------------
# Last updated 11/22/2021 by aburky@princeton.edu
# ----------------------------------------------------------------------------------------------------------------------


def rf_depth_conversion(rf_data, dt, ray_parameter, time_shift, dz, model):

    p = ray_parameter

    if model == 'iasp91':
        iasp91_file = open("IASP91.tvel")
        iasp91 = np.loadtxt(iasp91_file, delimiter=" ", skiprows=2)
        z = iasp91[:, 0][0:121]
        vp = iasp91[:, 1][0:121]
        vs = iasp91[:, 2][0:121]
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
    rf_data = rf_data / np.max(rf_data)

    # Remove time shift
    shidx = int(time_shift/dt)
    rf_data = rf_data[shidx:-1]

    # Depth convert the data!
    ttzm_r = np.round(ttzm / dt)
    depthRF = rf_data[ttzm_r.astype(int)]

    return depthRF
