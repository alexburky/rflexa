import obspy
import numpy as np
from rfDepcon import rf_depth_conversion
import scipy.io as scio

# ---------------------------------------------------------------------------------------------------------------------
# stackRFs.py
#
# This script depth converts and stacks receiver functions, and saves the result to a .mat file for subsequent
# plotting with a MATLAB script.
#
# ---------------------------------------------------------------------------------------------------------------------
# Last updated 9/8/2020 by aburky@princeton.edu
# ---------------------------------------------------------------------------------------------------------------------

# Specify QC constraints
snr_z = 4
snr_r = 4
nu = 0.2
fit = 80
gcarc = 29

# Specify range of gaussian widths for stacking
# gaussian_width = np.arange(0.5, 1.55, 0.05)

# Alternatively, specify range of filter corners for stacking
fc_high = np.unique(np.around(1.0/np.arange(1, 50, 1), 3))

# Specify the directory containing the receiver function data
# data_directory = '/Users/aburky/IFILES/NETWORKS/IU/BBSR/00/RFUNCS/UNFILTERED/'

# Loop over directories and generate depth converted stacks
# for gw in gaussian_width:
for fc in fc_high:
    data_directory = '/home/aburky/IFILES/NETWORKS/IU/BBSR/00/RFUNCS/FILTERED_0.02_' + str(fc) + '/'
    print(data_directory)
    rf_dir = data_directory + '/*.sac'
    st = obspy.read(rf_dir)
    j = 0
    rf_depth = []
    stack = []
    gcarcs = []
    for i in range(0, len(st)):
        if st[i].meta.sac.user0 > snr_z and st[i].meta.sac.user1 > snr_r and st[i].meta.sac.user2 > fit \
         and st[i].meta.sac.user3 > nu and st[i].meta.sac.gcarc >= gcarc:
            rf_depth.append(j)
            stack.append(0)
            P = st[i].meta.sac.user9 / 111.1949
            rf_depth[j] = rf_depth_conversion(st[i].data, st[i].meta.delta, P, 0.1, 'iasp91')
            stack[0] = stack[0] + rf_depth[j]
            gcarcs.append(j)
            gcarcs[j] = st[i].meta.sac.gcarc
            j = j + 1

    # Normalize amplitude of the stack
    stack[0] = stack[0] / j

    # Save the resulting stack to a dictionary
    stack_dict = {'stackLF.0.02.HF.' + str(fc): stack[0],
                  'nEventsLF.0.02.HF.' + str(fc): j,
                  'gcarcsLF.0.02.HF.' + str(fc): gcarcs}

    # Save the dictionary to a .mat file
    scio.savemat('stackLF.0.02.HF.' + str(fc) + '.mat', stack_dict)

    # gw_string = ''.join(str(gw).split('.'))
    # print('Stacking data for Receiver Functions with Gaussian Width: ' + gw_string)
    # rf_dir = data_directory + 'GW' + gw_string + '/*.sac'
    # st = obspy.read(rf_dir)
    # j = 0
    # rf_depth = []
    # stack = []
    # gcarcs = []
    # for i in range (0, len(st)):
    #     if st[i].meta.sac.user0 > snr_z and st[i].meta.sac.user1 > snr_r and st[i].meta.sac.user2 > fit \
    #      and st[i].meta.sac.user3 > nu and st[i].meta.sac.gcarc >= gcarc:
    #         rf_depth.append(j)
    #         stack.append(0)
    #         P = st[i].meta.sac.user9 / 111.1949
    #         rf_depth[j] = rf_depth_conversion(st[i].data, st[i].meta.delta, P, 0.1, 'iasp91')
    #         stack[0] = stack[0] + rf_depth[j]
    #         gcarcs.append(j)
    #         gcarcs[j] = st[i].meta.sac.gcarc
    #         j = j + 1
    # # Normalize amplitude of stack
    # stack[0] = stack[0] / j
    #
    # # Save the resulting stack to a dictionary
    # stack_dict = {'stack_GW' + gw_string: stack[0],
    #               'nEvents_GW' + gw_string: j,
    #               'gcarcs_GW' + gw_string: gcarcs}
    #
    # # Save the dictionary to a .mat file
    # scio.savemat('stack_GW' + gw_string + '.mat', stack_dict)

