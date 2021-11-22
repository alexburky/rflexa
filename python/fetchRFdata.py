from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
import numpy as np
import seisutils as su
import os
import time

# --------------------------------------------------------------------------------------------------------------
# newFetch.py
#
# This is an updated version of fetchRFdata.py which allows the user to download the data without performing
# any pre-processing. The instrument response information is saved in both 'RESP' and 'SACPZ' output formats.
# Since downloading data is often slower than processing it, all subsequent processing is done with subsequent
# scripts.
#
# TO DO: Revise this so that it can get all locations - save them to specific directories
#
# --------------------------------------------------------------------------------------------------------------
# Last updated 11/22/2021 by aburky@princeton.edu
# --------------------------------------------------------------------------------------------------------------


def fetch_rf_data(network, location, channel, data_directory, output_units, minimum_magnitude,
                  maximum_magnitude, station):

    # Track execution time for logging purposes
    t1 = time.time()

    ntwk = network
    stat = station
    loc = location
    chan = channel

    # Define the client that hosts the desired data
    client = Client("IRIS")

    # # Define directory where seismic data will be saved as SAC files
    # if output_units == 'counts':
    #     # sac_dir = data_directory + ntwk + '/' + stat + '/' + loc + '/RFQUAKES_COUNTS/'
    # elif output_units == 'displacement':
    #     # sac_dir = data_directory + ntwk + '/' + stat + '/' + loc + '/RFQUAKES_DISP/'
    # elif output_units == 'velocity':
    #     # sac_dir = data_directory + ntwk + '/' + stat + '/' + loc + '/RFQUAKES_VEL/'
    # elif output_units == 'acceleration':
    #     # sac_dir = data_directory + ntwk + '/' + stat + '/' + loc + '/RFQUAKES_ACC/'
    # else:
    #     print('ERROR: Invalid output units. Acceptable options are \'displacement,\' \'velocity,\' or \'counts\'')
    #     quit()

    # # For now: delete the directory if it exists...
    # if os.path.exists(sac_dir):
    #     print('Directory exists. Terminiating process...')
    #     return
    #     # shutil.rmtree(sac_dir)
    #
    # if not os.path.exists(sac_dir):
    #     os.makedirs(sac_dir)

    # Define amount of data desired (minutes)
    duration = 60

    # Log potential errors to a .log file
    # logFileName = sac_dir + ntwk + '.' + stat + '.log'

    # Fetch station information for data retrieval
    if loc == "ALL":
        try:
            inv = client.get_stations(network=ntwk, station=stat, channel=chan, level="response")
        except Exception as error:
            print('Error fetching station information with the IRIS client...')
            print(str(error))
            return
            # with open(logFileName, "a") as log:
            #     log.write(str(error))
            #     log.write('Error fetching station information with the IRIS client...')
            # return
    else:
        try:
            inv = client.get_stations(network=ntwk, station=stat, loc=loc, channel=chan, level="response")
        except Exception as error:
            print('Error fetching station information with the IRIS client...')
            print(str(error))
            return
            # with open(logFileName, "a") as log:
            #     log.write(str(error))
            #     log.write('Error fetching station information with the IRIS client...')
            # return

    # Save the pole zero files
    nstats = len(inv.networks[0])
    # resp_t0 = []
    # resp_tf = []
    # pre_filt = []
    for i in range(0, nstats):
        nresp = len(inv.networks[0].stations[i].channels)
        # Tag the PZ files and SAC files with a number indicating the period of operation
        for j in range(0, nresp):
            loc = inv.networks[0].stations[i].channels[j].location_code
            if loc == '':
                sac_dir = data_directory + ntwk + '/' + stat + '/NULL/RFQUAKES_COUNTS/'
            else:
                sac_dir = data_directory + ntwk + '/' + stat + '/' + loc + '/RFQUAKES_COUNTS/'

            # Make the directory if it doesn't exist
            if not os.path.exists(sac_dir):
                os.makedirs(sac_dir)

            fileName = sac_dir + "SAC_PZs_" + ntwk + '_' + stat + '_' + inv.networks[0].stations[i].channels[j].code + \
                       '.' + str(j)
            with open(fileName, "a") as pzFile:
                pzFile.write('* **********************************\n')
                pzFile.write('* NETWORK   (KNETWK): ' + inv.networks[0].code + '\n')
                pzFile.write('* STATION    (KSTNM): ' + inv.networks[0].stations[i].code + '\n')
                pzFile.write('* LOCATION   (KHOLE): ' + inv.networks[0].stations[i].channels[j].location_code + '\n')
                pzFile.write('* CHANNEL   (KCMPNM): ' + inv.networks[0].stations[i].channels[j].code + '\n')
                pzFile.write('* CREATED           : ' + str(UTCDateTime.now()).split('.')[0] + '\n')
                pzFile.write('* START             : ' +
                             str(inv.networks[0].stations[i].channels[j].start_date).split('.')[0] + '\n')
                pzFile.write('* END               : ' +
                             str(inv.networks[0].stations[i].channels[j].end_date).split('.')[0] + '\n')
                pzFile.write('* DESCRIPTION       : ' + inv.networks[0].stations[i].site.name + '\n')
                pzFile.write('* LATITUDE          : %0.6f\n' % inv.networks[0].stations[i].latitude)
                pzFile.write('* LONGITUDE         : %0.6f\n' % inv.networks[0].stations[i].longitude)
                pzFile.write('* ELEVATION         : %0.1f\n' % inv.networks[0].stations[i].channels[j].elevation)
                pzFile.write('* DEPTH             : %0.1f\n' % inv.networks[0].stations[i].channels[j].depth)
                pzFile.write('* DIP               : %0.1f\n' %
                             (90.0 - np.abs(inv.networks[0].stations[i].channels[j].dip)))
                pzFile.write('* AZIMUTH           : %0.1f\n' % inv.networks[0].stations[i].channels[j].azimuth)
                pzFile.write('* SAMPLE RATE       : %0.1f\n' % inv.networks[0].stations[i].channels[j].sample_rate)
                pzFile.write('* INPUT UNIT        : M\n')
                pzFile.write('* OUTPUT UNIT       : COUNTS\n')
                pzFile.write('* INSTTYPE          : ' +
                             inv.networks[0].stations[i].channels[j].sensor.description +'\n')
                pzFile.write('* INSTGAIN          : %e (M/S)\n' %
                             inv.networks[0].stations[i].channels[j].response.get_paz().stage_gain)
                pzFile.write('* COMMENT           : \n')
                pzFile.write('* SENSITIVITY       : %e (M/S)\n' %
                             inv.networks[0].stations[i].channels[j].response.instrument_sensitivity.value)
                pzFile.write('* A0                : %e\n' %
                             inv.networks[0].stations[i].channels[j].response.get_paz().normalization_factor)
                pzFile.write('* **********************************\n')

                # Save the poles, zeros, and constant
                nzeros = 3
                zeros = inv.networks[0].stations[i].channels[j].response.get_paz().zeros
                nz = np.nonzero(zeros)
                pzFile.write('ZEROS   ' + str(len(nz[0]) + nzeros) + '\n')
                pzFile.write("        %+e   %+e\n" % (0, 0))
                pzFile.write("        %+e   %+e\n" % (0, 0))
                pzFile.write("        %+e   %+e\n" % (0, 0))
                if len(nz[0]) != 0:
                    for k in range(0, len(nz[0])):
                        pzFile.write("        %+e   %+e\n" % (np.real(zeros[nz[0][k]]), np.imag(zeros[nz[0][k]])))

                poles = inv.networks[0].stations[i].channels[j].response.get_paz().poles
                pzFile.write('POLES   ' + str(len(poles)) + '\n')
                for k in range(0, len(poles)):
                    pzFile.write("        %+e   %+e\n" %
                                 (np.real(inv.networks[0].stations[i].channels[j].response.get_paz().poles[k]),
                                  np.imag(inv.networks[0].stations[i].channels[j].response.get_paz().poles[k])))

                pzFile.write('CONSTANT        %e' %
                             (inv.networks[0].stations[i].channels[j].response.get_paz().normalization_factor *
                              inv.networks[0].stations[i].channels[j].response.instrument_sensitivity.value))
                # pzFile.write(inv.networks[0].stations[i].channels[j].response.get_sacpz())

    # Loop over time-periods during which the station was operational and fetch data
    for i in range(0, nstats):
        for j in range(0, nresp):
            # Location code
            loc = inv.networks[0].stations[i].channels[j].location_code
            if loc == '':
                sac_dir = data_directory + ntwk + '/' + stat + '/NULL/RFQUAKES_COUNTS/'
            else:
                sac_dir = data_directory + ntwk + '/' + stat + '/' + loc + '/RFQUAKES_COUNTS/'
            logFileName = sac_dir + ntwk + '.' + stat + '.log'

            if inv.networks[0].stations[i].channels[j].end_date > UTCDateTime.now():
                t0 = inv.networks[0].stations[i].channels[j].start_date
                tf = UTCDateTime.now()
            else:
                t0 = inv.networks[0].stations[i].channels[j].start_date
                tf = inv.networks[0].stations[i].channels[j].end_date
            # Get station coordinates for event selection
            stla = inv.networks[0].stations[i].latitude
            stlo = inv.networks[0].stations[i].longitude
            # Fetch relevant events in time-window during which station was operational
            try:
                catalog = client.get_events(starttime=t0, endtime=tf, minmagnitude=minimum_magnitude,
                                            maxmagnitude=maximum_magnitude, latitude=stla, longitude=stlo, minradius=30,
                                            maxradius=90)
            except Exception as error:
                with open(logFileName, "a") as log:
                    log.write(str(error))
                    log.write('Error fetching event catalog...')
                continue

            nEvents = len(catalog.events)
            # Initialize list of events used for bulk request
            bulk = []
            # Fill 'bulk' with desired event information
            for k in range(0, nEvents):
                teq = catalog.events[k].origins[0].time
                chan = inv.networks[0].stations[i].channels[j].code
                bulk.append((ntwk, stat, loc, chan, teq, teq + duration*60))

            # Fetch the data!
            if output_units == 'counts':
                try:
                    st = client.get_waveforms_bulk(bulk)
                except Exception as error:
                    with open(logFileName, "a") as log:
                        log.write(str(error))
                        log.write('Unable to complete fetch request for: ' + stat + '.' + loc + '.' + chan)
                    continue
            else:
                try:
                    st = client.get_waveforms_bulk(bulk, attach_response=True)
                except Exception as error:
                    with open(logFileName, "a") as log:
                        log.write(str(error))
                        log.write('Unable to complete fetch request for: ' + stat + '.' + loc + '.' + chan)
                    continue

            # Do some file-formatting and optional minor pre-processing
            for k in range(0, len(st)):
                teq = st[k].meta.starttime

                # Prepare filename for saving
                evchan = st[k].meta.channel
                evid = st[k].meta.starttime.isoformat().replace('-', '.').replace('T', '.').replace(':', '.').split('.')[:-1]
                evid.extend([ntwk, stat, loc, evchan, str(j), 'SAC'])
                evid = ".".join(evid)
                # Add station specific metadata to SAC files
                st[k].stats.sac = {}
                st[k].stats.sac.stla = stla
                st[k].stats.sac.stlo = stlo
                # Channel orientation (CMPAZ)
                azid = [ntwk, stat, loc, evchan]
                azid = ".".join(azid)
                st[k].stats.sac.cmpaz = inv.get_orientation(azid, teq)["azimuth"]

                # Add event-specific metadata to SAC files (surely there must be a faster way to do this...?)
                for l in range(0, nEvents):
                    if catalog.events[l].origins[0].time - 5 <= st[k].meta.starttime <= \
                            catalog.events[l].origins[0].time + 5:
                        st[k].stats.sac.evla = catalog.events[l].origins[0].latitude
                        if st[k].stats.sac.evla is None:
                            with open(logFileName, "a") as log:
                                log.write('Couldn''t find event latitude for: ' + evid + '\n')
                            st[k].stats.sac.evla = 0.0
                        st[k].stats.sac.evlo = catalog.events[l].origins[0].longitude
                        if st[k].stats.sac.evlo is None:
                            with open(logFileName, "a") as log:
                                log.write('Couldn''t find event longitude for: ' + evid + '\n')
                            st[k].stats.sac.evlo = 0.0
                        st[k].stats.sac.evdp = catalog.events[l].origins[0].depth
                        if st[k].stats.sac.evdp is None:
                            with open(logFileName, "a") as log:
                                log.write('Couldn''t find event depth for: ' + evid + '\n')
                            st[k].stats.sac.evdp = 0.0
                        st[k].stats.sac.mag = catalog.events[l].magnitudes[0].mag
                        if st[k].stats.sac.mag is None:
                            with open(logFileName, "a") as log:
                                log.write('Couldn''t find event magnitude for: ' + evid + '\n')
                            st[k].stats.sac.mag = 0.0
                        # Calculate great circle distance and back-azimuth
                        gcarc, baz = su.haversine(stla, stlo, st[k].stats.sac.evla, st[k].stats.sac.evlo)
                        st[k].stats.sac.gcarc = gcarc
                        st[k].stats.sac.baz = baz
                        # Get theoretical P arrival time, and assign to header 'T0'
                        model = TauPyModel(model="iasp91")
                        phases = ["P"]
                        arrivals = model.get_travel_times(source_depth_in_km=st[k].stats.sac.evdp/1000.0,
                                                          distance_in_degree=gcarc, phase_list=phases)
                        st[k].stats.sac.t0 = arrivals[0].time

                        # Save the Pole Zero file index in 'USER0' Header
                        st[k].stats.sac.user0 = j

                        # Save the P-wave ray parameter in 'USER9' Header
                        st[k].stats.sac.user9 = arrivals[0].ray_param*(np.pi/180)

                        # Write the data to a SAC file
                        st[k].write(sac_dir + evid, format='SAC')

    elapsed = time.time() - t1
    print('Time required to complete fetch request: ' + str(elapsed))
    # with open(logFileName, "a") as log:
    #     log.write('Time required to complete fetch request: ' + str(elapsed))



