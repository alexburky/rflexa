from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
import numpy as np
import seisutils as su
import os
import shutil

# ------------------------------------------------------------------------------------------
# fetchData.py
# This script fetches seismic data for events within a give time window using ObsPy
# and saves the data into .sac files.
# TO DO: Add an 'example' where you fetch a small handful of events and then run them
# through the entire receiver function analysis for reproducibility and teaching purposes
#
# - Add an option that allows the user to not remove the instrument response
#
# ------------------------------------------------------------------------------------------
# Last updated 8/18/2020 by aburky@princeton.edu
# ------------------------------------------------------------------------------------------


# Make the whole thing a function!
def fetchRFdata(network, station, location, channel, data_directory, minimum_magnitude=6.9, maximum_magnitude=7.0,
                remove_response=True):
    # Define network, station, location, and channel codes to fetch data from
    ntwk = network
    stat = station
    loc = location
    chan = channel

    # Define the client that hosts the desired data
    client = Client("IRIS")
    # Define path to directory where seismic data will be saved as SAC files
    # sac_dir = "/mnt/usb/aburky/IFILES/STATIONS/" + ntwk + "_" + stat + "/RFQUAKES/"
    if remove_response:
        sac_dir = data_directory + ntwk + "/" + stat + "/" + loc + "/RFQUAKES/"
    else:
        sac_dir = data_directory + ntwk + "/" + stat + "/" + loc + "/RFQUAKES_COUNTS/"

    if os.path.exists(sac_dir):
        overwrite = input('Earthquake directory exists! Would you like to overwite it? [y/n]')
        if overwrite == 'y':
            shutil.rmtree(sac_dir)
        elif overwrite == 'n':
            print('You have chosen not to overwrite the directory. Terminating process...')
            quit()
        else:
            print('ERROR: Invalid response. Terminating process...')
            quit()
    if not os.path.exists(sac_dir):
        os.makedirs(sac_dir)

    # Define amount of data desired (minutes)
    duration = 60

    # Fetch station information for data retrieval
    if loc == "":
        inv = client.get_stations(network=ntwk, station=stat, channel=chan, level="response")
    else:
        inv = client.get_stations(network=ntwk, station=stat, loc=loc, channel=chan, level="response")

    nstats = len(inv.networks[0])
    resp_t0 = []
    resp_tf = []
    pre_filt = []
    for i in range(0, nstats):
        nresp = int(len(inv.networks[0].stations[i].channels)/3)
        for j in range(0, nresp):
            # Start time of station operation for given channels
            resp_t0.append(inv.networks[0].stations[i].channels[j].start_date)
            # End time of station operation for given channels
            resp_tf.append(inv.networks[0].stations[i].channels[j].end_date)
            fs = inv.networks[0].stations[i].channels[j].sample_rate
            # Get the instrument response and corresponding frequencies
            resp, freq = inv.networks[0].stations[i].channels[j].response.get_evalresp_response(1.0/fs, round(fs/1e-5),
                                                                                                output='VEL')
            # Find frequencies where instrument response is 'flat'
            dresp = np.diff(np.log10(abs(resp)))
            test = np.isclose(0, dresp, atol=3e-4)
            idx = np.where(test == True)
            idx = idx[0][0]
            f1 = freq[idx]
            f2 = f1*2.0
            f3 = fs/2.0
            pre_filt.append((f1, f2, f3, fs))

    # Loop over time-periods during which station was operational
    for i in range(0, nstats):
        if i == nstats - 1:
            t0 = inv.networks[0].stations[i].start_date
            tf = UTCDateTime.now()
        else:
            t0 = inv.networks[0].stations[i].start_date
            tf = inv.networks[0].stations[i].end_date
        # Get station coordinates for event selection
        stla = inv.networks[0].stations[i].latitude
        stlo = inv.networks[0].stations[i].longitude
        # Fetch relevant events in time-window during which station was operational
        catalog = client.get_events(starttime=t0, endtime=tf, minmagnitude=minimum_magnitude,
                                    maxmagnitude=maximum_magnitude, latitude=stla,
                                    longitude=stlo, minradius=30, maxradius=90)
        nevents = len(catalog.events)
        # Initialize list of events used for bulk request
        bulk = []
        # Fill 'bulk' with desired event information
        for j in range(0, nevents):
            teq = catalog.events[j].origins[0].time
            bulk.append((ntwk, stat, loc, chan, teq, teq+duration*60))
        # Fetch the data!
        st = client.get_waveforms_bulk(bulk, attach_response=True)
        # Do some minor pre-processing and file-formatting
        for j in range(0, len(st)):
            teq = st[j].meta.starttime
            # Check which instrument response to use for given event
            for k in range(0, len(resp_t0)):
                if resp_t0[k] <= teq <= resp_tf[k]:
                    pf = pre_filt[k]
            # Remove instrument response
            if remove_response:
                st[j].remove_response(pre_filt=pf, output="DISP", water_level=70, zero_mean=True, taper=True,
                                      taper_fraction=0.05)
            # Prepare filename for saving
            evchan = st[j].meta.channel
            evid = st[j].meta.starttime.isoformat().replace('-', '.').replace('T', '.').replace(':', '.').split('.')[:-1]
            evid.extend([ntwk, stat, loc, evchan, 'sac'])
            evid = ".".join(evid)
            # Add station-specific metadata to SAC files
            st[j].stats.sac = {}
            st[j].stats.sac.stla = stla
            st[j].stats.sac.stlo = stlo
            # Channel orientation (CMPAZ)
            azid = [ntwk, stat, loc, evchan]
            azid = ".".join(azid)
            st[j].stats.sac.cmpaz = inv.get_orientation(azid, teq)["azimuth"]
            # Add event-specific metadata to SAC files
            for k in range(0, nevents):
                if catalog.events[k].origins[0].time - 5 <= st[j].meta.starttime <= catalog.events[k].origins[0].time + 5:
                    print('Match!', catalog.events[k].origins[0].time, st[j].meta.starttime)
                    st[j].stats.sac.evla = catalog.events[k].origins[0].latitude
                    if st[j].stats.sac.evla is None:
                        print('Couldn''t find event latitude, setting to zero.')
                        st[j].stats.sac.evla = 0.0
                    st[j].stats.sac.evlo = catalog.events[k].origins[0].longitude
                    if st[j].stats.sac.evlo is None:
                        print('Couldn''t find event longitude, setting to zero.')
                        st[j].stats.sac.evlo = 0.0
                    st[j].stats.sac.evdp = catalog.events[k].origins[0].depth
                    if st[j].stats.sac.evdp is None:
                        print('Couldn''t find event depth, setting to zero.')
                        st[j].stats.sac.evdp = 0.0
                    st[j].stats.sac.mag = catalog.events[k].magnitudes[0].mag
                    if st[j].stats.sac.mag is None:
                        print('Couldn''t find event magnitude, setting to zero.')
                        st[j].stats.sac.mag = 0.0
                    # Calculate great circle distance and back-azimuth
                    gcarc, baz = su.haversine(stla, stlo, st[j].stats.sac.evla, st[j].stats.sac.evlo)
                    st[j].stats.sac.gcarc = gcarc
                    st[j].stats.sac.baz = baz
                    # Get theoretical P arrival time, and assign to header 'T0'
                    model = TauPyModel(model="iasp91")
                    phases = ["P"]
                    arrivals = model.get_travel_times(source_depth_in_km=st[j].stats.sac.evdp/1000.0,
                                                      distance_in_degree=gcarc, phase_list=phases)
                    st[j].stats.sac.t0 = arrivals[0].time
                    st[j].stats.sac.user9 = arrivals[0].ray_param*(np.pi/180)
                    print(arrivals[0].time, evid, catalog.events[k].origins[0].time)
                    # Write the trace to a SAC file
                    st[j].write(sac_dir + evid, format='SAC')

    print('Data fetching complete!')
