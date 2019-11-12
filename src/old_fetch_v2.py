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
# ------------------------------------------------------------------------------------------
# Last updated 11/07/2019 by aburky@princeton.edu
# ------------------------------------------------------------------------------------------

# Define network, station, location, and channel codes to fetch data from
ntwk = "IU"
stat = "BBSR"
loc = "00"
chan = "BH*"
# Define the client that hosts the desired data
client = Client("IRIS")
# Define path to directory where seismic data will be saved as SAC files
# sac_dir = "/Users/aburky/PycharmProjects/bermudaRFs/data/rfQuakes/"
sac_dir = "/mnt/usb/aburky/IFILES/STATIONS/IU_BBSR/RFQUAKES/"
if os.path.exists(sac_dir):
    # Maybe add an interface/dialogue that checks with user if they would like to overwrite folder?
    shutil.rmtree(sac_dir)
    os.makedirs(sac_dir)
if not os.path.exists(sac_dir):
    os.makedirs(sac_dir)
# Define amount of data desired (minutes)
duration = 60

# Fetch station information for data retrieval
inv = client.get_stations(network=ntwk, station=stat, loc=loc, channel=chan, level="response")
nstats = len(inv.networks[0])
resp_t0 = []
resp_tf = []
pre_filt = []
for i in range(0, nstats):
    nresp = len(inv.networks[0].stations[i].channels)/3
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
    catalog = client.get_events(starttime=t0, endtime=tf, minmagnitude=5.5, maxmagnitude=7.5, latitude=stla,
                                longitude=stlo, minradius=30, maxradius=90)
    nevents = len(catalog.events)
    # Initialize list of events used for bulk request
    bulk = []
    teq_list = []
    aftershocks = []
    # Fill 'bulk' with desired event information
    for j in range(0, nevents):
        teq_list.append(catalog.events[j].origins[0].time)
        bulk.append((ntwk, stat, loc, chan, teq_list[j], teq_list[j]+duration*60))
        if 0 < j < 250:
            for k in range(1, j):
                t1 = bulk[j][4]
                t2 = bulk[(j-k)][4]
                if t1 <= t2 <= t1+90*60:
                    print('Removing aftershock from request:', t1, t2)
                    aftershocks.append(j)
                    # aftershocks.append(int(j-k))
        if j > 250:
            for k in range(1, 250):
                t1 = bulk[j][4]
                t2 = bulk[(j-k)][4]
                if t1 <= t2 <= t1+90*60:
                    print('Removing aftershock from request:', t1, t2)
                    # aftershocks.append(int(j-k))
                    aftershocks.append(j)
    print(aftershocks)
    ashocks = list(set(aftershocks))
    ashocks.sort(reverse=True)
    print(ashocks)
    # Check 'bulk' for aftershocks/temporally close events and remove from list
    # aftershocks = []
    # for j in range(len(bulk)):
    #    for k in range(j+1, len(bulk)):
    #        t1 = bulk[j][4]
    #        t2 = bulk[k][4]
    #        if t2 <= t1 <= t2+duration*60:
    #            print('Removing aftershock from request:', t1, t2)
    #            aftershocks.append(j)
    for j in range(len(ashocks)):
        del bulk[ashocks[j]]
        del catalog.events[ashocks[j]]
    nevents = len(catalog.events)
    # Fetch the data!
    st = client.get_waveforms_bulk(bulk, attach_response=True)
    l = 0
    for j in range(0, len(st)):
        teq = st[j].meta.starttime
        # Check which instrument response to use for given event
        for k in range(0, len(resp_t0)):
            if resp_t0[k] <= teq <= resp_tf[k]:
                pf = pre_filt[k]
                # print("Filter using:", pf)
        # Remove instrument response
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
        nchans = len(st)/3
        if l < nevents:
            st[j].stats.sac.evla = catalog.events[-(l+1)].origins[0].latitude
            st[j].stats.sac.evlo = catalog.events[-(l+1)].origins[0].longitude
            st[j].stats.sac.evdp = catalog.events[-(l+1)].origins[0].depth
            # Calculate great circle distance and back-azimuth
            gcarc, baz = su.haversine(stla, stlo, st[j].stats.sac.evla, st[j].stats.sac.evlo)
            st[j].stats.sac.gcarc = gcarc
            st[j].stats.sac.baz = baz
            # Get theoretical P arrival time, and assign to header 'T0'
            model = TauPyModel(model="iasp91")
            phases = ["P"]
            arrivals = model.get_travel_times(source_depth_in_km=st[j].stats.sac.evdp/1000.0, distance_in_degree=gcarc,
                                              phase_list=phases)
            st[j].stats.sac.t0 = arrivals[0].time
            print(arrivals[0].time, evid, catalog.events[-(l + 1)].origins[0].time, l)
            # Write the trace to a SAC file
            st[j].write(sac_dir + evid, format='SAC')
            l += 1
        else:
            l = 0
            st[j].stats.sac.evla = catalog.events[-(l+1)].origins[0].latitude
            st[j].stats.sac.evlo = catalog.events[-(l+1)].origins[0].longitude
            st[j].stats.sac.evdp = catalog.events[-(l+1)].origins[0].depth
            st[j].stats.sac.mag = catalog.events[-(l+1)].magnitudes[0].mag
            # Calculate great circle distance and back-azimuth
            gcarc, baz = su.haversine(stla, stlo, st[j].stats.sac.evla, st[j].stats.sac.evlo)
            st[j].stats.sac.gcarc = gcarc
            st[j].stats.sac.baz = baz
            # Get theoretical P arrival time, and assign to header 'T0'
            model = TauPyModel(model="iasp91")
            phases = ["P"]
            arrivals = model.get_travel_times(source_depth_in_km=st[j].stats.sac.evdp/1000.0, distance_in_degree=gcarc,
                                              phase_list=phases)
            st[j].stats.sac.t0 = arrivals[0].time
            print(arrivals[0].time, evid, catalog.events[-(l+1)].origins[0].time, l)
            # Write the trace to a SAC file
            st[j].write(sac_dir + evid, format='SAC')
            l += 1
