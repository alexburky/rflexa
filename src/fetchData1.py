import obspy
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
import matplotlib
import numpy as np
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import seisutils as su

# ------------------------------------------------------------------------------------------
# fetchData.py
# This code fetches seismic data for events within a give time window using ObsPy
# and saves the data into .sac files.
# TO DO: Add an 'example' where you fetch a small handful of events and then run them
# through the entire receiver function analysis for reproducibility and teaching purposes
# ------------------------------------------------------------------------------------------
# Last updated 11/05/2019 by aburky@princeton.edu
# ------------------------------------------------------------------------------------------

# # Set matplotlib to use LaTeX
# matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# matplotlib.rc('text', usetex=True)

# Define network and station codes to fetch data from
ntwk = "IU"
stat = "BBSR"
loc = "00"
chan = "BH*"

# Define the client that hosts the desired data
client = Client("IRIS")

# Fetch station information for data retrieval
inv = client.get_stations(network=ntwk, station=stat, loc=loc)
nstats = len(inv.networks[0])

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
    catalog = client.get_events(starttime=t0, endtime=tf, minmagnitude=8, maxmagnitude=9, latitude=stla,
                                longitude=stlo, minradius=30, maxradius=90)
    nevents = len(catalog.events)
    # Initialize list of events used for bulk request
    bulk = []
    # Fill 'bulk' with desired event information
    for j in range(0, nevents):
        teq = catalog.events[j].origins[0].time
        bulk.append((ntwk, stat, loc, chan, teq, teq+60*60))

    print(bulk)


# nevents = len(catalog.events)
#
# # Initialize 'bulk' as empty list
# bulk = []
# # Fill 'bulk' with desired event information
# for i in range(0, nevents):
#     teq = catalog.events[i].origins[0].time
#     bulk.append(("IU", "BBSR", "00", "BH*", teq, teq+60*60))
#
# # Fetch data and remove instrument response
# st = client.get_waveforms_bulk(bulk, attach_response=True)
# pre_filt = [0.004, 0.008, 10, 20]
# st.remove_response(pre_filt=pre_filt, output="DISP", water_level=70, plot=False, zero_mean=True, taper=True,
#                    taper_fraction=0.05)
