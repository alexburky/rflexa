import obspy
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib
import numpy as np
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------------------
# fetchData.py
# This code fetches seismic data for events within a give time window using ObsPy
# and saves the data into .sac files.
# TO DO: Add an 'example' where you fetch a small handful of events and then run them
# through the entire receiver function analysis for reproducibility and teaching purposes
# ------------------------------------------------------------------------------------------
# Last updated 10/31/2019 by aburky@princeton.edu
# ------------------------------------------------------------------------------------------

# Define the client that hosts the desired data
client = Client("IRIS")

t0 = UTCDateTime("2017-01-01")
te = UTCDateTime("2018-10-30")

# Get the coordinates of the station to be used for analysis
inv = client.get_stations(network="IU", station="BBSR", starttime=t0, endtime=te, loc="00", channel="BH*",
                          level="channel")
coords = inv[0].get_coordinates("IU.BBSR.00.BHZ", te)
stla = coords["latitude"]
stlo = coords["longitude"]

# Redefine start time of event search to correspond to start time of seismic station

catalog = client.get_events(starttime=t0, endtime=te, minmagnitude=7, latitude=stla, longitude=stlo,
                            minradius=30, maxradius=90)

# This command extracts the origin time of a given event
# Looks like we'll need to loop over the size of the catalog to get the origins of all the events
t = catalog.events[0].origins[0].time
print(t)

# Fetch a trace!
st = client.get_waveforms("IU", "BBSR", "00", "BH*", t, t + 60*60, attach_response=True)
pre_filt = [0.004, 0.008, 10, 20]
st.remove_response(pre_filt=pre_filt, output="DISP", water_level=70, plot=False, zero_mean=True, taper=True,
                   taper_fraction=0.05)
st.detrend(type='simple')
print(st)

trz = st[2]

# Construct time vector
delta = trz.stats.delta
t = np.arange(0, trz.stats.npts*delta, delta)

# Read in old data to compare
bhz = obspy.read('../data/oldsac/2018*sac')
trz1 = bhz[0]
tz = np.arange(0, trz1.stats.npts*trz1.stats.delta, trz1.stats.delta)

plt.plot(t, trz.data, 'k', linewidth=0.25)
plt.plot(tz-0.75, trz1.data, 'r', linewidth=0.25)
# plt.xlim(350, 700)
plt.show()

# print(catalog.events[0].origins[0].time)
# catalog.plot()

# inventory = client.get_stations(network="IU", station="BBSR", starttime=t0, endtime=te, loc="00", channel="BH*",
#                                level="response")
# print(inventory[0].Channels)
# print(inventory[0].get_contents())

# response = inventory[0].get_response
# print(response)

# response = inventory[0].get_response("IU.BBSR.00.BHZ", te)
# print(response)