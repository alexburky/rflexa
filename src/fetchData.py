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
# Last updated 10/31/2019 by aburky@princeton.edu
# ------------------------------------------------------------------------------------------

# Set matplotlib to use LaTeX
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

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

catalog = client.get_events(starttime=t0, endtime=te, minmagnitude=6, maxmagnitude=7, latitude=stla, longitude=stlo,
                            minradius=30, maxradius=90)

# This command extracts the origin time of a given event
# Looks like we'll need to loop over the size of the catalog to get the origins of all the events
teq = catalog.events[1].origins[0].time

# Fetch a trace!
st = client.get_waveforms("IU", "BBSR", "00", "BH*", teq, teq + 60*60, attach_response=True)
pre_filt = [0.004, 0.008, 10, 20]
st.remove_response(pre_filt=pre_filt, output="DISP", water_level=70, plot=False, zero_mean=True, taper=True,
                   taper_fraction=0.05)
st.detrend(type='simple')
print(st)

tr1 = st[0]
tr2 = st[1]
trz = st[2]

# Construct time vector
delta = trz.stats.delta
t = np.arange(0, trz.stats.npts*delta, delta)

# Figure out how to rotate data!
from obspy.clients.iris import Client
client = Client()
# Get some event parameters for rotation
evla = catalog.events[1].origins[0].latitude
evlo = catalog.events[1].origins[0].longitude
evdp = catalog.events[1].origins[0].depth/1000
print(client.distaz(stla, stlo, evla, evlo))

# Use haversine formula to get back-azimuth
gcarc, baz = su.haversine(stla, stlo, evla, evlo)

# Rotate data!
tr1_azimuth = inv.get_orientation("IU.BBSR.00.BH1", teq)["azimuth"]
ndat, edat = su.seisne(tr1.data, tr2.data, tr1_azimuth)
trr, trt = su.seisrt(ndat, edat, baz)

# Get theoretical P arrival for windowing
model = TauPyModel(model="iasp91")
phases = ["P"]
arrivals = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=gcarc, phase_list=phases)
p = arrivals[0]

# Read in old data to compare
sac_data = obspy.read('../data/oldsac/2018*sac')
bh1 = sac_data[0]
bh2 = sac_data[1]
bhz = sac_data[2]
t_sac = np.arange(0, bhz.stats.npts*bhz.stats.delta, bhz.stats.delta)

# Plot traces
# -----------------------------------------------------------------
f = plt.figure(1, figsize=(12, 10))

p1 = f.add_subplot(311)
p1.plot(t, trr, 'k', linewidth=0.25)
# p1.plot(t_sac-0.75, bh1.data, 'r--', linewidth=0.5)
# p1.set_xlim(350, 700)
p1.set_title('IU.BBSR.00.BHR')

p2 = f.add_subplot(312)
p2.plot(t, trt, 'k', linewidth=0.25)
# p2.plot(t_sac-0.75, bh2.data, 'r--', linewidth=0.5)
# p2.set_xlim(350, 700)
p2.set_title('IU.BBSR.00.BHT')

p3 = f.add_subplot(313)
p3.plot(t, trz.data, 'k', linewidth=0.25)
p3.axvline(p.time, color='r', linewidth=0.25, label='P')
p3.legend()
# p3.plot(t_sac-0.75, bhz.data, 'r--', linewidth=0.5)
# p3.set_xlim(350, 700)
p3.set_title('IU.BBSR.00.BHZ')

f.subplots_adjust(hspace=0.5)
plt.show()
# -----------------------------------------------------------------

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