from obspy.clients.fdsn import Client
from fetchRFdata import fetch_rf_data
from computeRFs import computeRFs

client = Client("IRIS")

network = "TA"
station = "F36A"
location = "ALL"
channel = "BH*"

data_directory = "/Users/aburky/IFILES/NETWORKS/"

fetch_rf_data(network=network, station=station, location=location, channel=channel, data_directory=data_directory,
              output_units='velocity', minimum_magnitude=6.9, maximum_magnitude=7.0)

print('Data retrieval successful!')

# print('Data retrieval successful for station: ' + station)

# for station in stations:
computeRFs(network=network, station=station, location=location, data_directory=data_directory, input_units='velocity',
           gaussian_width=1.0, low_cut=0.02, high_cut=0.2)

# Not calculating the SNR in this fetch script currently, need to update this!

