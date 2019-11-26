from obspy.clients.fdsn import Client
from fetchRFdata import fetchRFdata
from computeRFs_new import computeRFs

client = Client("IRIS")
network = "IU"
location = "00"
channel = "BH*"
data_directory = "/mnt/usb/aburky/IFILES/NETWORKS/"

# inv = client.get_stations(network=network)

# for stat in inv[0].stations:
#    station = stat.code
#     try:
station = "BBSR"
#fetchRFdata(network=network, station=station, location=location, channel=channel, data_directory=data_directory,
#            minimum_magnitude=5.5, maximum_magnitude=8.0)
computeRFs(network=network, station=station, location=location, data_directory=data_directory, gaussian_width=1.0,
           low_cut=0.1, high_cut=1)
#    except:
#        print('Data retrieval failed, do some debugging!')

# stations = ["AAE", "BBSR"]
# for stat in stations:
#     try:
#         inv = client.get_stations(network="IU", station=stat, location="00")
#         print(inv)
#     except:
#         print('No data' + stat)
