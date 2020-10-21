from obspy.clients.fdsn import Client
from fetchRFdata import fetchRFdata
from computeRFs_new import computeRFs

client = Client("IRIS")
network = "IU"
# station = "PAYG"
station = "BBSR"
location = "00"
# location = ""
channel = "BH*"
data_directory = "/Users/aburky/IFILES/NETWORKS/"

# stations = ["G01", "G02", "G02B", "G03", "G04", "G04B", "G05", "G06", "G07", "G08", "G09", "G10"]

# stations = ["BIG2", "BYUH", "CCHM", "DLAH", "HPAH", "KCCH", "LHSM", "MRKH", "NGOK", "PHRM"]

# for station in stations:
    # try:
    # computeRFs(network=network, station=station, location=location, data_directory=data_directory,
    #             gaussian_width=1.0)
    # except:
    #     print('Receiver function computation failed, do some debugging!')

fetchRFdata(network=network, station=station, location=location, channel=channel, data_directory=data_directory,
            output_units='velocity', minimum_magnitude=6.9, maximum_magnitude=7.0)

print('Data retrieval successful!')

# print('Data retrieval successful for station: ' + station)

# for station in stations:
computeRFs(network=network, station=station, location=location, data_directory=data_directory, input_units='velocity',
           gaussian_width=1.0)

# Not calculating the SNR in this fetch script currently, need to update this!

