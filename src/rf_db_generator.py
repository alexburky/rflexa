from obspy.clients.fdsn import Client
from fetchRFdata import fetchRFdata
from computeRFs_new import computeRFs
import multiprocessing
from functools import partial

client = Client("IRIS")
network = "IU"
station = "BBSR"
location = "00"
channel = "BH*"
data_directory = "/Users/aburky/IFILES/NETWORKS/"

# inv = client.get_stations(network=network)

# for stat in inv[0].stations:
#    station = stat.code
#     try:

#fetchRFdata(network=network, station=station, location=location, channel=channel, data_directory=data_directory,
#            minimum_magnitude=5.5, maximum_magnitude=8.0)
if __name__ == '__main__':
    pool = multiprocessing.Pool()
    pool = multiprocessing.Pool(processes=2)
    inputs = [0.0213, 0.0208, 0.0204, 0.0182, 0.0167, 0.0154, 0.0143, 0.0133, 0.0125, 0.0118, 0.0111, 0.0105, 0.01]
    gaussian_width = 1.0
    high_cut = 1.0
    par_computeRFs = partial(computeRFs, network, station, location, data_directory, gaussian_width, high_cut)
    pool.map(par_computeRFs, inputs)
    # p1 = Process(target=computeRFs(network=network, station=station, location=location, data_directory=data_directory,
    #                                gaussian_width=1.0, low_cut=0.0277, high_cut=1.0))
    # p2 = Process(target=computeRFs(network=network, station=station, location=location, data_directory=data_directory,
    #                                gaussian_width=1.0, low_cut=0.0270, high_cut=1.0))
#    p3 = Process(target=computeRFs(network=network, station=station, location=location, data_directory=data_directory,
#                                   gaussian_width=1.0, low_cut=0.0263, high_cut=1.0))
#    p4 = Process(target=computeRFs(network=network, station=station, location=location, data_directory=data_directory,
#                                   gaussian_width=1.0, low_cut=0.0256, high_cut=1.0))
    # p1.start()
    # p2.start()
    # Join them?
    # p1.join()
    # p2.join()
    # p3.join()
    # p4.join()

# computeRFs(network=network, station=station, location=location, data_directory=data_directory, gaussian_width=1.0,
#            low_cut=0.0285, high_cut=1.0)
#    except:
#        print('Data retrieval failed, do some debugging!')

# stations = ["AAE", "BBSR"]
# for stat in stations:
#     try:
#         inv = client.get_stations(network="IU", station=stat, location="00")
#         print(inv)
#     except:
#         print('No data' + stat)
