from obspy.clients.fdsn import Client

client = Client("IRIS")
print(client)

inv = client.get_stations(network="IU")
print(inv)
