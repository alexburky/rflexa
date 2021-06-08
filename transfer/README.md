# rflexa: transfer

This repository contains MATLAB and Python codes for removing the instrument response from seismic data, using SAC_PZ poles and zeros files or RESP files. It also includes links to helpful resources for a more complete introduction to the topic.

## MATLAB Version
The main functionality provided by the MATLAB branch of rflexa: transfer, are three utility functions:
`parsePZ.m`
`parseRESP.m`
`transfer.m`
which can all be found within the `matlab` folder. The first two functions, `parsePZ.m` and `parseRESP.m`, are required to extract the poles, zeros, and constant from SAC_PZ or RESP files. The third function, `transfer.m`, performs the instrument response removal, and is benchmarked against the [SAC transfer function](https://ds.iris.edu/files/sac-manual/commands/transfer.html) of the same name.

In addition to these functions, there are four scripts which can be used to produce the four figures shown in the 2021 EduQuakes paper, "Instrument response removal and the 2020 M3.1 Marlboro, New Jersey, earthquake". The data required to make these figures can be found within the `data` directory, which contains the seismograms and the instrument response information.

## Python Version
The Python branch of rflexa: transfer is intended to mirror the MATLAB branch, and has three utility functions defined within the `transfer.py` code:
`parsePZ`
`parseRESP`
`transfer`
which can all be found within the `python` folder.