# rflexa

rflexa is a complete receiver function workflow tool,
with the purpose of assisting the observational seismologist *from
station to subsurface*. Although initially built for receiver function
analysis, the individual modules can be extracted for use and
integrated into your existing workflow. rflexa exists as two independent
branches, a MATLAB version, and a Python version.

## Python Version
This guide is not intended to be a complete introduction to Python,
but if you are new to the Python programming language I highly
recommend you go through the following steps before continuing any
further with rflexa:

1. Install the [Anaconda](https://www.anaconda.com/distribution/) package manager.  

2. Depending on how you prefer to write code, do one of the following:  

* If you like to work from the command-line/terminal, read about setting
up
[conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html)
to work from the command-line.

* If you like to work within an IDE, search Google for a tutorial on how
to setup conda in your respective IDE. I recommend
[PyCharm](https://www.jetbrains.com/pycharm/download/#section=mac), as
it is built entirely for Python development and comes with a thorough
[conda setup
guide](https://www.jetbrains.com/help/pycharm/conda-support-creating-conda-virtual-environment.html). The
Community version is free and has been sufficient for my needs.

rflexa is built on top of the
[ObsPy](https://github.com/obspy/obspy/wiki/Installation-via-Anaconda)
package, and requires that you have ObsPy installed in order to run
properly. If you already have ObsPy installed, you can continue to the
download and use the source code from the folders above.

## MATLAB Version
For a thorough explanation of the functions included in the MATLAB version of rflexa,
consult the documentation [here](https://github.com/alexburky/rflexa).