# rflexa

rflexa is a complete receiver function workflow tool,
with the purpose of assisting the observational seismologist *from
station to subsurface*. Although initially built for receiver function
analysis, the individual modules can be extracted for use and
integrated into your existing workflow. rflexa exists as two independent
branches, a MATLAB version, and a Python version.

[Installing the TauP Toolkit](###the-taup-toolkit)

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

## The TauP Toolkit
The TauP Toolkit is a useful command line utility which allows you to quickly calculate theoretical
travel times for a wide variety of seismic phases. To begin installing it, go to
the [TauP Website](https://www.seis.sc.edu/taup/) and download the latest release
(as of this writing, TauP-2.4.5).

If you have experience with the command line (also known as the Terminal), then you'll most likely
be able to complete the installation by following Section A of the official TauP documentation. If you're
not very comfortable with the command line, then continue reading!

### MacOS/Linux Users

#### Step 1: Download the TauP Toolkit
If you haven't already done so, download the latest version of the [TauP Toolkit](https://www.seis.sc.edu/taup/) and save it to your `~/Downloads` folder.

#### Step 2: Open a Terminal Window
Once the Terminal successfully launches, execute the following command to navigate to your `~/Downloads` folder:
```bash
cd ~/Downloads
```
you will make frequent use of the Terminal as a seismologist, so it's best if you start
getting familiar with it now! The `cd` command means *change directory.* In Unix world,
a folder is referred to as a *directory*, so this command is equivalent to opening up
Finder and clicking on your **Downloads** folder.

#### Step 3: Unzip your TauP Installation
Now that we are in the `~/Downloads` folder, go ahead and run the following command in the Terminal:
```bash
gunzip TauP*
```
if this fails for some reason, try running the following commands, one at a time, until one works:
```bash
tar -xvf TauP*
unzip TauP*
```
You should now have a new folder in your `~/Downloads` folder called `TauP-X.X.X` where `X` refers to
the version number. Before we contine, lets unpack the commands we just executed a bit.

The commands `gunzip`, `tar -xvf`, and `unzip` all have the effect of taking a `.tar` or `.zip` file and
inflating them. When you are sharing a large piece of software, or any large folder, it is good practice
to *compress* the folder so that the file size isn't quite as large. Your installation of the TauP Toolkit
came as a compressed file, so we needed to decompress it by running the above command.

The second part of our command, referred to as the *argument*, was `TauP*`. The `*` in this command is known
as a *wildcard*. This means that, after executing our command, the computer checks to see if there are any files
that start with `TauP` in the current directory. If it finds a match, it executes the command on this file, and
then proceeds to check if there are any more files that match. If, for example, you had downloaded 5 copies of
the TauP Toolkit, this command would effectively decompress *all 5 of them*.

#### Step 4: Move TauP to your Home Directory
Next, we are going to move the TauP folder somewhere where we can easily find and access it. Since some of you
may be working on computers without root access, we will put it in our home directory. You can do this by
executing:
```bash
mv TauP* ~/
```
to check that this was successful, execute:
```bash
cd ~/
```
followed by:
```bash
ls
```
which *lists* all of the files and folders in the current directory. The first command, `mv`, had the effect
of moving the TauP folder to the home directory `~/`. The second then took us from the `~/Downloads` folder to
our home directory, which is universally referred to as `~/`. Before moving forward, make sure the `ls` command
showed the TauP folder.

#### Step 5: Add TauP to your $PATH
Now, we need to add TauP to our $PATH variable. You can do this by opening your `~/.bash_profile` file
using a command line editor such as [emacs](https://www.gnu.org/software/emacs/tour/). If you have emacs
up and running, the following will open your `~/.bash_profile` file:
```bash
emacs -nw ~/.bash_profile
```
once the file opens, add the following lines to the bottom:
```bash
# Added for the TauP Toolkit
export TAUPHOME=/Users/$USERNAME/Taup-X.X.X
export PATH=$TAUPHOME/bin:$PATH
```
where you replace `$USERNAME` by your username, which you can find by typing `pwd` in your home directory, and Taup-X.X.X with the version of TauP that you downloaded. These commands are creating a new *environment variable* named
`TAUPHOME`, and adding it to your `$PATH`. In Unix, your `$PATH` is a set of directories which your command line
can access from anywhere, so files in those directories, like the TauP executables, can be executed from anywhere.

#### Step 6: Source your ~/.bash_profile
The last step before we can run the TauP Toolkit is applying the changes we just made to our `~/.bash_profile`
file by *sourcing* it with the following:
```bash
source ~/.bash_profile
```
After doing this, verify that everything has worked by executing:
```bash
taup_time
```
if everything was successful, an interactive input should appear asking you to specify model parameters for
a travel time calculation. Congratulations! You now have a useful command line utility ready to invoke at a
moment's notice.
