# WaterBot 
a beta version 

Written by Jane Yin
# A set of python scripts that provides a quick testing of the bulk properties for batch of new water models


Waterbot is a small program that allows you to evaluate the bulk properties of new water models. It can assess multiple water models 
in one click, and can compute five critical bulk properties from your MD simulations of water: density, enthalpy of vaporization, dielectric 
constant, isothermal compressibility and thermal expansion coefficient. I am still working on expanding its funtionality, but right now it 
only suports perturbing the van der waals parameters of rigid three site water model (TIP3P).

In order to run Waterbot, you need to first install AMBER (for equilibrating and generating the topology file), OpenMM and MDtraj. 
You may also need a stand-alone version of Parmed, if it is not included in your AmberTools package. The production phase will be performed
using OpenMM and the h5 format trajectories will be generated for analysis using MDtraj program.

The lastest version of OpenMM is 7.1. You can quickly install it using this command:

conda install -c omnia openmm

Note that the 7.1 version requires CUDA 8.0 if you usually run simulations using Nvidia GPUs.

You can install MDtraj by tying the following command:

conda install -c omnia mdtraj



To print the help message of WaterBot, simply issue the command:

python waterbot.py  or python waterbot.py -help
   
To run simulations and compute water properties, you only need to type:

python2 waterbot.py  -i aterbot.in    -p water_param.dat   -o waterbot.out

The "-i" flag indicates the Waterbot input fiel, which stores the user defined options. A sample file waterbot.in is provided here.

The "-p" flag indicates a parameter file, which should contains the radius and epsilon parameters of your water models.
A sample file water_param.dat is provided.

The "-o" flag indidates the name of the output file, which reports the values of five properties mentioned above, as well as a water score.
A sample file janewaterbot.out is provided. The "-o" flag is optional: if the name of the output file is not specified, the default filename
waterbot.out will be used.

Then inside the waterbot input file (waterbot.in):

You need to specify the Amber version you use. if you use amber16, then "Amber16 = yes"; if you still use some older versions of Amber,
such as amber14 or amber12, then use "amber16 = no".

The "quiet" entry is for you to decide whether a detailed input file for each water model will be printed out (quiet = no), which contains
the cumulative results after each iteration, so that you can check the convergence of each property. Otherwise, only the summary file
which contains the final resutls after all iterations will be generated (quiet = yes).

A water score is an arbitrary measurement which reflects the deviation of the computed properties from the experimental data 
for a given water model. A higher score means smaller deviations. The score for TIP3P water is about 16 (sometimes depends on how long
you run). In the sample file waterbot.in, five properties were given an equal weight of 20%. You could put more weight on certain 
properties, or you could put zero weight on proerties that are not siginicant to you.

The geometry and non-bonded parameters of TIP3P water is hand coded in waterbot.py for now, so are the the experimental properties 
of the bulk water in liquid phase at 298.15 K. Please check carefully to make sure those are the numbers you desire. If you want to 
test water properties at different temperatures, you may need to modify the experimental values in the scripts.
If you have any feedback about this program, drop me a few lines to jiyin >_< ucsd.edu. I wil be also glad to answer any questions you have.

Acknowledgements

Many thanks to Niel Henriksen who provided great advice on how to compute enthalpy of vaporization, and David Slochower who provided siginicant help on my Python programming. And our project leader, Michael Gilson, of course. 

