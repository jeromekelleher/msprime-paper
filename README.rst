************************
Source for msprime paper
************************

This repo contains all source code required to create the msprime paper.

**Note that this is research software that was built for a particular
purpose. It is not high quality software and probably contains bugs.
Use at your own risk!**

+++++++++++++
Prerequisites
+++++++++++++

To run the code to compare simulations on a clean Ubuntu 15.10 system,
we need the following::

    $ sudo apt-get install python-dev python-pip libgsl0-dev libhdf5-serial-dev \
        pkg-config python-numpy python-pandas python-matplotlib libboost-dev \
        default-jre
    $ pip install -r requirements.txt --user


+++++++++++++++++++
Running Simulations
+++++++++++++++++++

Creating the plots to compare simulations works in three stages. First, we run the
simulations and output the replicate data; second, we collect the collect the replicate
data and put it into a central file; and finally, we create the plots.

All work running simulations and plotting the results is done in the file
``src/plots.py``.

First, we must download and build the other simulation packages that we'll
be using::

    $ make -C simulators

If everything goes according to plan, we should have working versions of
the simulators we need once this completes.

To run the simulations, use::

    $ python src/plots.py simulate

