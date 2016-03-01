************************
Source for msprime paper
************************

This repo contains all source code required to create the msprime paper.

**Note that this is research software that was built for a particular
purpose. It is not high quality software and probably contains bugs.
Use at your own risk!**

The ``src/algorithms.py`` file contains simple implementations of
all the algorithms in the paper.

-------------
Prerequisites
-------------

First, clone the repo::

    $ git clone https://github.com/jeromekelleher/msprime-paper.git

To run the code to compare simulations on a clean Ubuntu 15.10 system,
we need the following::

    $ sudo apt-get install python-dev python-pip libgsl0-dev libhdf5-serial-dev \
        pkg-config python-numpy python-pandas python-matplotlib libboost-dev \
        default-jre asymptote r-cran-ape python-biopython
    $ sudo pip install -r requirements.txt

-----------------------
Paper and illustrations
-----------------------

To create the paper and illustrations, we first need to get the PLoS style
files::

    $ wget http://journals.plos.org/plosone/s/file?id=SyJ3/plos-latex-template.zip
    $ unzip file\?id\=SyJ3%2Fplos-latex-template.zip

The location of this file will probably change, but these commands worked at the
time of writing.

To compile the LaTeX document, run::

    $ make

All of the prerequisites should already be installed, and the data required
to create the plots is already present in the ``data`` directory.


-----------
Simulations
-----------

If you wish to run all the simulations for yourself, follow the three stages
outlined in this section. However, if you just wish to obtain the plots from
data used in the paper, run::

    $ python src/plots.py plot all

All the plots should now be in the figures directory.

+++++++++++++++++++
Running Simulations
+++++++++++++++++++

Creating the plots to compare simulations works in three stages. First, we run the
simulations and output the replicate data; second, we collect the replicate
data and put it into a central file; and finally, we create the plots.

All work running simulations and plotting the results is done in the file
``src/plots.py``.

First, we must download and build the other simulation packages that we'll
be using::

    $ make -C simulators

If everything goes according to plan, we should have working versions of
the simulators we need once this completes.

Before we run any simulations, we must first set the locale to avoid problems
with ``cosi2``::

    $ export LC_ALL=C

To run the simulations, use::

    $ python src/plots.py simulate <dataset ID> -n <num_replicates> -p <num_processes>

The first argument to simulate is a dataset identifier. There are three
datasets, ``alg_stats``, ``tree_simulator_stats`` and
``haplotype_simulator_stats``. There is also a special identifier ``all``, which
will run all available jobs sequentially. Therefore, we might run::

    $ python src/plots.py simulate all -n 1 -p 8

This will run all available simulations with 1 replicate each (so, very inaccurate
timings will be obtained) over 8 threads.

Obtaining accurate timing information from modern processors is quite
challenging, as clock frequencies are dynamically rescaled constantly during
workloads. This can lead to significant skews in timings if, for example, one timing
was obtained when the processor core was running at peak frequency and another
was obtained when it was at minimum frequency. A crude (but effective) way to
work around this is ensure that the CPUs on the system are always busy
when measurements are being taken. For example, if we have a machine with
20 threads then we first set off a stress job at low priority to ensure
that all cores hit their equilibrium frequency::

    $ nice -n 19 stress -c 20

This job will use up any spare CPU resources, and ensure that the clock
frequencies of different cores don't vary unpredictably. If using
an Intel processor it is also a good idea to turn off Turbo Boost
while the measurements are being taken.

++++++++++++++++++
Processing results
++++++++++++++++++

Once a dataset has been simulated, we can then assemble it into a
more convenient form for plotting. This is done by running::

    $ python src/plots.py process <dataset ID>

In the same way as the ``simulate`` command, we provide a dataset ID
as the argument, and we can pass the special ``all`` value to indicate
that we want to process all data sets. After running the ``process``
command the raw replicate data in the ``data/results__NOBACKUP__``
directory is compressed into the corresponding data files in
the data directory. These are stored in JSON format for convenience.

The data used in the paper is included here in the data
directory, and will be over written if you run ``process``.

++++++++++++++++
Plotting results
++++++++++++++++

To plot the results, simply run::

    $ python src/plots.py plot all

Similarly to the other commands, we provide a plot ID to the plot command
to determine which figure we should draw.

--------
Examples
--------

To run the examples in which we compare the sizes of files in various
formats and benchmark different Newick parsers, run::

    $ make -C examples

