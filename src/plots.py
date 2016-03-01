"""
Simulations to produce the plots in the paper.
"""
from __future__ import print_function
from __future__ import division

import os
import sys
import json
import time
import random
import pickle
import argparse
import tempfile
import subprocess
import multiprocessing

import numpy as np
import pandas as pd

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
# Colour wheel from http://tools.medialab.sciences-po.fr/iwanthue/
# Uses the 'intense' preset with k = 5
colour_wheel = [
    "#C252AE",
    "#72953D",
    "#BB5A40",
    "#4E9091",
    "#7D6FA9"
]
font = {'size': 15}

matplotlib.rcParams['axes.color_cycle'] = colour_wheel
matplotlib.rc('font', **font)

from matplotlib import pyplot

import msprime

# Default to human-ish values
DEFAULT_RECOMBINATION_RATE = 1e-8
DEFAULT_EFFECTIVE_POPULATION_SIZE = 1e4
DEFAULT_MUTATION_RATE = 1e-8

TMPFILE_PREFIX = "msprime_sim_"


def harmonic_number(n):
    """
    Returns the nth Harmonic number of the numpy array n
    """
    ret = np.zeros(n.shape[0])
    for j, x in enumerate(n):
        ret[j] = np.sum(1 / np.arange(1, x + 1))
    return ret


def process_worker(work):
    sim, results_dir = work
    sys.stdout.flush()
    datum = sim.run()
    results_file = os.path.join(results_dir, sim.get_results_file())
    with open(results_file, "w") as f:
        pickle.dump(datum, f)


def get_mean_records_per_tree(tree_sequence):
    """
    Returns the mean number of records per tree in the specified tree
    sequence. This is equal to the total number of records divided by
    the number of distinct left-values, or also the total number of
    records divided by the number of trees.
    """
    num_trees = 0
    mean_out = 0
    for _, records_out, records_in in tree_sequence.diffs():
        num_trees += 1
        mean_out += len(records_out)
    return mean_out / num_trees


class Simulator(object):
    """
    Class that represents a simulator for the coalescent with recombination.
    """
    def __init__(self, sample_size, num_loci):
        self.sample_size = sample_size
        self.num_loci = num_loci
        self.timeout = None
        self.recombination_rate = DEFAULT_RECOMBINATION_RATE
        self.effective_population_size = DEFAULT_EFFECTIVE_POPULATION_SIZE
        self.mutation_rate = DEFAULT_MUTATION_RATE
        self.verbosity = 0
        self.num_replicates = 1
        self.simulation_id = 0
        self.generate_trees = True
        self.generate_haplotypes = False

    def set_generate_trees(self, generate_trees):
        self.generate_trees = generate_trees

    def set_generate_haplotypes(self, generate_haplotypes):
        self.generate_haplotypes = generate_haplotypes

    def get_results_file(self):
        return "{0}.pkl".format(self.simulation_id)

    def set_num_replicates(self, num_replicates):
        self.num_replicates = num_replicates

    def set_simulation_id(self, simulation_id):
        self.simulation_id = simulation_id

    def run(self):
        pass

    def get_scaled_recombination_rate(self):
        v = self.recombination_rate * (self.num_loci - 1)
        R = 4 * self.effective_population_size * v
        return R

    def get_scaled_mutation_rate(self):
        v = self.mutation_rate * self.num_loci
        rho = 4 * self.effective_population_size * v
        return rho

    def get_base_datum(self):
        datum = {
            "sample_size": self.sample_size,
            "num_loci": self.num_loci,
            "recombination_rate": self.recombination_rate,
            "R": self.get_scaled_recombination_rate(),
            "simulator": self.simulator_id,
            "num_replicates": self.num_replicates
        }
        return datum


class MsprimeSimulator(Simulator):
    """
    Class representing the msprime simulator. Uses the Python API.
    """
    simulator_id = "msprime"

    def run_replicate(self, j, treefile):
        sim = msprime.TreeSimulator(self.sample_size)
        sim.set_num_loci(self.num_loci)
        scaled_recombination_rate = self.get_scaled_recombination_rate()
        # We don't use ms's scaling.
        scaled_recombination_rate /= (self.num_loci - 1)
        sim.set_scaled_recombination_rate(scaled_recombination_rate)
        sim.set_max_memory("10G")
        sim.run()
        tree_sequence = sim.get_tree_sequence()
        if self.generate_trees:
            tree_sequence.dump(treefile)
        if self.generate_haplotypes:
            scaled_mutation_rate = self.get_scaled_mutation_rate()
            scaled_mutation_rate /= self.num_loci
            tree_sequence.generate_mutations(scaled_mutation_rate)
            with open(treefile, "w") as f:
                for h in tree_sequence.haplotypes():
                    print(h, file=f)
        self.tree_file_size[j] = os.path.getsize(treefile)
        self.used_memory[j] = sim.get_used_memory()
        self.num_trees[j] = sim.get_num_breakpoints()
        self.num_multiple_re_events = \
            sim.get_num_multiple_recombination_events()
        self.num_re_events[j] = sim.get_num_recombination_events()
        self.num_ca_events[j] = sim.get_num_common_ancestor_events()
        self.num_records[j] = tree_sequence.get_num_records()
        self.num_nodes[j] = tree_sequence.get_num_nodes()
        self.num_records_per_tree = get_mean_records_per_tree(tree_sequence)

    def run(self):
        self.tree_file_size = np.zeros(self.num_replicates)
        self.used_memory = np.zeros(self.num_replicates)
        self.num_trees = np.zeros(self.num_replicates)
        self.num_re_events = np.zeros(self.num_replicates)
        self.num_multiple_re_events = np.zeros(self.num_replicates)
        self.num_ca_events = np.zeros(self.num_replicates)
        self.num_records = np.zeros(self.num_replicates)
        self.num_nodes = np.zeros(self.num_replicates)
        self.num_records_per_tree = np.zeros(self.num_replicates)
        with tempfile.NamedTemporaryFile("w", prefix=TMPFILE_PREFIX) as f:
            scaled_recombination_rate = self.get_scaled_recombination_rate()
            scaled_recombination_rate /= (self.num_loci - 1)
            scaled_mutation_rate = self.get_scaled_mutation_rate()
            scaled_mutation_rate /= self.num_loci
            print(
                "running: API", self.sample_size, self.num_loci,
                scaled_recombination_rate,
                scaled_mutation_rate * self.generate_haplotypes)
            before = time.clock()
            for j in range(self.num_replicates):
                self.run_replicate(j, f.name)
            duration = time.clock() - before
        datum = {
            "cpu_time": duration,
            "memory": np.max(self.used_memory),
            "tree_file_size": np.mean(self.tree_file_size),
            "num_trees_mean": np.mean(self.num_trees),
            "num_trees_var": np.var(self.num_trees),
            "num_records_mean": np.mean(self.num_records),
            "num_nodes_mean": np.mean(self.num_nodes),
            "num_records_per_tree_mean": np.mean(self.num_records_per_tree),
            "re_events_mean": np.mean(self.num_re_events),
            "multiple_re_events_mean": np.mean(self.num_multiple_re_events),
            "ca_events_mean": np.mean(self.num_ca_events),
        }
        datum.update(self.get_base_datum())
        return datum


class MsSimulator(Simulator):
    """
    Class representing Hudson's ms simulator.
    """
    simulator_id = "ms"
    executable = ["simulators/ms"]
    replicate_start_token = "//"
    tree_start_token = "["
    segsites_start_token = "segsites"
    # These are Cosi2 specific
    re_events_start_token = "stat nrecomb"
    ca_events_start_token = "stat ncoal"

    def run_command(self, args):
        """
        Runs the specified command returning an iterator over
        the lines in the output file. Stores its user time, system time and
        memory usage statistics.
        """
        with tempfile.TemporaryFile() as stderr:
            print("running:", " ".join(args))
            a = ["/usr/bin/time", "-f%M %S %U"] + args
            bufsize = 32 * 2**10  # Make a large buffer
            p = subprocess.Popen(
                a, bufsize=bufsize, stdout=subprocess.PIPE, stderr=stderr)
            self.tree_file_size = 0
            self.num_trees = np.zeros(self.num_replicates)
            self.segsites = np.zeros(self.num_replicates)
            # These are only used for Cosi2
            self.num_re_events = np.zeros(self.num_replicates)
            self.num_ca_events = np.zeros(self.num_replicates)
            j = -1
            p.poll()
            while p.returncode is None:
                for l in p.stdout:
                    self.tree_file_size += len(l)
                    if l.startswith(self.replicate_start_token):
                        j += 1
                    elif l.startswith(self.tree_start_token):
                        self.num_trees[j] += 1
                    elif l.startswith(self.segsites_start_token):
                        self.segsites[j] = int(l.split()[1])
                    elif l.startswith(self.re_events_start_token):
                        self.num_re_events[j] = int(l.split()[2])
                    elif l.startswith(self.ca_events_start_token):
                        self.num_ca_events[j] = int(l.split()[2])
                p.poll()
            stderr.seek(0)
            if p.returncode != 0:
                raise ValueError("process error: {0}:{1}".format(
                    p.returncode, stderr.read()))
            if j + 1 != self.num_replicates:
                print("ERROR!!!!", j, self.num_replicates, "::", args)
            assert j + 1 == self.num_replicates
            split = stderr.read().split()
            # From the time man page:
            # M: Maximum resident set size of the process during its lifetime,
            #    in Kilobytes.
            # S: Total number of CPU-seconds used by the system on behalf of
            #    the process (in kernel mode), in seconds.
            # U: Total number of CPU-seconds that the process used directly
            #    (in user mode), in seconds.
            self.max_memory = int(split[0]) * 1024
            self.system_time = float(split[1])
            self.user_time = float(split[2])

    def get_args(self):
        rho = self.get_scaled_recombination_rate()
        theta = self.get_scaled_mutation_rate()
        args = self.executable + [
            str(self.sample_size), str(self.num_replicates), "-r", str(rho),
            str(self.num_loci)]
        if self.generate_trees:
            args += ["-T"]
        if self.generate_haplotypes:
            args += ["-t", str(theta)]
        return args

    def _prepare_datum(self):
        datum = {
            "cpu_time": self.system_time + self.user_time,
            "memory": self.max_memory,
            "tree_file_size": self.tree_file_size / self.num_replicates,
            "num_trees_mean": np.mean(self.num_trees),
            "num_trees_var": np.var(self.num_trees),
            "segsites_mean": np.mean(self.segsites),
            "segsites_var": np.var(self.segsites),
            "ca_events_mean": np.mean(self.num_ca_events),
            "ca_events_var": np.var(self.num_ca_events),
            "re_events_mean": np.mean(self.num_re_events),
            "re_events_var": np.var(self.num_re_events),
        }
        datum.update(self.get_base_datum())
        return datum

    def run(self):
        args = self.get_args()
        try:
            self.run_command(args)
        except Exception as e:
            print("EXCEPTION occured:", e)
            print("command = ", " ".join(args))
            raise e
        return self._prepare_datum()


class MscompatSimulator(MsSimulator):
    """
    Class representing the mspms msprime simulator. Assumes mspms is
    available in PATH
    """
    simulator_id = "mspms"
    executable = ["mspms"]

    def get_args(self):
        args = super(MscompatSimulator, self).get_args()
        # Set the max memory to something more sensible
        args += ['-M', '10G']
        return args


class ScrmApproxSimulator(MsSimulator):
    """
    Class representing the SCRM simulator running the SMC' approximation.
    """
    simulator_id = "scrm"
    executable = ["simulators/scrm"]
    def __init__(self, sample_size, num_loci, approx_length=0):
        super(ScrmApproxSimulator, self).__init__(sample_size, num_loci)
        self._approx_length = approx_length

    def get_args(self):
        args = super(ScrmApproxSimulator, self).get_args()
        args += ['-l', str(self._approx_length)]
        return args


class MacsSimulator(MsSimulator):
    """
    Class representing the MaCS simulator.
    """
    simulator_id = "MaCS"
    executable = ["simulators/wrap_macs.sh"]
    replicate_start_token = "Graph Builder begin"
    tree_start_token = "NEWICK_TREE"

    def get_args(self):
        rho = self.get_scaled_recombination_rate()
        theta = self.get_scaled_mutation_rate()
        r = rho / (self.num_loci - 1)
        t = theta / (self.num_loci - 1)
        args = self.executable + [
            str(self.sample_size), str(self.num_loci),
            "-h", "1", "-r", str(r), "-i", str(self.num_replicates)]
        if self.generate_trees:
            args += ["-T"]
        if self.generate_haplotypes:
            args += ["-t", str(t)]
        return args


class MsmsSimulator(MsSimulator):
    """
    Class representing Ewing's msms simulator. Takes care of running the
    simulations and outputing the results.
    """
    executable = ["java", "-Xmx16G", "-jar", "simulators/msms.jar"]
    simulator_id = "msms"


class Cosi2Simulator(MsSimulator):
    """
    Class representing the Cosi2 simulator which simulates the exact
    coalescent but only supports outputting haplotypes.
    """
    simulator_id = "cosi2"

    def run(self):
        assert self.generate_haplotypes
        assert not self.generate_trees
        with tempfile.NamedTemporaryFile(prefix="cosi") as param_file, \
                tempfile.NamedTemporaryFile(prefix="cosi") as recomb_file:
            print("0", DEFAULT_RECOMBINATION_RATE, file=recomb_file)
            print("length", self.num_loci, file=param_file)
            print("mutation_rate", DEFAULT_MUTATION_RATE, file=param_file)
            print("recomb_file", recomb_file.name, file=param_file)
            print("pop_define 1 1", file=param_file)
            print(
                "pop_size 1", DEFAULT_EFFECTIVE_POPULATION_SIZE,
                file=param_file)
            print("sample_size 1", self.sample_size, file=param_file)
            param_file.flush()
            recomb_file.flush()
            args = [
                "simulators/cosi2",  "-p", param_file.name, "-m", "-T",
                "-n", str(self.num_replicates)]
            try:
                self.run_command(args)
            except Exception as e:
                print("EXCEPTION occured:", e)
                print("command = ", " ".join(args))
                param_file.seek(0)
                print("params = \n", param_file.read())
                raise e
        return self._prepare_datum()


class Dataset(object):
    """
    Class representing a dataset backing a plot. A single dataset can
    be used for several plots.
    """
    default_num_replicates = 1000

    def run_simulations(self, num_replicates, num_processes):
        self.make_results_dir()
        for sim_id, s in enumerate(self.simulators):
            s.set_num_replicates(num_replicates)
            s.set_simulation_id(sim_id)

        if num_processes == 1:
            for s in self.simulators:
                process_worker((s, self.get_results_dir()))
        else:
            work = [(s, self.get_results_dir()) for s in self.simulators]
            random.shuffle(work)
            p = multiprocessing.Pool(num_processes)
            p.map(process_worker, work)

    def process_results(self):
        """
        Processes the data in the results dir and creates the data file
        ready for plotting.
        """
        data = []
        results_dir = self.get_results_dir()
        for datafile in os.listdir(results_dir):
            path = os.path.join(results_dir, datafile)
            with open(path) as f:
                data.append(pickle.load(f))
        df = pd.DataFrame(data)
        df.to_json(self.get_data_file())

    def get_data_file(self):
        """
        Returns the data file used for this dataset.
        """
        return "data/{0}.dat".format(self.identifier)

    def read_data_file(self):
        """
        Returns a pandas data frame of the data for this figure.
        """
        # This is a new feature in Pandas - do the input ourselves.
        # df = pd.read_json(self.get_data_file())
        with open(self.get_data_file()) as f:
            d = json.load(f)
        df = pd.DataFrame(d)
        return self.normalise(df)

    def normalise(self, df):
        """
        Normalise any values in the dataframe.
        """
        df.memory /= 1024 * 1024 * 1024
        df.tree_file_size /= 1024 * 1024 * 1024
        df.cpu_time /= df.num_replicates
        return df

    def get_results_dir(self):
        """
        Returns the directory used to keep intermediate results.
        """
        return "data/results__NOBACKUP__/{0}".format(self.identifier)

    def make_results_dir(self):
        """
        Makes a directory to hold the result files for this figure.
        """
        d = self.get_results_dir()

        if not os.path.exists(d):
            os.makedirs(d)

    def scaled_recombination_rate_to_loci(self, R):
        """
        Converts the specified scaled recombination rates to number of loci.
        """
        Ne = DEFAULT_EFFECTIVE_POPULATION_SIZE
        rho = DEFAULT_RECOMBINATION_RATE
        num_loci = (1 + R / (4 * Ne * rho)).astype(int)
        return num_loci


class SmallNDataset(Dataset):
    """
    Dataset holding simulations in which we compare simulators for
    small sample sizes with good approximation values.
    """
    default_num_replicates = 100
    identifier = "small_n"
    fixed_n = 20

    def __init__(self):
        self.simulators = []
        n = self.fixed_n
        num_loci = np.linspace(100, 100 * 10**6, 20).astype(int)
        for m in num_loci:
            self.simulators.append(MscompatSimulator(n, m))
            sim = ScrmApproxSimulator(n, m)
            sim.simulator_id = "scrm_smc"
            self.simulators.append(sim)
            sim = ScrmApproxSimulator(n, m, "100r")
            sim.simulator_id = "scrm_100"
            self.simulators.append(sim)
            sim = ScrmApproxSimulator(n, m, "500r")
            sim.simulator_id = "scrm_500"
            self.simulators.append(sim)


class MixedRandNDataset(Dataset):
    """
    Superclass of datasets in which we obtain data for by fixing
    n and varying m, and vice versa.
    """
    fixed_n = 1000
    fixed_m = 50 * 10**6

    def get_parameters(self):
        num_loci = np.linspace(100, 100 * 10**6, 20).astype(int)
        for m in num_loci:
            yield self.fixed_n, m
        sample_size = np.logspace(1, 5, 20).astype(int)
        for n in sample_size:
            yield n, self.fixed_m


class AlgorithmStatsDataset(MixedRandNDataset):
    """
    Class representing the dataset where we examine fundamental statistics
    of Hudson's algorithm.
    """
    default_num_replicates = 100
    identifier = "alg_stats"

    def __init__(self):
        self.simulators = []
        for n, m in self.get_parameters():
            self.simulators.append(MsprimeSimulator(n, m))


class SimulatorStatsDataset(Dataset):
    """
    Superclass of datasets in which we obtain data for by fixing
    n and varying m, and vice versa.
    """
    fixed_n = 1000
    fixed_m = 50 * 10**6

    def get_parameters(self):
        num_loci = np.linspace(100, 100 * 10**6, 20).astype(int)
        for m in num_loci:
            yield self.fixed_n, m
        sample_size = np.linspace(1000, 10**5, 20).astype(int)
        for n in sample_size:
            yield n, self.fixed_m


class TreeSimulatorStatsDataset(SimulatorStatsDataset):
    """
    Class representing the dataset where we examine collect statistics
    on various coalescent simulators.
    """
    default_num_replicates = 100
    identifier = "tree_simulator_stats"

    def __init__(self):
        self.simulators = []
        for n, m in self.get_parameters():
            self.simulators.append(MsprimeSimulator(n, m))
            if n == self.fixed_n:
                self.simulators.append(ScrmApproxSimulator(n, m))
                self.simulators.append(MscompatSimulator(n, m))
                if m <= 20 * 10**6:
                    self.simulators.append(MsmsSimulator(n, m))
                    self.simulators.append(MacsSimulator(n, m))
            else:
                # We don't use msms at all because it needs too much
                # memory at 50 megabases. We also skip MaCS because it's
                # too slow.
                if n < 10000:
                    self.simulators.append(ScrmApproxSimulator(n, m))
                if n < 50000:
                    self.simulators.append(MscompatSimulator(n, m))


class HaplotypeSimulatorStatsDataset(SimulatorStatsDataset):
    """
    Class representing the dataset where we examine collect statistics
    on various coalescent simulators outputting haplotype data.
    """
    default_num_replicates = 100
    identifier = "haplotype_simulator_stats"

    def __init__(self):
        self.simulators = []
        msprime_simulators = set()
        for n, m in self.get_parameters():
            sim = MsprimeSimulator(n, m)
            self.simulators.append(sim)
            msprime_simulators.add(sim)
            self.simulators.append(MscompatSimulator(n, m))
            if m == self.fixed_m:
                if n < 15000:
                    self.simulators.append(ScrmApproxSimulator(n, m))
                if n < 50000:
                    self.simulators.append(Cosi2Simulator(n, m))
            else:
                self.simulators.append(ScrmApproxSimulator(n, m))
                self.simulators.append(Cosi2Simulator(n, m))
        for sim in self.simulators:
            sim.set_generate_trees(False)
            if sim not in msprime_simulators:
                sim.set_generate_haplotypes(True)


class Figure(object):
    """
    Superclass of all Figures.
    """
    def __init__(self):
        self.dataset = self.dataset_class()
        self.data = self.dataset.read_data_file()

    def save_plot(self, figure):
        """
        Saves the plot to the destination files.
        """
        formats = ["pdf", "eps"]
        files = [
            "figures/{0}.{1}".format(self.identifier, f) for f in formats]
        for f in files:
            figure.savefig(f)
        pyplot.clf()

    def get_fixed_m_data(self):
        df = self.data[self.data.num_loci == self.dataset_class.fixed_m]
        return df.sort("sample_size")

    def get_fixed_n_data(self):
        df = self.data[self.data.sample_size == self.dataset_class.fixed_n]
        df = df[df.num_loci != self.dataset_class.fixed_m]
        return df.sort("R")


class DualPanelFigure(Figure):
    """
    Superclass of figures that plot some value with varying sequence length
    and sample size.
    """
    plot_legend = False
    legend_location = "upper left"
    y_scale_factor = 1
    y_label = ""
    sample_size_scale_factor = 1
    sample_size_label = "Sample size"
    plot_points = False
    plot_lines = True
    specific_labels = {}
    skip_simulators = set()

    def get_colour(self, simulator_id):
        """
        Returns the colour for the specified simulator.
        """
        d = {
            MsprimeSimulator.simulator_id: 0,
            MscompatSimulator.simulator_id: 1,
            MsmsSimulator.simulator_id: 2,
            ScrmApproxSimulator.simulator_id: 3,
            MacsSimulator.simulator_id: 4,
            Cosi2Simulator.simulator_id: 4,
        }
        return colour_wheel[d[simulator_id]]

    def get_marker(self, simulator_id):
        """
        Returns the colour for the specified simulator.
        """
        d = {
            MsprimeSimulator.simulator_id: "^",
            MscompatSimulator.simulator_id: "o",
            MsmsSimulator.simulator_id: "s",
            ScrmApproxSimulator.simulator_id: "+",
            MacsSimulator.simulator_id: "x",
            Cosi2Simulator.simulator_id: "d"
        }
        return d[simulator_id]

    def tweak_plot(self, fig, ax1, ax2):
        pass

    def plot_theoretical_loci_value(self, ax1, line, df, loci_scale):
        pass

    def plot_theoretical_sample_value(self, ax1, line, df, sample_scale):
        pass

    def plot(self):
        fig, (ax1, ax2) = pyplot.subplots(1, 2, sharey=True, figsize=(8, 5.5))
        df = self.get_fixed_n_data()
        R_scale = 10**3
        loci_scale = 10**6
        label_map = {
            "cosi2": "cosi2",
            "scrm": "scrm",
            "MaCS": "MaCS",
            "msms": "msms",
            "mspms": "msprime",
            "msprime": "msprime"
        }
        linestyle = "-" if self.plot_lines else ""
        for sim in ["mspms", "msprime"]:
            if sim in self.specific_labels:
                label_map[sim] = self.specific_labels[sim]

        lines = []
        labels = []
        for sim in sorted(df.simulator.unique()):
            if sim not in self.skip_simulators:
                colour = self.get_colour(sim)
                marker = self.get_marker(sim) if self.plot_points else ""
                s = df[df.simulator == sim]
                s = s.sort("R")
                v = s[self.plotted_value] / self.y_scale_factor
                l, = ax1.plot(
                    s.num_loci / loci_scale, v, color=colour, marker=marker,
                    linestyle=linestyle)
                lines.append(l)
                labels.append(label_map[sim])
        self.plot_theoretical_loci_value(ax1, l, s, loci_scale)

        if self.plot_legend:
            ax1.legend(
                lines, labels, loc=self.legend_location, numpoints=1,
                fontsize="small")
        ax1.set_xlabel("Megabases")
        ax1.set_ylabel(self.y_label)

        df = self.get_fixed_m_data()
        s = df[df.simulator == "msprime"]
        self.plot_theoretical_sample_value(
            ax2, l, s, self.sample_size_scale_factor)
        for sim in df.simulator.unique():
            if sim not in self.skip_simulators:
                colour = self.get_colour(sim)
                marker = self.get_marker(sim) if self.plot_points else ""
                s = df[df.simulator == sim]
                s = s.sort("sample_size")
                v = s[self.plotted_value] / self.y_scale_factor
                ax2.plot(
                    s.sample_size / self.sample_size_scale_factor,
                    v, color=colour, marker=marker,
                    linestyle=linestyle)

        ax2.set_xlabel(self.sample_size_label)
        self.tweak_plot(fig, ax1, ax2)

        ty = ax1.twiny()
        ty.set_xlabel("$\\rho \\times \\, 10^3$")

        def xconv(m):
            Ne = DEFAULT_EFFECTIVE_POPULATION_SIZE
            rho = DEFAULT_RECOMBINATION_RATE
            R = 4 * Ne * rho * (m * loci_scale - 1)
            return R / R_scale
        ty.set_xlim(xconv(x) for x in ax1.get_xlim())

        fig.tight_layout()
        fig.text(0.19, 0.97, "sample size = 1000")
        fig.text(0.60, 0.97, "sequence length = 50Mb")
        self.save_plot(fig)


class SmallNFigure(Figure):
    """
    Figure showing a comparison of scrm and msprime for a small sample size
    and using scrm's improved approximations.
    """
    dataset_class = SmallNDataset

    def plot(self):
        df = self.data
        fig, ax1 = pyplot.subplots()
        R_scale = 10**3
        loci_scale = 10**6
        label_map = {
            "scrm_smc": "scrm -l 0 (SMC')",
            "scrm_100": "scrm -l 100r",
            "scrm_500": "scrm -l 500r",
            "mspms": "msprime (Newick)"
        }
        colour_map = {
            "scrm_smc": colour_wheel[2],
            "scrm_100": colour_wheel[1],
            "scrm_500": colour_wheel[3],
            "mspms": colour_wheel[0],
        }
        marker_map = {
            "mspms": "o",
            "scrm": "+",
            "scrm_smc": "x",
            "scrm_100": "d",
            "scrm_500": "^",
        }

        y_scale_factor = 1
        lines = []
        labels = []
        for sim in ["mspms", "scrm_500", "scrm_100", "scrm_smc"]:
            colour = colour_map[sim]
            marker = marker_map[sim]
            s = df[df.simulator == sim]
            s = s.sort("R")
            v = s[self.plotted_value] / y_scale_factor
            l, = ax1.plot(
                s.num_loci / loci_scale, v, color=colour, marker=marker)
            lines.append(l)
            labels.append(label_map[sim])
        ax1.legend(
            lines, labels, loc="upper left", numpoints=1,
            fontsize="small")
        ax1.set_xlabel("Megabases")
        ax1.set_ylabel(self.y_label)

        ty = ax1.twiny()
        ty.set_xlabel("$\\rho \\times \\, 10^3$")

        def xconv(m):
            Ne = DEFAULT_EFFECTIVE_POPULATION_SIZE
            rho = DEFAULT_RECOMBINATION_RATE
            R = 4 * Ne * rho * (m * loci_scale - 1)
            return R / R_scale
        ty.set_xlim(xconv(x) for x in ax1.get_xlim())
        fig.tight_layout()
        fig.text(0.41, 0.97, "sample size = 20")
        self.save_plot(fig)


class SmallNSimulationTimeFigure(SmallNFigure):
    identifier = "small_n_simulation_time"
    y_label = "CPU Time (seconds)"
    plotted_value = "cpu_time"


class SmallNSimulationMemoryFigure(SmallNFigure):
    identifier = "small_n_simulation_memory"
    y_label = "Memory (Gigabytes)"
    plotted_value = "memory"


class NumEventsFigure(DualPanelFigure):
    """
    The figure plotting the number of coancestry events.
    """
    dataset_class = AlgorithmStatsDataset
    identifier = "num_events"
    plotted_value = "re_events_mean"
    y_scale_factor = 10**6
    y_label = "Recombination events $\\times 10^6$"

    def plot_theoretical_loci_value(self, ax, line, df, loci_scale):
        z = np.polyfit(df.R, df.re_events_mean, 2)
        print("Estimated quadratic = ", z)
        f = np.poly1d(z)
        prediction = f(df.R)
        # Now we need to convert back to loci and scale for plotting.
        y = prediction / self.y_scale_factor
        x = df.num_loci / loci_scale
        ax.plot(x, y, "o", color=line.get_color())

    def tweak_plot(self, fig, ax1, ax2):
        ax1.set_ylim(-0.5, 17)
        ax1.set_xlim(-5, 105)
        ax2.set_xlim(9, 1e5 + 5000)
        ax2.set_xscale("log")


class SimulatorComparisonFigure(DualPanelFigure):
    """
    Superclass of figures used to compare different simulators.
    """
    plot_legend = True
    plot_points = True
    sample_size_scale_factor = 1000
    sample_size_label = "Sample size$\\times 10^3$"


class TreeSimulatorComparisonFigure(SimulatorComparisonFigure):
    """
    Superclass of figures in which we compae tree simulations.
    """
    specific_labels = {
        "mspms": "msprime (Newick)",
        "msprime":  "msprime (HDF5)"
    }
    dataset_class = TreeSimulatorStatsDataset


class TreeSimulationTimeFigure(TreeSimulatorComparisonFigure):
    """
    The figure plotting the running time for various simulators.
    """
    identifier = "tree_simulation_time"
    plotted_value = "cpu_time"
    y_label = "CPU Time (seconds)"
    legend_location = "upper right"

    def tweak_plot(self, fig, ax1, ax2):
        ax1.set_xlim(-5, 105)
        ax1.set_ylim(-5, 250)
        ax2.set_xlim(-5, 105)


class TreeSimulationMemoryFigure(TreeSimulatorComparisonFigure):
    """
    The figure plotting the running memory for various simulators.
    """
    identifier = "tree_simulation_memory"
    plotted_value = "memory"
    y_label = "Memory (Gigabytes)"
    legend_location = "upper right"

    def tweak_plot(self, fig, ax1, ax2):
        ax2.set_ylim(0, 5)


class TreeSimulationNumTreesFigure(TreeSimulatorComparisonFigure):
    """
    The figure plotting the number of trees for various simulators.
    This is a sanity check to make sure we've got the right parameters.
    """
    specific_labels = {}
    skip_simulators = set(["mspms"])
    identifier = "tree_simulation_num_trees"
    plotted_value = "num_trees_mean"
    y_label = "Num trees"
    y_scale_factor = 10**5
    y_label = "Breakpoints $\\times 10^5$"
    plot_points = True
    plot_lines = False

    def plot_theoretical_loci_value(self, ax, line, df, loci_scale):
        R = df.R
        v = R * harmonic_number(df.sample_size - 1)
        # # Now we need to convert back to loci and scale for plotting.
        y = v / self.y_scale_factor
        x = df.num_loci / loci_scale
        ax.plot(x, y, "-k")

    def plot_theoretical_sample_value(self, ax, line, df, sample_scale):
        R = df.R
        v = R * harmonic_number(df.sample_size - 1)
        # # Now we need to convert back to loci and scale for plotting.
        y = v / self.y_scale_factor
        x = df.sample_size / sample_scale
        ax.plot(x, y, "-k")


class HaplotypeSimulatorComparisonFigure(SimulatorComparisonFigure):
    """
    Superclass of figures in which we compare haplotype simulations.
    """
    dataset_class = HaplotypeSimulatorStatsDataset
    specific_labels = {
        "mspms": "msprime (Text)",
        "msprime":  "msprime (HDF5)"
    }


class HaplotypeSimulationTimeFigure(HaplotypeSimulatorComparisonFigure):
    """
    The figure plotting the running time for various simulators.
    """
    identifier = "haplotype_simulation_time"
    plotted_value = "cpu_time"
    y_label = "CPU Time (seconds)"

    def tweak_plot(self, fig, ax1, ax2):
        ax2.set_ylim(0, 400)


class HaplotypeSimulationMemoryFigure(HaplotypeSimulatorComparisonFigure):
    """
    The figure plotting the running memory for various simulators.
    """
    identifier = "haplotype_simulation_memory"
    plotted_value = "memory"
    y_label = "Memory (Gigabytes)"
    plot_points = True

    def tweak_plot(self, fig, ax1, ax2):
        ax2.set_ylim(0, 6)


class HaplotypeSimulationNumEventsFigure(HaplotypeSimulatorComparisonFigure):
    """
    The figure plotting the number of recombination events. This is a sanity
    check for correctness.
    """
    identifier = "haplotype_simulation_num_events"
    plotted_value = "re_events_mean"
    y_label = "RE events $\\times 10^6$"
    plot_points = True
    y_scale_factor = 10**6
    skip_simulators = set(["mspms", "scrm"])


class HaplotypeSimulationSegsitesFigure(HaplotypeSimulatorComparisonFigure):
    """
    The figure plotting the number of segregating sites for various simulators.
    This is a sanity check to make sure we've got the right parameters.
    """
    identifier = "haplotype_simulation_segsites"
    plotted_value = "segsites_mean"
    y_scale_factor = 10**5
    y_label = "Segsites $\\times 10^5$"
    plot_points = True


def run_simulate(cls, args):
    f = cls()
    n = cls.default_num_replicates if args.n == -1 else args.n
    f.run_simulations(n, args.processes)


def run_process(cls, args):
    f = cls()
    f.process_results()


def run_plot(cls, args):
    f = cls()
    f.plot()


def main():
    datasets = [
        AlgorithmStatsDataset,
        TreeSimulatorStatsDataset,
        HaplotypeSimulatorStatsDataset,
        SmallNDataset,
    ]
    plots = [
        NumEventsFigure,
        SmallNSimulationTimeFigure,
        SmallNSimulationMemoryFigure,
        TreeSimulationTimeFigure,
        TreeSimulationMemoryFigure,
        TreeSimulationNumTreesFigure,
        HaplotypeSimulationSegsitesFigure,
        HaplotypeSimulationNumEventsFigure,
        HaplotypeSimulationTimeFigure,
        HaplotypeSimulationMemoryFigure,
    ]
    key_map = dict([(d.identifier, d) for d in datasets + plots])
    parser = argparse.ArgumentParser(
        description=(
            "Run simulations, processes output data and "
            "generates plots."))
    subparsers = parser.add_subparsers()

    simulate_parser = subparsers.add_parser('simulate')
    simulate_parser.add_argument(
        '-n', help="number of replicates", type=int, default=-1)
    simulate_parser.add_argument(
        "--processes", '-p', help="number of processes",
        type=int, default=1)
    simulate_parser.add_argument(
        'key', metavar='KEY', type=str, nargs=1,
        help='the dataset identifier')
    simulate_parser.set_defaults(func=run_simulate)

    process_parser = subparsers.add_parser('process')
    process_parser.add_argument(
        'key', metavar='KEY', type=str, nargs=1,
        help='the simulation identifier')
    process_parser.set_defaults(func=run_process)

    plot_parser = subparsers.add_parser('plot')
    plot_parser.add_argument(
        'key', metavar='KEY', type=str, nargs=1,
        help='the plot dentifier')
    plot_parser.set_defaults(func=run_plot)

    args = parser.parse_args()
    k = args.key[0]
    if k == "all":
        classes = datasets
        if args.func == run_plot:
            classes = plots
        for key, cls in key_map.items():
            if cls in classes:
                print("Running:", key)
                args.func(cls, args)
    else:
        cls = key_map[k]
        args.func(cls, args)

if __name__ == "__main__":
    main()
