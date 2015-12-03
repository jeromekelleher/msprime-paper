"""
Miscellaneous code to run benchmarks and operate on mpsrime files.
"""
from __future__ import print_function
from __future__ import division

import argparse
import os.path
import tempfile
import time

import msprime
import Bio
from Bio import Phylo
import dendropy
import ete2
import numpy as np


def run_compress(args):
    tree_sequence = msprime.load(args.infile)
    tree_sequence.dump(args.outfile, zlib_compression=True)

def run_mutate(args):
    tree_sequence = msprime.load(args.infile)
    tree_sequence.generate_mutations(args.mutation_rate, random_seed=1)
    tree_sequence.dump(args.outfile)
    print("Generated ", tree_sequence.get_num_mutations(), "mutations")

def run_benchmark_trees(args):
    megabytes = 1024 * 1024
    gigabytes = megabytes * 1024
    terabytes = gigabytes * 1024
    max_trees = args.num_trees
    tree_sequence = msprime.load(args.infile)
    file_size = os.path.getsize(args.infile)
    # Now benchmark the native reading
    before = time.clock()
    total_trees = 0
    for t in tree_sequence.trees():
        total_trees += 1
    duration = time.clock() - before
    print("Read", total_trees, "trees in ", duration, "seconds")

    # Write out the first max_trees in Newick to a temporary file.
    with tempfile.NamedTemporaryFile() as f:
        num_trees = 0
        for _, ns in tree_sequence.newick_trees(10):
            print(ns, file=f)
            num_trees += 1
            if num_trees == max_trees:
                break
        f.flush()
        size = f.tell()
        estimated_size = size * total_trees / max_trees
        print(
            "Estimated total Newick size = ", estimated_size / terabytes, "TB")
        print("Appoximate difference = ", estimated_size / file_size)

        # BioPython
        f.seek(0)
        before = time.clock()
        num_trees = 0
        for tree in Phylo.parse(f, "newick"):
            num_trees += 1
        assert num_trees == max_trees
        avg_duration = (time.clock() - before) / max_trees
        print(
            "Read newick trees in appox {} seconds per tree with BioPython {}".format(
            avg_duration, Bio.__version__))
        days = total_trees * avg_duration / (60 * 60 * 24)
        print("\tEstimated time is {} days".format(days))

        # Dendropy
        f.seek(0)
        before = time.clock()
        num_trees = 0
        for line in f:
            t = dendropy.Tree.get_from_string(line, schema="newick")
            num_trees += 1
        assert num_trees == max_trees
        avg_duration = (time.clock() - before) / max_trees
        print(
            "Read newick trees in appox {} seconds per tree with Dendropy {}".format(
            avg_duration, dendropy.__version__))

        # ETE
        f.seek(0)
        before = time.clock()
        num_trees = 0
        for line in f:
            t = ete2.Tree(line)
            num_trees += 1
        assert num_trees == max_trees
        avg_duration = (time.clock() - before) / max_trees
        print(
            "Read newick trees in appox {} seconds per tree with ETE {}".format(
            avg_duration, ete2.__version__))


def run_write_ped(args):
    tree_sequence = msprime.load(args.infile)
    # Write the ped file.
    s = tree_sequence.get_num_mutations()
    indices = np.arange(0, s)
    filename = "{}.ped".format(args.prefix)
    with open(filename, "w") as f:
        for j, h in enumerate(tree_sequence.haplotypes(), 1):
            family_id = j
            ind_id = j
            paternal_id = 0
            maternal_id = 0
            sex = 1
            # 2 == case, 1 == control
            phenotype = 2 if j <= args.num_cases else 1
            print(
                family_id, ind_id, paternal_id, maternal_id, sex, phenotype,
                end=" ", file=f)
            # Neat trick for converting numerical strings to numpy array:
            # interpret data as raw chars, and subtrace ascii code.
            nph = np.fromstring(h, np.int8) - ord('0')
            # Plink genotypes are 1/2, not 0/1
            nph += 1
            genotypes = np.zeros(4 * s, dtype="S1")
            genotypes[indices * 4] = nph
            genotypes[indices * 4 + 1] = " "
            genotypes[indices * 4 + 2] = nph
            genotypes[indices * 4 + 3] = " "
            genotypes.tofile(f)
            # Remove trailing space
            f.seek(-1, os.SEEK_CUR)
            print(file=f)
    # Make a map file. We're using haploid data, so we put it on the
    # male X chromosome.
    filename = "{}.map".format(args.prefix)
    with open(filename, "w") as f:
        for j in range(tree_sequence.get_num_mutations()):
            print("X", "rs{}".format(j), 0, j + 1, file=f)


def run_gwas(args):
    tree_sequence = msprime.load(args.infile)
    num_cases = args.num_cases
    site = 0
    n = tree_sequence.get_sample_size()
    cases = list(range(1, num_cases + 1))
    with open(args.outfile, "w") as output:
        for tree in tree_sequence.trees(cases):
            for pos, node in tree.mutations():
                num_leaves = tree.get_num_leaves(node)
                cases_with_mut = tree.get_num_tracked_leaves(node)
                controls_with_mut = tree.get_num_leaves(node) - cases_with_mut
                f_cases = cases_with_mut / num_cases
                f_controls = controls_with_mut / (n - num_cases)
                if num_leaves >= n / 2:
                    # The mutation is the major allele
                    a1 = 1
                    a2 = 2
                    fa = 1 - f_cases
                    fu = 1 - f_controls
                else:
                    # The mutation is the minor allele
                    a1 = 2
                    a2 = 1
                    fa = f_cases
                    fu = f_controls
                case_odds = fa / (1 - fa)
                control_odds = fu / (1 - fu)
                odds_ratio = "NA" if control_odds == 0 else case_odds / control_odds
                print(
                    site + 1, a1, fa, fu, a2, odds_ratio, sep="\t",
                    file=output)


def main():
    parser = argparse.ArgumentParser(
        description="Command line interface examples.")
    # This is required to get uniform behaviour in Python2 and Python3
    subparsers = parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    compress_parser = subparsers.add_parser(
        "compress",
        help="Compress the input and write to the output.")
    compress_parser.add_argument("infile")
    compress_parser.add_argument("outfile")
    compress_parser.set_defaults(runner=run_compress)

    mutate_parser = subparsers.add_parser(
        "mutate",
        help="Add mutations")
    mutate_parser.add_argument("infile")
    mutate_parser.add_argument("outfile")
    mutate_parser.add_argument(
        "--mutation-rate", "-u", type=float, default=0,
        help=(
            "The rate at which mutations occur within a single locus "
            "in units of 4N generations"))
    mutate_parser.set_defaults(runner=run_mutate)

    benchmark_trees_parser = subparsers.add_parser(
        "benchmark-trees",
        help="Benchmarks tree processing.")
    benchmark_trees_parser.add_argument("infile")
    benchmark_trees_parser.add_argument(
        "--num-trees", "-n", type=int, default=10)
    benchmark_trees_parser.set_defaults(runner=run_benchmark_trees)

    write_ped_parser = subparsers.add_parser(
        "write-ped",
        help="Converts history file to PED format.")
    write_ped_parser.add_argument("infile")
    write_ped_parser.add_argument("prefix")
    write_ped_parser.add_argument(
        "--num-cases", "-n", type=int, default=10)
    write_ped_parser.set_defaults(runner=run_write_ped)

    gwas_parser = subparsers.add_parser(
        "gwas",
        help="Runs an association test and writes the results to stdout")
    gwas_parser.add_argument("infile")
    gwas_parser.add_argument("outfile")
    gwas_parser.add_argument(
        "--num-cases", "-n", type=int, default=10)
    gwas_parser.set_defaults(runner=run_gwas)

    args = parser.parse_args()
    args.runner(args)


if __name__ == "__main__":
    main()
