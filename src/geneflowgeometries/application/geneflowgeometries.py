#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
## Copyright (c) 2023 Adrian Ortiz-Velez.
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * The names of its contributors may not be used to endorse or promote
##       products derived from this software without specific prior written
##       permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL ADRIAN ORTIZ-VELEZ BE LIABLE FOR ANY DIRECT,
## INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
## BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
## LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
## OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
## ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

# arun code comments and readme
# kelley gihub

# get set restriction

# json config file********


# OOP class methods
# pass config and output writers******

# log handlers*

# runtime context
# loger
# output handler
#

# demes_sequences class
# simulation a method for demes
# or
# simulation a class that is used as inhertiance and acts on demes seq


import os
import pathlib  # barcharts next
import sys
import argparse
import datetime
import random


from geneflowgeometries.simulate import (
    ancestral_deme_sequences,
    continuous_trait_evolution,
)
from geneflowgeometries.calculate import analyze


class InvalidMigrationRateException(Exception):
    "Raised when the input value is greater than 1/k, where k is the number of demes"
    pass


def main():
    parser = argparse.ArgumentParser(description=None)

    parser.add_argument(
        "-S",
        "--simulate-sequences-or-continous",
        action="store",
        default="sequences",
        help="Choose to simulate sequences or continous variable [default=%(default)s].",
    )
    parser.add_argument(
        "-G",
        "--geometry",
        action="store",
        default="complete graph",
        help="Network/graph to simulate. Options: Complete, chain, or FILE [default=%(default)s].",
    )
    parser.add_argument(
        "-N",
        "--number-of-chromosomes-per-deme",
        action="store",
        default=100,
        help="Number of chromosomes per deme [default=%(default)s].",
    )
    parser.add_argument(
        "-P",
        "--number-of-chromosomes-per-invdividual",
        action="store",
        default=2,
        help="Number of ploidy by organism [default=%(default)s].",
    )
    parser.add_argument(
        "-k",
        "--number-of-demes",
        action="store",
        default=4,
        help="Number of demes in simulation[default=%(default)s].",
    )
    parser.add_argument(
        "-m",
        "--migration-rate",
        action="store",
        default=0,
        help="Migration rate across demes. Max = 1/k, where k is number of demes [default=%(default)s].",
    )
    parser.add_argument(
        "-mut",
        "--mutation-rate",
        action="store",
        default=0.0004,
        help="Mutation rate for gene simulaiton [default=%(default)s].",
    )
    parser.add_argument(
        "-L",
        "--sequence-length",
        action="store",
        default=1000,
        help="Sequence length for gene simulaiton [default=%(default)s].",
    )
    parser.add_argument(
        "-mean",
        "--mean",
        action="store",
        default=0,
        help="Starting mean for simulating continous varible[default=%(default)s].",
    )
    parser.add_argument(
        "-std",
        "--standard-deviation",
        action="store",
        default=1,
        help="Starting standard deviation for simulating continous varible[default=%(default)s].",
    )
    parser.add_argument(
        "-simT",
        "--simulation-time",
        action="store",
        default=2000,
        help="Number of generations to run simulation for [default=%(default)s].",
    )
    parser.add_argument(
        "-simK",
        "--number-of-simulations",
        action="store",
        default=200,
        help="Number of trials of simulations to run [default=%(default)s].",
    )
    parser.add_argument(
        "-id",
        "--log-id",
        action="store",
        default="log-DATE",
        help="Filename for log files [default=%(default)s].",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        action="store",
        default="output",
        help="Prefix for output files [default=%(default)s].",
    )
    parser.add_argument(
        "-s",
        "--seed",
        action="store",
        default="random",
        help="Seed for reproducibilty [default=%(default)s].",
    )
    parser.add_argument(
        "-cf",
        "--config",
        action="store",
        default="FILE",
        help="Config file [default=%(default)s].",
    )


###########################LOGGER
    args = parser.parse_args()
    print("Hello, simulation begining")
    now = str(datetime.datetime.now())

    print(now.split(" ")[0])
    date = now.split(" ")[0]
    main_log_row = now.split(" ")[0] + ","
###############################



#**********************************CONFIG

    import json
    #Load the configuration file 
    with open(args.config) as f:
        config = f.read()
    


    config_dict = json.loads(config)
    print(config_dict)
    print(config_dict.keys())




    # Pass args
    simulate_what = config_dict['simulator']['simulate_sequences_or_continous']
    print("Simulating", simulate_what)
    Geometry = args.geometry
    number_of_chromosomes = int(args.number_of_chromosomes_per_deme)
    number_of_ploidy = int(args.number_of_chromosomes_per_invdividual)
    number_of_demes = int(args.number_of_demes)
    migration_rate = float(args.migration_rate)

    try:
        if migration_rate > 1 / number_of_demes:
            raise InvalidMigrationRateException
        else:
            print("Valid migration rate")
    except InvalidMigrationRateException:
        print("Exception occurred: Invalid migration rate")
        raise SystemExit(1)

    mutation_rate = float(args.mutation_rate)
    sequence_length = int(args.sequence_length)
    simT = int(args.simulation_time)
    start_mean = float(args.mean)
    start_std = float(args.standard_deviation)
    snapshot_times = [
        1,
        number_of_chromosomes,
        5 * number_of_chromosomes,
        10 * number_of_chromosomes,
        20 * number_of_chromosomes,
    ]
    # random number
    if args.seed == "random":
        seed_range = 2**32 - 1
        seed = int(random.uniform(0, seed_range))
    else:
        seed = int(args.seed)
#*********************************   
        
        
        
#############################
    main_log_row += simulate_what + ","
    main_log_row += Geometry + ","
    main_log_row += str(number_of_chromosomes) + ","
    main_log_row += str(number_of_ploidy) + ","
    main_log_row += str(number_of_demes) + ","
    main_log_row += str(migration_rate) + ","
    main_log_row += str(mutation_rate) + ","
    main_log_row += str(sequence_length) + ","
    main_log_row += str(simT) + ","
    main_log_row += str(snapshot_times) + ","
######################



    # make config class dictionary
    if simulate_what == "sequences":
        (filenames_snapshot, tag) = ancestral_deme_sequences.ancestral_deme_sequences(
            config_dict)

        for i, fastafiles in enumerate(filenames_snapshot[0]):
            filenames = [filenames_snapshot[0][i], filenames_snapshot[1][i]]
            print("Analysis for ", filenames)

            (sequence_dataframe, data_matrix) = analyze.parse(filenames)

            analyze.nei_fst(sequence_dataframe, data_matrix, tag)

            analyze.weir_goudet_population_specific_fst(
                sequence_dataframe, data_matrix, tag
            )

            analyze.by_deme_pairwise_fst(sequence_dataframe, data_matrix, tag)

            analyze.wright_fis(sequence_dataframe, data_matrix, tag)

            # analyze.pairwise_fst(sequence_dataframe, data_matrix,tag)
    else:
        print(seed)
        csvfile = continuous_trait_evolution.continuous_trait_evoluion(config_dict)

    f = open("main.log", "a")
    f.write(main_log_row)
    f.close()


if __name__ == "__main__":
    main()
