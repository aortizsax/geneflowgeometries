#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
## Copyright (c) 2023 Adrian Ortiz-Velez and Jeet Sukumaran.
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
## DISCLAIMED. IN NO EVENT SHALL ADRIAN ORTIZ-VELEZ OR JEET SUKUMARAN BE LIABLE
## FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
## OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################
# mutate funciton not nnlog(n)
# actully knml only for loo thourgh mutation sites

# directive from last week get ancsequences to OOP
# additionally continuous trait,
# argparse
#   no more casting
#   override config file
#


# my to do
# modularization of writing files(close, have note),
# analysis OOP similar to PDM class??
# log files;
#   can write and choose level
#   need ';' terrminator and multilog handling
#


import os
import pathlib  # barcharts next
import sys
import argparse
import datetime
import random
import numpy as np

from geneflowgeometries.simulate.Simulator import Simulator
from geneflowgeometries.calculate import analyze

import logging
from geneflowgeometries.log.multi_handlers_logger import setup_logger
from geneflowgeometries.config_parse.Config import read_config, read_args

from geneflowgeometries.config_parse.Config import Config
from configparser import ConfigParser


def main():
    args = sys.argv[1:]
    # setting the log level on the root logger must happen BEFORE any output
    logging_argparse = argparse.ArgumentParser(prog=__file__, add_help=False)
    logging_argparse.add_argument(
        "-l", "--log-level", default="WARNING", help="set log level"
    )
    logging_args, _ = logging_argparse.parse_known_args(args)

    try:
        logging.basicConfig(level=logging_args.log_level)
    except ValueError:
        logging.error("Invalid log level: {}".format(logging_args.log_level))
        sys.exit(1)

    logger = logging.getLogger(__name__)
    logger.info(
        "Log level set: {}".format(logging.getLevelName(logger.getEffectiveLevel()))
    )

    # parse values from a configuration file if provided and use those as the
    # default values for the argparse arguments
    config_argparse = argparse.ArgumentParser(prog=__file__, add_help=False)
    config_argparse.add_argument(
        "-c", "--config-file", help="path to configuration file"
    )
    config_args, _ = config_argparse.parse_known_args(args)

    defaults = {
        "simulate_sequences_or_continous": "sequences",
        "geometry": "complete graph",
        "number_of_chromosomes_per_deme": 100,
        "number_of_chromosomes_per_invdividual": 2,
        "number_of_demes": 4,
        "migration_rate": 0.1,
        "mutation_rate": 0.0004,
        "sequence_length": 1000,
        "simulation_time": 2000,
        "number_of_simulations": 200,
        "log_id": "log_DATE",
        "output_prefix": "output",
        "seed": 12345,
        "mean": 0,
        "standard_deviation": 1,
    }

    if config_args.config_file:
        logger.info("Loading configuration: {}".format(config_args.config_file))
        try:
            config_parser = ConfigParser()
            with open(config_args.config_file) as f:
                config_parser.read_file(f)
            config_parser.read(config_args.config_file)
        except OSError as err:
            logger.error(str(err))
            sys.exit(1)

        defaults.update(dict(config_parser.items("SIMULATOR")))

    # parse the program's main arguments using the dictionary of defaults and
    # the previous parsers as "parent' parsers
    parsers = [logging_argparse, config_argparse]
    main_parser = argparse.ArgumentParser(prog=__file__, parents=parsers)
    main_parser.set_defaults(**defaults)

    main_parser.add_argument(
        "-S",
        "--simulate-sequences-or-continous",
        choices=["sequences", "continuous"],
        action="store",
        help="Choose to simulate sequences or continous variable [default=%(default)s].",
    )
    main_parser.add_argument(
        "-G",
        "--geometry",
        choices=["complete graph", "chain graph", "random connectivity?"],
        action="store",
        help="Network/graph to simulate. [default=%(default)s].",
    )
    main_parser.add_argument(
        "-N",
        "--number-of-chromosomes-per-deme",
        action="store",
        type=int,
        help="Number of chromosomes per deme [default=%(default)s].",
    )
    main_parser.add_argument(
        "-P",
        "--number-of-chromosomes-per-invdividual",
        choices=range(1, 5),
        action="store",
        type=int,
        help="Number of ploidy by organism [default=%(default)s].",
    )
    main_parser.add_argument(
        "-k",
        "--number-of-demes",
        action="store",
        type=int,
        help="Number of demes in simulation[default=%(default)s].",
    )
    main_parser.add_argument(
        "-m",
        "--migration-rate",
        action="store",
        type=float,
        help="Migration rate across demes. Max = 1/k, where k is number of demes [default=%(default)s].",
    )
    main_parser.add_argument(
        "-mut",
        "--mutation-rate",
        action="store",
        type=float,
        help="Mutation rate for gene simulaiton [default=%(default)s].",
    )
    main_parser.add_argument(
        "-L",
        "--sequence-length",
        action="store",
        type=int,
        help="Sequence length for gene simulaiton [default=%(default)s].",
    )
    main_parser.add_argument(
        "-mean",
        "--mean",
        action="store",
        type=float,
        help="Starting mean for simulating continous varible[default=%(default)s].",
    )
    main_parser.add_argument(
        "-std",
        "--standard-deviation",
        action="store",
        type=float,
        help="Starting standard deviation for simulating continous varible[default=%(default)s].",
    )
    main_parser.add_argument(
        "-simT",
        "--simulation-time",
        action="store",
        type=int,
        help="Number of generations to run simulation for [default=%(default)s].",
    )
    main_parser.add_argument(
        "-simK",
        "--number-of-simulations",
        action="store",
        type=int,
        help="Number of trials of simulations to run [default=%(default)s].",
    )
    main_parser.add_argument(
        "-id",
        "--log-id",
        action="store",
        help="Filename for log files [default=%(default)s].",
    )
    main_parser.add_argument(
        "-o",
        "--output-prefix",
        action="store",
        help="Prefix for output files [default=%(default)s].",
    )
    main_parser.add_argument(
        "-s",
        "--seed",
        action="store",
        help="Seed for reproducibilty [default=%(default)s].",
    )

    main_args = main_parser.parse_args(args)

    # where did the value of each argument come from?
    logger.info("Seed: {}".format(main_args.seed))
    logger.info(
        "simulate_sequences_or_continous: {}".format(
            main_args.simulate_sequences_or_continous
        )
    )
    logger.info("geometry: {}".format(main_args.geometry))
    logger.info(
        "number_of_chromosomes_per_deme: {}".format(
            main_args.number_of_chromosomes_per_deme
        )
    )
    logger.info(
        "number_of_chromosomes_per_invdividual: {}".format(
            main_args.number_of_chromosomes_per_invdividual
        )
    )
    logger.info("number_of_demes: {}".format(main_args.number_of_demes))
    logger.info("migration_rate: {}".format(main_args.migration_rate))
    logger.info("mutation_rate: {}".format(main_args.mutation_rate))
    logger.info("sequence_length: {}".format(main_args.sequence_length))
    logger.info("simulation_time: {}".format(main_args.simulation_time))
    logger.info("number_of_simulations: {}".format(main_args.number_of_simulations))
    logger.info("log_id: {}".format(main_args.log_id))
    logger.info("output_prefix: {}".format(main_args.output_prefix))
    logger.info("mean: {}".format(main_args.mean))
    logger.info("standard_deviation: {}".format(main_args.standard_deviation))
    print("Hello, simulation begining")

    #    if args.config:
    #        configuration = read_config(args)

    #    #need override https://gist.github.com/gene1wood/9217725; args, unparsed = parser.parse_known_args()
    #    else:
    configuration = read_args(main_args)

    # inialize Simulator instance
    simulator = Simulator(configuration) #get rid of self.configuration...

    # to check
    setup_logger()
    logging.info("Beginning")  # for check

    if configuration.trait == "sequences":
        # SIMULATE
        simulator.simulate_multiple_trials_sequences()

        # analysis to OOP; be able to just pass simulator obj to Analysis(simulator)
        filenames_snapshot = (simulator.sequence_files, simulator.metadata_files)
        tag = simulator.resolution

        # ANALYSIS
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

    else:
        # SIMULATE
        simulator.simulate_multiple_trials_continuous()

        print(simulator.continuous_trial_files)
        simulator.plot_continuous()
        # pass to analysis OPP or something synomims to PDM


if __name__ == "__main__":
    main()
