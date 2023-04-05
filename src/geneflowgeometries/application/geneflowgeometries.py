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


#directive from last week get ancsequences to OOP
#additionally continuous trait, argparse, 


#my to do 
#override args
#analysis OOP



import os
import pathlib  # barcharts next
import sys
import argparse
import datetime
import random
import numpy as np


from geneflowgeometries.simulate import (
    _AncestralDemeSequences,
    _ContinuousTraitEvolution,
)
from geneflowgeometries.simulate.Simulator import Simulator
from geneflowgeometries.calculate import analyze

import logging
from geneflowgeometries.log.multi_handlers_logger import setup_logger
from geneflowgeometries.config_parse.parse import read_config, read_args

from geneflowgeometries.config_parse.Config import Config
from configparser import ConfigParser



def main():
    args = sys.argv[1:]
    # setting the log level on the root logger must happen BEFORE any output
    logging_argparse = argparse.ArgumentParser(prog=__file__, add_help=False)
    logging_argparse.add_argument('-l', '--log-level', default='WARNING',
                                  help='set log level')
    logging_args, _ = logging_argparse.parse_known_args(args)

    try:
        logging.basicConfig(level=logging_args.log_level)
    except ValueError:
        logging.error("Invalid log level: {}".format(logging_args.log_level))
        sys.exit(1)

    logger = logging.getLogger(__name__)
    logger.info("Log level set: {}"
                .format(logging.getLevelName(logger.getEffectiveLevel())))

    # parse values from a configuration file if provided and use those as the
    # default values for the argparse arguments
    config_argparse = argparse.ArgumentParser(prog=__file__, add_help=False)
    config_argparse.add_argument('-c', '--config-file',
                                 help='path to configuration file')
    config_args, _ = config_argparse.parse_known_args(args)

    defaults = {
        "simulate_sequences_or_continous":"sequences",
        "geometry":"complete graph",
        "number_of_chromosomes_per_deme":100,
        "number_of_chromosomes_per_invdividual":2,
        "number_of_demes":4,
        "migration_rate":0.1,
        "mutation_rate":0.0004,
        "sequence_length":1000,
        "simulation_time":2000,
        "number_of_simulations":200,
        "log_id":"log_DATE",
        "output_prefix":"output",
        "seed":12345,
        "mean":0,
        "standard_deviation":1,
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

        defaults.update(dict(config_parser.items('SIMULATOR')))

    # parse the program's main arguments using the dictionary of defaults and
    # the previous parsers as "parent' parsers
    parsers = [logging_argparse, config_argparse]
    main_parser = argparse.ArgumentParser(prog=__file__, parents=parsers)
    main_parser.set_defaults(**defaults)
    main_parser.add_argument('-1', '--option1')
    main_parser.add_argument('-2', '--option2')


 
    
    main_parser.add_argument(
        "-S",
        "--simulate-sequences-or-continous",
        choices=["sequences","continuous"],
        action="store",
        default="sequences",
        help="Choose to simulate sequences or continous variable [default=%(default)s].",
    )
    main_parser.add_argument(
        "-G",
        "--geometry",
        choices=["complete graph","chain graph","random connectivity?"],
        action="store",
        default="complete graph",
        help="Network/graph to simulate. [default=%(default)s].",
    )
    main_parser.add_argument(
        "-N",
        "--number-of-chromosomes-per-deme",
        action="store",
        default=100,
        type=int,
        help="Number of chromosomes per deme [default=%(default)s].",
    )
    main_parser.add_argument(
        "-P",
        "--number-of-chromosomes-per-invdividual",
        choices=range(1,5),
        action="store",
        default=2,
        type=int,
        help="Number of ploidy by organism [default=%(default)s].",
    )
    main_parser.add_argument( 
        "-k",
        "--number-of-demes",
        action="store",
        default=4,
        type=int,
        help="Number of demes in simulation[default=%(default)s].",
    )
    main_parser.add_argument(
        "-m",
        "--migration-rate",
        action="store",
        default=0,
        type=float,
        help="Migration rate across demes. Max = 1/k, where k is number of demes [default=%(default)s].",
    )
    main_parser.add_argument(
        "-mut",
        "--mutation-rate",
        action="store",
        default=0.0004,
        type=float,
        help="Mutation rate for gene simulaiton [default=%(default)s].",
    )
    main_parser.add_argument(
        "-L",
        "--sequence-length",
        action="store",
        default=1000,
        type=int,
        help="Sequence length for gene simulaiton [default=%(default)s].",
    )
    main_parser.add_argument(
        "-mean",
        "--mean",
        action="store",
        default=0,
        type=float,
        help="Starting mean for simulating continous varible[default=%(default)s].",
    )
    main_parser.add_argument(
        "-std",
        "--standard-deviation",
        action="store",
        default=1,
        type=float,
        help="Starting standard deviation for simulating continous varible[default=%(default)s].",
    )
    main_parser.add_argument(
        "-simT",
        "--simulation-time",
        action="store",
        default=2000,
        type=int,
        help="Number of generations to run simulation for [default=%(default)s].",
    )
    main_parser.add_argument(
        "-simK",
        "--number-of-simulations",
        action="store",
        default=200,
        type=int,
        help="Number of trials of simulations to run [default=%(default)s].",
    )
    main_parser.add_argument(
        "-id",
        "--log-id",
        action="store",
        default="log-DATE",
        help="Filename for log files [default=%(default)s].",
    )
    main_parser.add_argument(
        "-o",
        "--output-prefix",
        action="store",
        default="output",
        help="Prefix for output files [default=%(default)s].",
    )
    main_parser.add_argument(
        "-s",
        "--seed",
        action="store",
        default="random",
        help="Seed for reproducibilty [default=%(default)s].",
    )
    
    

    main_args = main_parser.parse_args(args)
    
    # where did the value of each argument come from?
    logger.info("Option 1: {}".format(main_args.option1))
    logger.info("Option 2: {}".format(main_args.option2))
    logger.info("Option 2: {}".format(main_args.seed))

    
    print("Hello, simulation begining")





#    if args.config:
#        Configuration = read_config(args)

#    #need override https://gist.github.com/gene1wood/9217725; args, unparsed = parser.parse_known_args()
#    else:
    Configuration = read_args(main_args)


    # inialize Simulator instance
    simulator = Simulator(Configuration)

    # to check
    setup_logger()
    logging.info("Beginning")  # for check

    if Configuration.trait == "sequences":
        # SIMULATE
        simulator.ancestralDemeSequences()
        
        
        
        #analysis to OOP; be able to just pass simulator obj to Analysis(simulator)
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
        simulator.continuousTraitEvolution()
        
        print(simulator.continuous_data_filename)
        #pass to analysis OPP or something synomims to PDM 


if __name__ == "__main__":
    main()
