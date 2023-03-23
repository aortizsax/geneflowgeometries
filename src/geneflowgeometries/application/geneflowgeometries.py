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

# make config class
# get set restriction

#make this use OOP classes 


import os
import pathlib  # barcharts next
import sys
import argparse
import datetime
import random
import numpy as np


from geneflowgeometries.simulate import (
    _AncestralDemeSequences,
    _ContinuousTraitEvolution
    )
from geneflowgeometries.simulate.Simulator import Simulator
from geneflowgeometries.calculate import analyze

import logging
from geneflowgeometries.log.multi_handlers_logger import setup_logger
from geneflowgeometries.config_parse.parse import read_config, read_args

from geneflowgeometries.config_parse.Config import Config



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

    args = parser.parse_args()
    print("Hello, simulation begining")

    if args.config != "FILE":
        (config_dict, simulate_what) = read_config(args)

    else:
        (snapshot_times, config_dict, simulate_what) = read_args(args)


    #START get set outfile_prefix
    outfile_prefix = np.array(list(config_dict["simulator"].values()))
    outfile_prefix = "_".join(outfile_prefix)
    outfile_prefix = outfile_prefix.replace(" ", "-")
    #END get set outfile_prefix
    
    config_dict['outfile_prefix'] = outfile_prefix
    
    simulator = Simulator(config_dict)
    
    setup_logger()
    logging.info("Beginning") #for check

    if simulate_what == "sequences":
        
        #SIMULATE
        (filenames_snapshot, tag) = simulator.ancestralDemeSequences()


        #ANALYSIS
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
        #SIMULATE
        csvfile = simulator.continuousTraitEvolution()


if __name__ == "__main__":
    main()
