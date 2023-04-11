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

import json
import random
import datetime

from geneflowgeometries.config_parse.Config import Config


class InvalidMigrationRateException(Exception):
    "Raised when the input value is greater than 1/k, where k is the number of demes"
    pass


def read_config(args):
    # Load the configuration file
    with open(args.config) as f:
        config = f.read()

    date = datetime.datetime.now()  # pull up in scripts
    date = str(date).split(" ")[0]
    config_dict = {"simulator": {}}
    config_dict["simulator"]["date"] = date

    json_dict = json.loads(config)

    config_dict["simulator"].update(json_dict["simulator"])

    simulate_what = config_dict["simulator"]["simulate_sequences_or_continous"]

    configuration = Config(config_dict)

    return configuration


def read_args(args):
    # **********************************args to module
    # Pass args
    print(args)
    simulate_what = args.simulate_sequences_or_continous
    print("Simulating", simulate_what)

    Geometry = args.geometry
    number_of_chromosomes = args.number_of_chromosomes_per_deme
    number_of_ploidy = args.number_of_chromosomes_per_invdividual
    number_of_demes = args.number_of_demes
    migration_rate = args.migration_rate

    try:
        if migration_rate > 1 / number_of_demes:
            raise InvalidMigrationRateException
        else:
            print("Valid migration rate")
    except InvalidMigrationRateException:
        print("Exception occurred: Invalid migration rate")
        raise SystemExit(1)

    mutation_rate = args.mutation_rate
    sequence_length = args.sequence_length
    simT = args.simulation_time
    start_mean = args.mean
    start_std = args.standard_deviation
    snapshot_times = [
        1,
        number_of_chromosomes,
        5 * number_of_chromosomes,
        10 * number_of_chromosomes,
        20 * number_of_chromosomes,
    ]

    # random number
    if args.seed == "random":
        seed_range = 2 ** 32 - 1
        seed = int(random.uniform(0, seed_range))
    else:
        seed = int(args.seed)

    date = datetime.datetime.now()  # pull up in scripts
    date = str(date).split(" ")[0]

    config_dict = {
        "simulator": {
            "date": date,
            "simulate_sequences_or_continous": simulate_what,
            "geometry": Geometry,
            "number_of_chromosomes_per_deme": number_of_chromosomes,
            "number_of_chromosomes_per_invdividual": number_of_ploidy,
            "number_of_demes": number_of_demes,
            "migration_rate": migration_rate,
            "mutation_rate": mutation_rate,
            "sequence_length": sequence_length,
            "mean": start_mean,
            "standard_deviation": start_std,
            "simulation_time": simT,
            "number_of_simulations": args.number_of_simulations,
            "log_id": "log_DATE",
            "output_prefix": "output",
            "seed": seed,
        }
    }

    configuration = Config(config_dict)

    return configuration
