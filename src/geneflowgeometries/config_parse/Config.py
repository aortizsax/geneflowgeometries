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

import numpy as np
import datetime

#combine config class and config parser

################################################################################
### Config
class Config:
    """
    Bin
    """

    ############################################################################
    ## Life cycle
    def __init__(self, config_dict):
        """
        """

        self.config_dict = config_dict

        # pass to local variables
        self.geometry = config_dict["simulator"]["geometry"]
        self.trait = config_dict["simulator"]["simulate_sequences_or_continous"]
        self.ghost_population_model = config_dict["simulator"]["ghost_population_model"]
        self.number_of_chromosomes = config_dict["simulator"][
            "number_of_chromosomes_per_deme"
        ]
        self.number_of_ploidy = config_dict["simulator"][
            "number_of_chromosomes_per_invdividual"
        ]
        self.number_of_demes = config_dict["simulator"]["number_of_demes"]
        self.migration_rate = config_dict["simulator"]["migration_rate"]
        self.migration_directionality_ratio = config_dict["simulator"]["migration_directionality_ratio"]
        self.restriction_rate = 1 - self.migration_rate
        self.number_generations = config_dict["simulator"]["simulation_time"]
        self.number_simulations = config_dict["simulator"]["number_of_simulations"]
        self.snapshot_times = [
       # 1,
       # int(0.2 * self.number_of_chromosomes),
       # int(0.4 * self.number_of_chromosomes),
        int(0.6 * self.number_of_chromosomes),
        int(0.8 * self.number_of_chromosomes),
        self.number_of_chromosomes,
        5 * self.number_of_chromosomes,
        8 * self.number_of_chromosomes,
        10 * self.number_of_chromosomes,
        12 * self.number_of_chromosomes,
        14 * self.number_of_chromosomes,
        16 * self.number_of_chromosomes,
        18 * self.number_of_chromosomes,
        20 * self.number_of_chromosomes,
    ]
        self.seed = config_dict["simulator"]["seed"]

        if self.trait == "continuous":
            self.start_mean = config_dict["simulator"]["mean"]
            self.start_std = config_dict["simulator"]["standard_deviation"]
        else:
            self.mutation_rate = config_dict["simulator"]["mutation_rate"]
            self.sequence_length = config_dict["simulator"]["sequence_length"]

        self.seed = config_dict["simulator"]["seed"]

        outfile_prefix = np.array(list(config_dict["simulator"].values()))
        outfile_prefix = "_".join(outfile_prefix)
        outfile_prefix = outfile_prefix.replace(" ", "-")
        self.outfile_prefix = outfile_prefix
        
        
        
        
        
        
        
        
        
        
        
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
    "Raised when the input value is greater than 1/k (complete graph) or 1/3 (chain graph), where k is the number of demes"
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
    #print(args)
    simulate_what = args.simulate_sequences_or_continous
    print("Simulating", simulate_what)

    geometry = args.geometry
    ghost_population_model = args.ghost_population_model
    number_of_chromosomes = args.number_of_chromosomes_per_deme
    number_of_ploidy = args.number_of_chromosomes_per_invdividual
    number_of_demes = args.number_of_demes
    migration_rate = args.migration_rate

    try:
        if geometry == "complete graph":
            if migration_rate > 1 / number_of_demes:
                raise InvalidMigrationRateException    
            else:
                print("Valid migration rate")
        elif geometry == "chain graph":
            if migration_rate > 1 / 3:
                raise InvalidMigrationRateException    
            else:
                print("Valid migration rate")
    except InvalidMigrationRateException:
        print("Exception occurred: Invalid migration rate")
        raise SystemExit(1)
        
    migration_directionality_ratio = args.migration_directionality_ratio
    mutation_rate = args.mutation_rate
    sequence_length = args.sequence_length
    simT = args.simulation_time
    start_mean = args.mean
    start_std = args.standard_deviation
#    snapshot_times = [
#        1,
#        int(0.2 * number_of_chromosomes),
#        int(0.4 * number_of_chromosomes),
#        int(0.6 * number_of_chromosomes),
#        int(0.8 * number_of_chromosomes),
#        number_of_chromosomes,
#        5 * number_of_chromosomes,
#        8 * number_of_chromosomes,
#        10 * number_of_chromosomes,
#        12 * number_of_chromosomes,
#        14 * number_of_chromosomes,
#        16 * number_of_chromosomes,
#        18 * number_of_chromosomes,
#        20 * number_of_chromosomes,
#    ]

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
            "geometry": geometry,
            "ghost_population_model":ghost_population_model,
            "number_of_chromosomes_per_deme": number_of_chromosomes,
            "number_of_chromosomes_per_invdividual": number_of_ploidy,
            "number_of_demes": number_of_demes,
            "migration_rate": migration_rate,
            "migration_directionality_ratio":migration_directionality_ratio,
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
