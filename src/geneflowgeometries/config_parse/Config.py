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
        self.number_of_chromosomes = config_dict["simulator"][
            "number_of_chromosomes_per_deme"
        ]
        self.number_of_ploidy = config_dict["simulator"][
            "number_of_chromosomes_per_invdividual"
        ]
        self.number_of_demes = config_dict["simulator"]["number_of_demes"]
        self.migration_rate = config_dict["simulator"]["migration_rate"]
        self.restriction_rate = 1 - self.migration_rate
        self.number_generations = config_dict["simulator"]["simulation_time"]
        self.snapshot_times = [
            1,
            self.number_of_chromosomes,
            5 * self.number_of_chromosomes,
            10 * self.number_of_chromosomes,
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
