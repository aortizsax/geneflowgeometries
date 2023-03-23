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
#
# Simulation
#
###########################################################################
## Reader Implementation

# Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import bernoulli, norm, poisson
from scipy import stats
from scipy.stats import kurtosis, skew
import seaborn as sns
import math
import pandas as pd
import string
import datetime

import random

import logging
from geneflowgeometries.run_time.write_and_log import write_log_continuous


# Classes
class Continuous_trait_deme:
    def __init__(self, mean, std, deme):
        self.mean = mean
        self.std = std
        self.deme = deme

    def __str__(self):
        return "Deme:" + self.deme + ",Mean:" + self.mean + ",Std:" + self.std


alphabet = string.ascii_lowercase


def continuous_trait_evoluion(self):
    # pass to local variables
    config_dict = self.config_dict
    outfile_prefix = config_dict['outfile_prefix']
    Geometry = config_dict["simulator"]["geometry"]
    number_of_chromosomes = config_dict["simulator"]["number_of_chromosomes_per_deme"]
    number_of_demes = config_dict["simulator"]["number_of_demes"]
    migration_rate = config_dict["simulator"]["migration_rate"]
    number_generations = config_dict["simulator"]["simulation_time"]
    start_mean = config_dict["simulator"]["mean"]
    start_std = config_dict["simulator"]["standard_deviation"]
    seed = config_dict["simulator"]["seed"]

    FASTAFILENAMES = []  # csv means, std
    CSVFILENAMES = []  # pop data
    demes = [[]]
    labels = []

    logging.info("Simulating")

    # set seed
    random.seed(seed)
    np.random.seed(seed)

    for i in range(number_of_demes):
        labels.append(alphabet[i])
        demes[0].append(Continuous_trait_deme(start_mean, start_std, alphabet[i]))

    csv_header = ",generation,"
    csv_header += ",".join(labels)

    csv_list = [csv_header]  # add labels

    # EXPERIMENT
    for generation in range(number_generations):
        temp_demes_mean = [0] * number_of_demes
        temp_demes_std = [0] * number_of_demes
        # simulate draws foreward in time

        temp_draws = [0.0] * number_of_demes

        demes_mean = []
        demes_std = []

        for i, deme_traits in enumerate(demes[generation]):
            demes_mean.append(deme_traits.mean)
            demes_std.append(deme_traits.std)

        for i, deme_traits in enumerate(demes[generation]):
            # migration arrary
            migration_array = [migration_rate] * number_of_demes
            migration_array[i] = 1 - (migration_rate * (number_of_demes - 1))

            number_draw = int(migration_array[i] * number_of_chromosomes)
            temp_draws[i] = np.random.normal(
                loc=deme_traits.mean, scale=deme_traits.std, size=number_draw
            )

            migration_array = migration_array[:i] + migration_array[i + 1 :]

            temp_std = demes_std[:i] + demes_std[i + 1 :]
            for j, meanj in enumerate(demes_mean[:i] + demes_mean[i + 1 :]):
                number_draw = int(migration_array[j] * number_of_chromosomes)
                np.append(
                    temp_draws[i],
                    np.random.normal(loc=meanj, scale=temp_std[j], size=number_draw),
                )

            temp_demes_mean[i] = round(np.mean(temp_draws[i]), 2)
            temp_demes_std[i] = round(np.std(temp_draws[i]), 2)

        demes_mean = temp_demes_mean
        demes_std = temp_demes_std

        demes.append([])
        for i in range(number_of_demes):
            mean = demes_mean[i]
            std = demes_std[i]
            label = alphabet[i]
            demes[generation + 1].append(Continuous_trait_deme(mean, std, label))

        csv_list.append(str(generation) + "," + "".join(str(demes_mean))[1:-1])
    print(csv_list)

    CSVFILENAME = write_log_continuous(csv_list, outfile_prefix)

    return CSVFILENAME
