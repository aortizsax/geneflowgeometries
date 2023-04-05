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
import copy
import logging


from geneflowgeometries.run_time.write_and_log import write_log_sequences


# Classes
class Chromosome:
    def __init__(self, sequence, ancenstral_deme):
        self.sequence = sequence
        self.ancenstral_deme = ancenstral_deme

    def __str__(self):
        return self.sequence + self.ancenstral_deme


nuc_alphabet = ["A", "T", "C", "G"]

alphabet = string.ascii_lowercase


# self is a Simulator


def ancestral_deme_sequences(self):
    # pass to local variables
    config_dict = self.config_dict
    outfile_prefix = config_dict["outfile_prefix"]
    Geometry = config_dict["simulator"]["geometry"]
    number_of_chromosomes = config_dict["simulator"]["number_of_chromosomes_per_deme"]
    number_of_ploidy = config_dict["simulator"]["number_of_chromosomes_per_invdividual"]
    number_of_demes = config_dict["simulator"]["number_of_demes"]
    migration_rate = config_dict["simulator"]["migration_rate"]
    mutation_rate = config_dict["simulator"]["mutation_rate"]
    sequence_length = config_dict["simulator"]["sequence_length"]
    number_generations = config_dict["simulator"]["simulation_time"]
    snapshot_times = [
        1,
        number_of_chromosomes,
        5 * number_of_chromosomes,
        10 * number_of_chromosomes,
        20 * number_of_chromosomes,
    ]
    seed = config_dict["simulator"]["seed"]

    logger = logging.getLogger("geneflowgeometries.main.log")

    logging.info("Simulating")

    FASTAFILENAMES = []
    CSVFILENAMES = []
    demes = [[]] * number_of_demes
    labels = []

    # set seed
    random.seed(seed)
    np.random.seed(seed)

    for i, deme_sequences in enumerate(demes):
        labels.append(alphabet[i])

    # set orginating deme sequences
    start_seqs = []
    for i, deme_sequences in enumerate(demes):
        seq = ""
        for j in range(sequence_length):
            seq += random.choice(nuc_alphabet)

        start_seqs.append(seq)

    # intalize class
    for i, deme_sequences in enumerate(demes):
        demes[i] = [0] * number_of_chromosomes

        for k in range(number_of_chromosomes):
            deme_chromosome = Chromosome(start_seqs[i], labels[i])
            demes[i][k] = deme_chromosome

    # migration matrix
    migration_matrix = [0] * number_of_demes
    for i in range(number_of_demes):
        migration_matrix[i] = [migration_rate] * number_of_demes

    for i in range(number_of_demes):
        print(migration_matrix)
        migration_matrix[i][i] = 1 - (migration_rate * (number_of_demes - 1))

    # EXPERIMENT
    for generation in range(number_generations):
        temp_demes = []

        for i, temp_sequences in enumerate(demes):
            temp_draws = []
            # migration arrary
            migration_array = migration_matrix[i]
            for k in range(number_of_chromosomes):
                sample_demes = range(number_of_demes)

                from_deme_draw = np.random.choice(sample_demes, p=migration_array)

                from_sequence_draw = random.randint(0, number_of_chromosomes - 1)
                temp_draws.append(demes[from_deme_draw][from_sequence_draw])
            random.shuffle(temp_draws)
            temp_demes.append(temp_draws)

        # update demes/pops with temp
        demes = temp_demes

        # mutate function
        for i, deme_sequences in enumerate(demes):
            for j, chromosome in enumerate(deme_sequences):
                which_sites_mut = np.random.uniform(0, 1, sequence_length)
                which_sites_mut = which_sites_mut < mutation_rate
                which_sites_mut = np.where(which_sites_mut)[0]
                if any(which_sites_mut):
                    for l, mutate_allele in enumerate(
                        which_sites_mut
                    ):  # add find = j; [:j]+[j:]
                        allele = chromosome.sequence[mutate_allele]
                        possible_mutations = copy.deepcopy(nuc_alphabet)
                        possible_mutations.remove(allele)
                        demes[i][j].sequence = (
                            demes[i][j].sequence[:mutate_allele]
                            + random.choice(possible_mutations)
                            + demes[i][j].sequence[mutate_allele + 1 :]
                        )

        # change to runtime*****************************

        if (generation in snapshot_times) | (generation == number_generations - 1):
            (FASTAFILENAME, CSVFILENAME, Resolution) = write_log_sequences(
                demes, outfile_prefix, generation, labels, config_dict
            )
            FASTAFILENAMES.append(FASTAFILENAME)
            CSVFILENAMES.append(CSVFILENAME)

    return ((FASTAFILENAMES, CSVFILENAMES), Resolution)
