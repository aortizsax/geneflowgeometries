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
import json

#from geneflowgeometries.run_time.write_and_log import write_log_sequences
#from geneflowgeometries.run_time.write_and_log import write_log_continuous

from geneflowgeometries.config_parse.Config import Config


# add analysis? or is this like pdm||new object entirely


# Classes
################################################################################
### Chromosome
class Chromosome:
    def __init__(self, sequence, ancenstral_deme):
        self.sequence = sequence
        self.ancenstral_deme = ancenstral_deme

    def __str__(self):
        return self.sequence + self.ancenstral_deme


class Continuous_trait_deme:
    def __init__(self, mean, std, deme):
        self.mean = mean
        self.std = std
        self.deme = deme

    def __str__(self):
        return "Deme:" + self.deme + ",Mean:" + self.mean + ",Std:" + self.std

################################################################################
### Simulator
class Simulator:
    nuc_alphabet = ["A", "T", "C", "G"]

    alphabet = string.ascii_lowercase

    ############################################################################
    ## Life-cycle

    # configuation is a Config obj
    def __init__(self, Configuration):
        self.Configuration = Configuration
        # intialize demes
        self.demes = [[]] * self.Configuration.number_of_demes

        # set seed
        random.seed(self.Configuration.seed)
        np.random.seed(self.Configuration.seed)

        self.labels()
        self.migrationMatrix()

        return

    ############################################################################
    ## Initialze Demes
    def labels(self):
        self.labels = []
        for i, deme_sequences in enumerate(self.demes):
            self.labels.append(self.alphabet[i])
        return self.labels


    ############################################################################
    ## Initialize Sequnce Model
    
    def startingSequences(self):
        # set orginating deme sequences
        self.start_seqs = []
        for i, deme_sequences in enumerate(self.demes):
            seq = ""
            for j in range(self.Configuration.sequence_length):
                seq += random.choice(self.nuc_alphabet)

            self.start_seqs.append(seq)
        return
        
    def startingDemeSequences(self):
        # intalize class
        for i, deme_sequences in enumerate(self.demes):
            self.demes[i] = [0] * self.Configuration.number_of_chromosomes

            for k in range(self.Configuration.number_of_chromosomes):
                deme_chromosome = Chromosome(self.start_seqs[i], self.labels[i])
                self.demes[i][k] = deme_chromosome

        return

    def migrationMatrix(self):
        # migration matrix
        self.migration_matrix = [0] * self.Configuration.number_of_demes
        for i in range(self.Configuration.number_of_demes):
            self.migration_matrix[i] = [
                self.Configuration.migration_rate
            ] * self.Configuration.number_of_demes

        for i in range(self.Configuration.number_of_demes):
            print(self.migration_matrix)
            self.migration_matrix[i][i] = 1 - (
                self.Configuration.migration_rate
                * (self.Configuration.number_of_demes - 1)
            )

        return self.migration_matrix

    ############################################################################
    ## mutate
    def mutate(self):
        # mutate function
        for i, deme_sequences in enumerate(self.demes):
            for j, chromosome in enumerate(self.deme_sequences):
                which_sites_mut = np.random.uniform(0, 1, self.sequence_length)
                which_sites_mut = which_sites_mut < self.mutation_rate
                which_sites_mut = np.where(which_sites_mut)[0]
                if any(which_sites_mut):
                    for l, mutate_allele in enumerate(
                        which_sites_mut
                    ):  # add find = j; [:j]+[j:]
                        allele = chromosome.sequence[mutate_allele]
                        # got rid of deep copy
                        isAllele = self.nuc_alphabet.find(allele)
                        possible_mutations = self.nuc_alphabet[:isAllele]
                        possible_mutations = +self.nuc_alphabet[j + 1 :]
                        demes[i][j].sequence = (
                            demes[i][j].sequence[:mutate_allele]
                            + random.choice(possible_mutations)
                            + demes[i][j].sequence[mutate_allele + 1 :]
                        )

    ############################################################################
    ## Simulate
    def ancestralDemeSequences(self):
        self.sequence_files = []
        self.metadata_files = []
        self.startingSequences()
        self.startingDemeSequences()
        # EXPERIMENT
        for generation in range(self.Configuration.number_generations):
            self.generation = generation
            temp_demes = []

            for i, temp_sequences in enumerate(self.demes):
                temp_draws = []
                # migration arrary
                migration_array = self.migration_matrix[i]
                for k in range(self.Configuration.number_of_chromosomes):
                    sample_demes = range(self.Configuration.number_of_demes)

                    from_deme_draw = np.random.choice(sample_demes, p=migration_array)

                    from_sequence_draw = random.randint(
                        0, self.Configuration.number_of_chromosomes - 1
                    )
                    temp_draws.append(self.demes[from_deme_draw][from_sequence_draw])
                random.shuffle(temp_draws)
                temp_demes.append(temp_draws)

            # update demes/pops with temp
            self.demes = temp_demes

            self.writeSnapShots()

        return
        
    ############################################################################
    ## Initialize Continuous Model
    def startingContinuousTrait(self):
        #####
        for i in range(self.Configuration.number_of_demes):
            self.labels.append(self.alphabet[i])
            self.demes[0].append(
                Continuous_trait_deme(
                    self.Configuration.start_mean,
                    self.Configuration.start_std,
                    self.alphabet[i],
                )
            )

        return

    ############################################################################
    ## Simulate
    def continuousTraitEvolution(self):
        self.demes = [[]]  # * self.Configuration.number_of_demes
        self.startingContinuousTrait()

        self.csv_header = ",generation,"
        self.csv_header += ",".join(self.labels)


        self.continuous_data_list = [self.csv_header]  

        # EXPERIMENT
        for generation in range(self.Configuration.number_generations):
            self.generation = generation
            print("Generation:", generation)
            temp_demes_mean = [0] * self.Configuration.number_of_demes
            temp_demes_std = [0] * self.Configuration.number_of_demes

            temp_draws = [0.0] * self.Configuration.number_of_demes

            demes_mean = []
            demes_std = []

            for i, deme_traits in enumerate(self.demes[generation]):
                demes_mean.append(deme_traits.mean)
                demes_std.append(deme_traits.std)

            for i, deme_traits in enumerate(self.demes[generation]):
                # migration arrary
                migration_array = self.migration_matrix[i]

                number_draw = int(
                    migration_array[i] * self.Configuration.number_of_chromosomes
                )
                temp_draws[i] = np.random.normal(
                    loc=deme_traits.mean, scale=deme_traits.std, size=number_draw
                )

                migration_array = migration_array[:i] + migration_array[i + 1 :]

                temp_std = demes_std[:i] + demes_std[i + 1 :]
                for j, meanj in enumerate(demes_mean[:i] + demes_mean[i + 1 :]):
                    number_draw = int(
                        migration_array[j] * self.Configuration.number_of_chromosomes
                    )
                    
                    np.append(
                        temp_draws[i],
                        np.random.normal(
                            loc=meanj, scale=temp_std[j], size=number_draw
                        ),
                    )

                temp_demes_mean[i] = round(np.mean(temp_draws[i]), 2)
                temp_demes_std[i] = round(np.std(temp_draws[i]), 2)

            demes_mean = temp_demes_mean
            demes_std = temp_demes_std

            self.demes.append([])
            for i in range(self.Configuration.number_of_demes):
                mean = demes_mean[i]
                std = demes_std[i]
                label = self.alphabet[i]

                self.demes[generation + 1].append(
                    Continuous_trait_deme(mean, std, label)
                )

            self.continuous_data_list.append(str(generation) + "," + "".join(str(demes_mean))[1:-1])

        self.write_continuous_data()

        return 
        
    ############################################################################
    ## Write to files 
    def writeSnapShots(self):
        # change to runtime*****************************    
        
        #print("Generation:", self.generation)
        if (self.generation in self.Configuration.snapshot_times) | (
            self.generation == self.Configuration.number_generations - 1
        ):
            (sequence_file, metadata_file, self.Resolution) = self.write_sequence_data()

            self.sequence_files.append(sequence_file)
            self.metadata_files.append(metadata_file)

        return

    def log_analysis():
        logging.info("Analyzing")
        return


    def write_continuous_data(self):
        # write deme data for analysis downstream
        # open file in write mode
        
        self.continuous_data_filename = self.Configuration.outfile_prefix + ".csv"
        with open(self.continuous_data_filename, "w") as fp:
            for item in self.continuous_data_list:
                # write each item on a new line
                fp.write("%s\n" % item)
            print("Done", self.continuous_data_filename)

        return 


    def write_sequence_data(self):
        print("snapshot", self.generation)
        #self.Configuration
        taxon_labels = []
        csv_list = [",IDV,SUBPOP,ANC"]
        fasta_list = []
        for deme_count, deme in enumerate(self.demes):
            for chromose_count, chromosome in enumerate(deme):
                sequence_bases = chromosome.sequence

                fasta_list.append(
                    ">" + self.labels[deme_count] + ".chromsome." + str(chromose_count)
                )

                chromosome_count = self.labels[deme_count] + ".chromsome." + str(chromose_count)

                ind_count_label = (
                    self.labels[deme_count] + "|" + str(math.floor(chromose_count / 2))
                )
                # change to number ploidy

                csv_to_append = ",".join(
                    [
                        chromosome_count,
                        ind_count_label,
                        self.labels[deme_count],
                        chromosome.ancenstral_deme,
                    ]
                )
                csv_list.append(csv_to_append)

                fasta_list.append(sequence_bases)
            print()

       
        logging.info(json.dumps(self.Configuration.outfile_prefix))  # as same line term = ,
        logging.info(json.dumps(self.generation)) 
        
        
        # write deme data for analysis downstream
        # open file in write mode
        self.metadata_filename = self.Configuration.outfile_prefix + "_" + str(self.generation) + ".csv"

        with open(self.metadata_filename, "w") as fp:
            for item in csv_list:
                # write each item on a new line
                fp.write("%s\n" % item)
            print("Done", self.generation)
            
            
            

        # Write deme sequences for analyis downstream
        # open file in write mode
        self.sequence_filename = self.Configuration.outfile_prefix + "_" + str(self.generation) + ".fasta"
        with open(self.sequence_filename, "w") as fp:
            for item in fasta_list:
                # write each item on a new line
                fp.write("%s\n" % item)
            print("Done", self.generation)

        self.resolution = csv_list[0].split(",")[2]

        return (self.sequence_filename, self.metadata_filename, self.resolution)
        
        
        
        
        
        
        
        
        
