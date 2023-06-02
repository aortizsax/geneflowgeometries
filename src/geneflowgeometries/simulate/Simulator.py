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

"""
This module handles the core definition of the simulation data class as well as 
the analysis data class that performs pop gen and entropy calculations on the 
data.
"""

""" notes
Data class (config)
    no more self.configuration.....
method class (simulator)
    act on next data class
    More modulization split up mig_matrix set to case dependent mig_matrix setup
data class (sequqnces/deme probs)
    make this 
method class (analysis)
    move from neymetyaay?

"""

"""
Make test for objects and arg parser

"""


import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
import seaborn as sns

import math
import pandas as pd
import logging
import json
from geneflowgeometries.config_parse.Config import Config

# add analysis? or is this like pdm||new object entirely

# Classes
################################################################################
## Chromosome
class Chromosome:
    """
    """

    ############################################################################
    # Life-cycle
    def __init__(self, sequence, ancenstral_deme):
        """
        Construct a new ...
        
        Keyword arguments  
        
        Notes 
        -----
        
        Parameters
        ----------
        
        Returns
        -------
        """
        self.sequence = sequence
        self.ancenstral_deme = ancenstral_deme

    ############################################################################
    ### Representation

    def __str__(self):
        """
                Composes and returns and representation of the 

        Parameters
        ----------
        Returns
        -------
        str
            ---

        Example
        -------


        """
        return self.sequence + self.ancenstral_deme


class Continuous_trait_deme:
    """
    """

    ############################################################################
    ### Life-cycle
    def __init__(self, mean, std, deme):
        """
        Construct a new ...
        
        Keyword arguments  
        
        Notes 
        -----
        
        Parameters
        ----------
        
        Returns
        -------        
        """
        self.mean = mean
        self.std = std
        self.deme = deme

    ############################################################################
    ### Representation

    def __str__(self):
        """
        Composes and returns and representation of the 

        Parameters
        ----------
        Returns
        -------
        str
            ---

        Example
        -------


        """
        return "Deme:" + self.deme + ",Mean:" + self.mean + ",Std:" + self.std


################################################################################
## Simulator
class Simulator:
    """
    """

    import string

    alphabet = string.ascii_lowercase

    ############################################################################
    ### Life-cycle

#    def __init__(self, configuration):  # configuation is a Config obj
    def __init__(self, configuration):  # configuation is a Config obj
        """
        Construct a new ...
        
        Keyword arguments  
        
        Notes 
        -----
        
        Parameters
        ----------
        
        Returns
        -------
        """
        self.configuration = configuration #loop thru parameters or ????
        # intialize demes
        self.demes = [[]] * self.configuration.number_of_demes

        # intialize rand generator and set seed
        self.rng = default_rng(self.configuration.seed)

        self.set_labels()
        self.calculate_migration_matrix()


    def calculate_migration_matrix(self):
        """
        
        #"No ghost","Equi-bidirectional migration", "Equi-directional migration", "High bidirectional migration","High directional migration"
        
        """
        if self.configuration.ghost_population_model == "Equi-bidirectional migration": #complete
            self.migration_matrix = [0] * self.configuration.number_of_demes
            for i in range(self.configuration.number_of_demes): #setupper and lower triangle
                self.migration_matrix[i] = [
                    self.configuration.migration_rate
                ] * self.configuration.number_of_demes

            for i in range(self.configuration.number_of_demes): #set diagonals
                self.migration_matrix[i][i] = 1 - (
                    self.configuration.migration_rate
                    * (self.configuration.number_of_demes - 1)
                )
            print(self.migration_matrix)
            return None
            
        elif self.configuration.ghost_population_model == "Equi-directional migration": #complete
            self.migration_matrix = [0] * self.configuration.number_of_demes
            number_sampled = int(self.configuration.number_of_demes * 2 / 3)
            number_self = int(self.configuration.number_of_demes / 3)
            number_ghost = int(self.configuration.number_of_demes / 3)
            sampled_range = range(number_sampled)
            ghost_range = range(number_sampled,number_ghost)
            
            
               # for i in sampled_range:
                
            
            for i in range(self.configuration.number_of_demes): #setupper and lower triangle
                self.migration_matrix[i] = [
                    self.configuration.migration_rate 
                ] * self.configuration.number_of_demes

            for i in range(self.configuration.number_of_demes): #set diagonals
                self.migration_matrix[i][i] = 1 - (
                    self.configuration.migration_rate
                    * (self.configuration.number_of_demes - 1)
                )
            print(self.migration_matrix)
            return None        
         
        # migration matrix complete graph
        if self.configuration.geometry == 'complete graph':
            self.migration_matrix = [0] * self.configuration.number_of_demes
            for i in range(self.configuration.number_of_demes): #setupper and lower triangle
                self.migration_matrix[i] = [
                    self.configuration.migration_rate
                ] * self.configuration.number_of_demes

            for i in range(self.configuration.number_of_demes): #set diagonals
                self.migration_matrix[i][i] = 1 - (
                    self.configuration.migration_rate
                    * (self.configuration.number_of_demes - 1)
                )
        elif self.configuration.geometry == 'chain graph':
            # migration matrix chain graph
            if self.configuration.migration_directionality_ratio == 1:
                self.migration_matrix = [0] * self.configuration.number_of_demes
                for i in range(self.configuration.number_of_demes):
                    self.migration_matrix[i] = [0] * self.configuration.number_of_demes

                for i in range(self.configuration.number_of_demes):
                    if i == 0:
                        self.migration_matrix[i][0] = 1 - self.configuration.migration_rate
                        self.migration_matrix[i][1] = self.configuration.migration_rate
                        
                    elif i == (len(self.labels)-1):
                        self.migration_matrix[i][-2] = self.configuration.migration_rate
                        self.migration_matrix[i][-1] = 1 -self.configuration.migration_rate
                    else:
                        self.migration_matrix[i][i-1] = self.configuration.migration_rate
                        self.migration_matrix[i][i] = 1 - 2 * (self.configuration.migration_rate)
                        self.migration_matrix[i][i+1] = self.configuration.migration_rate
            else:
                out_mig_ratio = self.configuration.migration_directionality_ratio
                in_mig_ratio = 1 / self.configuration.migration_directionality_ratio 
                # migration matrix complete graph
                # migration matrix chain graph
                self.migration_matrix = [0] * self.configuration.number_of_demes
                for i in range(self.configuration.number_of_demes):
                    self.migration_matrix[i] = [0] * self.configuration.number_of_demes

                for i in range(self.configuration.number_of_demes):
                    if i == 0:
                        self.migration_matrix[i][0] = 1 - self.configuration.migration_rate * out_mig_ratio
                        self.migration_matrix[i][1] = self.configuration.migration_rate * out_mig_ratio
                        
                    elif i == (len(self.labels)-1):
                        self.migration_matrix[i][-2] = self.configuration.migration_rate 
                        self.migration_matrix[i][-1] = 1 -self.configuration.migration_rate 
                    else:
                        self.migration_matrix[i][i-1] = self.configuration.migration_rate 
                        self.migration_matrix[i][i] = 1 - (1+out_mig_ratio) * (self.configuration.migration_rate)
                        self.migration_matrix[i][i+1] = self.configuration.migration_rate * out_mig_ratio
                
        print(self.migration_matrix)
        return None

    ############################################################################
    ### Initialze Demes

    def set_labels(self):
        """
        """
        self.labels = []
        for i, deme_sequences in enumerate(self.demes):
            self.labels.append(self.alphabet[i])
        return self.labels

    ############################################################################ 
    ### Initialize Sequnce Model

    def set_starting_sequences(self):
        """
        """

        # set orginating deme sequences
        self.start_seqs = []
        for i, deme_sequences in enumerate(self.demes):
            seq = "".join(
                self.rng.choice(self.nuc_alphabet, self.configuration.sequence_length) 
            )
            self.start_seqs.append(seq)
        return None

    def set_starting_demes_sequences(self):

        """
        """
        # intalize class
        for i, deme_sequences in enumerate(self.demes):
            self.demes[i] = [0] * self.configuration.number_of_chromosomes

            for k in range(self.configuration.number_of_chromosomes):
                deme_chromosome = Chromosome(self.start_seqs[i], self.labels[i]) 
                self.demes[i][k] = deme_chromosome

        return None

    ############################################################################ 
    ### Mutate sequence

    def mutate_sequences(self):
        """
        """
        for i, deme_sequences in enumerate(self.demes):
            for j, chromosome in enumerate(deme_sequences):
                which_sites_mut = self.rng.random(size = self.configuration.sequence_length)

                which_sites_mut = which_sites_mut < self.configuration.mutation_rate
                which_sites_mut = np.where(which_sites_mut)[0]
                if any(which_sites_mut):
                    for l, mutate_allele in enumerate(which_sites_mut):  
                    # added find = j; [:j]+[j:]
                        allele = chromosome.sequence[mutate_allele]
                        # got rid of deep copy
                        isAllele = self.nuc_alphabet.index(allele)
                        possible_mutations = self.nuc_alphabet[:isAllele]
                        possible_mutations += self.nuc_alphabet[isAllele + 1 :]
                        mutation_allele = self.rng.choice(possible_mutations, size=1)[0]
                        self.demes[i][j].sequence = (
                            self.demes[i][j].sequence[:mutate_allele]
                            + mutation_allele
                            + self.demes[i][j].sequence[mutate_allele + 1 :]
                        )
        return None


    def mutate_continuous(self):
        """
        """
        for i, deme_sequences in enumerate(self.demes):
            for j, chromosome in enumerate(self.deme_sequences):
                which_sites_mut = self.rng.uniform([0, 1, self.sequence_length])

                which_sites_mut = which_sites_mut < self.mutation_rate
                which_sites_mut = np.where(which_sites_mut)[0]
                if any(which_sites_mut):
                    for l, mutate_allele in enumerate(
                        which_sites_mut
                    ):  # added find = j; [:j]+[j:]
                        allele = chromosome.sequence[mutate_allele]
                        # got rid of deep copy
                        isAllele = self.nuc_alphabet.find(allele)
                        possible_mutations = self.nuc_alphabet[:isAllele]
                        possible_mutations = +self.nuc_alphabet[j + 1 :]
                        mutation_allele = self.rng.choice(possible_mutations, size=1)
                        demes[i][j].sequence = (
                            demes[i][j].sequence[:mutate_allele]
                            + mutation_allele
                            + demes[i][j].sequence[mutate_allele + 1 :]
                        )
                        
##    def set_analysis_files(self):?????????
    ############################################################################
    ### Simulate
    
    def simulate_multiple_trials_sequences(self):
        """
        """
        #loop through simK number_trials
        self.sequence_files = []  # move to multple trial function
        self.metadata_files = []  # move to multple trial function

        for trial in range(self.configuration.number_simulations):
            print("Trial/Simulation:",trial)
            self.trial = trial
            self.simulate_ancestral_deme_sequences()
        
        return None

    def simulate_ancestral_deme_sequences(self):
        """
        """
        self.nuc_alphabet = ["A", "T", "C", "G"]### global with __nuc_alphabet__??

        #self.sequence_files = []  # move to multple trial function
        #self.metadata_files = []  # move to multple trial function

        # Intialize properties #dynamic for multiple experiments
        self.set_starting_sequences()
        self.set_starting_demes_sequences()

        # EXPERIMENT
        for generation in range(self.configuration.number_generations):
            self.generation = generation
            temp_demes = []

            for i, temp_sequences in enumerate(self.demes):
                self.current_deme_i = i
                temp_draws = []
                # migration arrary
                migration_array = self.migration_matrix[i]
                for k in range(self.configuration.number_of_chromosomes):
                    sample_demes = range(self.configuration.number_of_demes)

                    from_deme_draw = self.rng.choice(
                        sample_demes, p=migration_array
                    ) 
                    
                    from_sequence_draw = self.rng.integers(
                        self.configuration.number_of_chromosomes
                    )

                    temp_draws.append(self.demes[from_deme_draw][from_sequence_draw])
                self.rng.shuffle(temp_draws)
                temp_demes.append(temp_draws)
            self.mutate_sequences()

            # update demes/pops with temp
            self.demes = temp_demes

            self.write_snap_shots()

        return None

    ############################################################################
    ### Write sequence traits to files
    # out = StringIO()
    # and write out silimarly
    # change write to as_....(fasta, csvtable,etc)
#    def write_snap_shots(self):
    def write_snap_shots(self):
        """
        """
        # change to runtime*****************************

        # print("Generation:", self.generation)
        if (self.generation in self.configuration.snapshot_times) | (
            self.generation == self.configuration.number_generations - 1
        ):
            (sequence_file, metadata_file, self.Resolution) = self.write_sequence_data()

            self.sequence_files.append(sequence_file)
            self.metadata_files.append(metadata_file)

        return None

    def write_sequence_data(self):
        """
        """

        # initialize array of strings to write
        csv_list = [",IDV,SUBPOP,ANC"]
        fasta_list = []
        for deme_count, deme in enumerate(self.demes):
            for chromose_count, chromosome in enumerate(deme):
                sequence_bases = chromosome.sequence

                fasta_list.append(
                    ">" + self.labels[deme_count] + ".chromsome." + str(chromose_count)
                )

                chromosome_count = (
                    self.labels[deme_count] + ".chromsome." + str(chromose_count)
                )

                ind_count_label = (
                    self.labels[deme_count] + "." + str(math.floor(chromose_count / 2))
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

        logging.info(
            json.dumps(self.configuration.outfile_prefix)
        )  # as same line term = ,
        logging.info(json.dumps(self.generation))

        # write deme data for analysis downstream
        # open file in write mode
        self.metadata_filename = self.configuration.outfile_prefix
        self.metadata_filename += "_"
        self.metadata_filename += str(self.generation)
        self.metadata_filename += "_"
        self.metadata_filename += str(self.trial)
        self.metadata_filename += ".csv"

        with open(self.metadata_filename, "w") as fp:
            for item in csv_list:
                # write each item on a new line
                fp.write("%s\n" % item)

        # Write deme sequences for analyis downstream
        # open file in write mode
        self.sequence_filename = self.configuration.outfile_prefix
        self.sequence_filename += "_"
        self.sequence_filename += str(self.generation)
        self.sequence_filename += "_"
        self.sequence_filename += str(self.trial)
        self.sequence_filename += ".fasta"
        with open(self.sequence_filename, "w") as fp:
            for item in fasta_list:
                # write each item on a new line
                fp.write("%s\n" % item)

        self.resolution = csv_list[0].split(",")[2]

        return (self.sequence_filename, self.metadata_filename, self.resolution)

    ############################################################################
    ### Initialize Continuous Model

    def set_starting_demes_continuous(self):
        """
        """
        #####
        for i in range(self.configuration.number_of_demes):
            self.labels.append(self.alphabet[i])
            self.demes[0].append(
                Continuous_trait_deme(
                    self.configuration.start_mean,
                    self.configuration.start_std,
                    self.alphabet[i],
                )
            )

        return None

    ############################################################################
    ## Simulate continuous traits

    def simulate_multiple_trials_continuous(self):
        """
        """
        self.csv_header = "generation,"
        self.csv_header += ",".join(self.labels)
        #loop through simK number_trials
        self.continuous_trial_files = []  # move to multple trial function

        for trial in range(self.configuration.number_simulations):
            self.trial = trial
            print("Trial/Simulation:",trial)
            self.simulate_deme_continuous_trait()
        
        return None        
        
    def set_migration_network_demes(self):
     
        if self.configuration.geometry == 'complete graph':

            self.migration_network = self.demes_mean[:self.current_deme_i] + self.demes_mean[self.current_deme_i + 1 :]
       
        else: #chain grapsh

            if self.current_deme_i == 0: 
                                #boundary case
                self.migration_network = [self.demes_mean[1]]
            
            elif self.current_deme_i == (len(self.demes_mean)-1): 
                                #boundary case
                self.migration_network = [self.demes_mean[-2]]
                
            else: 

                # chain case
                self.migration_network = [self.demes_mean[self.current_deme_i-1],self.demes_mean[(self.current_deme_i + 1)]]
            
        
        return None

    def simulate_deme_continuous_trait(self):
        """
        """
        self.demes = [[]]  # * self.configuration.number_of_demes
        self.set_starting_demes_continuous()
        
        self.continuous_mean_list = [self.csv_header]
        self.continuous_std_list = [self.csv_header]

        #set first row in csv
        demes_mean = []
        demes_std = []
        for i, deme_traits in enumerate(self.demes[0]):
            demes_mean.append(deme_traits.mean)
            demes_std.append(deme_traits.std)

        self.continuous_mean_list.append(
            str(-1) + "," + "".join(str(demes_mean))[1:-1]
        )
        self.continuous_std_list.append(
            str(-1) + "," + "".join(str(demes_std))[1:-1]
        )
        #end set first row

        # EXPERIMENT
        for generation in range(self.configuration.number_generations):
            self.generation = generation 

            self.demes_mean = []
            self.demes_std = []
            for i, deme_traits in enumerate(self.demes[generation]):
                self.demes_mean.append(deme_traits.mean)
                self.demes_std.append(deme_traits.std)


            temp_demes_mean = [0] * self.configuration.number_of_demes
            temp_demes_std = [0] * self.configuration.number_of_demes

            temp_draws = [0.0] * self.configuration.number_of_demes

            for i, deme_traits in enumerate(self.demes[generation]):
                self.current_deme_i = i
                # migration arrary
                migration_array = self.migration_matrix[i]

                number_draw = int(
                    migration_array[i] * self.configuration.number_of_chromosomes
                )
                temp_draws[i] = self.rng.normal(
                    loc=deme_traits.mean, scale=deme_traits.std, size=number_draw
                )

                migration_array = migration_array[:i] + migration_array[i + 1 :]

                temp_std = demes_std[:i] + demes_std[i + 1 :]
                
                self.set_migration_network_demes()
                
                for j, meanj in enumerate(self.migration_network):
                #for j, meanj in enumerate(demes_mean[:i] + demes_mean[i + 1 :]):
                    number_draw = int(
                        migration_array[j] * self.configuration.number_of_chromosomes
                    )

                    np.append(
                        temp_draws[i],
                        self.rng.normal(loc=meanj, scale=temp_std[j], size=number_draw),
                    )

                temp_demes_mean[i] = round(np.mean(temp_draws[i]), 2)
                temp_demes_std[i] = round(np.std(temp_draws[i]), 2)

            self.demes_mean = temp_demes_mean
            self.demes_std = temp_demes_std

            self.demes.append([])
            for i in range(self.configuration.number_of_demes):
                mean = self.demes_mean[i]
                std = self.demes_std[i]
                label = self.alphabet[i]

                self.demes[generation + 1].append(
                    Continuous_trait_deme(mean, std, label)
                )

            self.continuous_mean_list.append(
                str(generation) + "," + "".join(str(self.demes_mean))[1:-1]
            )
            self.continuous_std_list.append(
                str(generation) + "," + "".join(str(self.demes_std))[1:-1]
            )

        self.write_continuous_data()

        return None























    ############################################################################
    ## Write Continuous trait data to file
    # out = StringIO()
    # and write out silimarly
    # change write to as_....(fasta, csvtable,etc)

    def write_continuous_data(self):
        """
        """
        # write deme data for analysis downstream
        # open file in write mode

        self.continuous_mean_filename = self.configuration.outfile_prefix + "_" + str(self.trial) + "mean.csv"
        self.continuous_trial_files.append(self.continuous_mean_filename)
        with open(self.continuous_mean_filename, "w") as fp:
            for item in self.continuous_mean_list:
                # write each item on a new line
                fp.write("%s\n" % item)

        self.continuous_std_filename = self.configuration.outfile_prefix + "_" + str(self.trial) + "std.csv"
        self.continuous_trial_files.append(self.continuous_std_filename)
        with open(self.continuous_std_filename, "w") as fp:
            for item in self.continuous_std_list:
                # write each item on a new line
                fp.write("%s\n" % item)


        return None

#    def calc_continuous_entropy(self):
    def plot_continuous_entropy(self):
        """
        Clean up and format plots 
        """
        alpha = 0.5
        continuous_dict = {}
        continuous_mean_difference_dict = {}
        continuous_mean_difference_list = []
        
        means_dict = {}
        for i, filename in enumerate(self.continuous_trial_files[::2]): 
            means = pd.read_csv(filename)
            print(filename[:-8])
            means_dict[filename[:-8]] = []
            pairwise_difference_row = []
            for ii, meani in enumerate(means.iloc[-1].tolist()[1:]):
                print('mean',ii,meani)
                means_dict[filename[:-8]].append(meani)
        
        for i, filename in enumerate(self.continuous_trial_files[1::2]): 
            std = pd.read_csv(filename)
            print(filename[:-7])
            
            timestamp_means = means_dict[filename[:-7]]
            cross_entropy_row = []
            for ii, stdi in enumerate(std.iloc[-1].tolist()[1:]):
                print(stdi)
                shannon_index = 0.5 
                log_factor = ((2 * np.pi)**0.5) * ( stdi ** 2)
                shannon_index += np.log(log_factor) #log2piexp()std^2
                print('H',shannon_index)
                for j, stdj in enumerate(std.iloc[-1].tolist()[1:]):
                    if ii < j:   
                    #if ii > j:     
                        print(stdi,stdj)
                        meani = timestamp_means[ii]
                        meanj = timestamp_means[j]
                        print(meani,meanj)  
                        varience_star_h = stdj**2 
                        varience_star_h += ((alpha - 1) * stdi**2)
                        factora = np.log(2*np.pi*stdj**2)
                        factorb = 1 / (1 - 0.99) * np.log(stdj**2/varience_star_h)
                        factorc = ((meani- meanj) **2)/(varience_star_h)
                        cross_entropy = 0.5 * (factora + factorb + factorc)
                        
                        print('h_alpha(f_1,f_2)',cross_entropy)
                        
                        cross_entropy_row.append(cross_entropy)
                        
        return None
            

        
    def plot_continuous(self):
        """
        Clean up and format plots 
        """
        continuous_dict = {}
        continuous_mean_difference_dict = {}
        continuous_mean_difference_list = []
        for i, filename in enumerate(self.continuous_trial_files[::2]): 
            means = pd.read_csv(filename)
                  
            
            pairwise_difference_row = []
            for ii, meani in enumerate(means.iloc[-1].tolist()[1:]):
                for j, meanj in enumerate(means.iloc[-1].tolist()[1:]):
                    if ii != j:            
                        pairwise_difference_row.append(abs(meani-meanj)) 
            
            continuous_mean_difference_dict[i] = pairwise_difference_row
           
            for difference in pairwise_difference_row:
                continuous_mean_difference_list.append(difference)
       
        #df = pd.DataFrame(continuous_mean_difference_dict)
        #print(df)
        #df.plot(kind='hist')
        #df.hist()
        #plt.show()
        print(self.configuration.outfile_prefix + "_abs.png")
        abs_pairwise_difference_plot_filename = self.configuration.outfile_prefix + "_abs.png"
        plt.hist(continuous_mean_difference_list, bins="auto")
        plt.xlabel("Change in mean by deme")
        plt.ylabel("Count")
        plt.savefig(abs_pairwise_difference_plot_filename,dpi=300,format='png')
        #plt.show()
        
        sqr_pairwise_difference_plot_filename = self.configuration.outfile_prefix + "_sqr.png"
        continuous_mean_square_difference_list = np.square(continuous_mean_difference_list)
        plt.hist(continuous_mean_square_difference_list, bins="auto")
        plt.xlabel("Square change in mean")
        plt.ylabel("Count")
        #plt.show()
        plt.savefig(sqr_pairwise_difference_plot_filename,dpi=300,format='png')
        
    
        return None

    def log_analysis():
        logging.info("Analyzing")  # to check
        return None
        
        
