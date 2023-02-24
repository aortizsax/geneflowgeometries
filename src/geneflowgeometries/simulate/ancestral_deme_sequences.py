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
## DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN BE LIABLE FOR ANY DIRECT,
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

import random 



alphabet_code = {'A':'0','C':'1','G':'2','T':'3'}
code_alphabet = {'0':'A', '1':'C','2':'G','3':'T'}
alphabet = string.ascii_lowercase


def ancestral_deme_sequences(Geometry, number_of_chromosomes,
                             number_of_ploidy, number_of_demes, migration_rate, 
                             restriction_rate, mutation_rate, sequence_length,
                             simT)
    
    demes = [[]] * number_of_demes
    labels=[]
    for i, deme_sequences in enumerate(demes):
        labels.append(alphabet[i])
        
    for i, deme_sequences in enumerate(demes):
        seq = ''
        for j in range(length_sequence):
            seq += str(random.randint(0,3))
        demes[i] = [seq] * N
    for seq in demes:
        print(seq[0])
        
    for trial in range(number_generations):
        temp_demes = []
        for i, deme_sequences in enumerate(demes):
            temp_pops.append(deme_sequences.copy())
        
        #simulate coal foreward in time
        for i, temp_sequences in enumerate(temp_pops):
            temp_draws = []
            for j, deme_sequences in enumerate(demes):
                for k in range(int(N / number_cells)):
                    draw = random.randint(0,N-1)
                    temp_draws.append(pops[j][draw])
                random.shuffle(temp_draws)
                temp_pops[i] = temp_draws

        #update demes/pops with temp
        pops = temp_pops
        
        
        #mutate function
        for i, deme_sequences in enumerate(demes):
            for j,sequence  in enumerate(deme_sequences):
                for k, allele in enumerate(sequence):
                    if random.randint(0,int(3*mu/4))==0:
                        print(i,j,demes[i][j])
                        print(demes[i][j][k],demes[i][j][:k], str(random.randint(0,3)), demes[i][j][k+1:])
                        demes[i][j] = demes[i][j][:k] + str(random.randint(0,3)) + demes[i][j][k+1:]
                        print(len(demes[i][j]))


    taxon_labels = [] 
    csv_list = [',IDV,SUBPOP']
    fasta_list = []
    for deme_count,deme in enumerate(demes):
        
        for chromose_count, sequence in enumerate(deme):
            sequence_bases = ''

            for base in sequence:
                sequence_bases += code_alphabet[base]

            fasta_list.append('>'+str(deme_count)
                              +'.chromsome.'+str(chromose_count))
            
            csv_list.append(str(deme_count)
                              +'.chromsome.'+str(chromose_count)+
                  ','+str(pop_count)+'.'+
                  str(math.floor(chromose_count/2))+
                           ','+str(pop_count))
            
            fasta_list.append(sequence_bases)
        print( )

    # write deme data for analysis downstream
    # open file in write mode
    with open('./drift_simulation.csv', 'w') as fp:
        for item in csv_list:
            # write each item on a new line
            fp.write("%s\n" % item)
        print('Done')
        

    # Write deme sequences for analyis downstream
    # open file in write mode
    with open('./drift_simulation.fasta', 'w') as fp:
        for item in fasta_list:
            # write each item on a new line
            fp.write("%s\n" % item)
        print('Done')
