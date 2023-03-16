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


# Classes
class Chromosome:
    def __init__(self, sequence, ancenstral_deme):
        self.sequence = sequence
        self.ancenstral_deme = ancenstral_deme

    def __str__(self):
        return self.sequence + self.ancenstral_deme

class Continuous_trait_deme:
    def __init__(self, mean, std,deme):
        self.mean = mean
        self.std = std
        self.deme = deme

    def __str__(self):
        return "Deme:" + self.deme+ ",Mean:"+self.mean+',Std:'+self.std
        
        
alphabet = string.ascii_lowercase


def continuous_trait_evoluion(
    Geometry,
    number_of_chromosomes,
    number_of_demes,
    migration_rate,
    number_generations,
    start_mean,
    start_std,
    ):
    
    FASTAFILENAMES = []  # csv means, std
    CSVFILENAMES = []  # pop data
    demes = [[]]
    labels = []
    for i in range(number_of_demes):
        labels.append(alphabet[i])
        demes[0].append(Continuous_trait_deme(start_mean,start_std,alphabet[i]))

    csv_list = [',generation,A,B,C,D']#add labels
    
    
    #EXPERIMENT
    for trial in range(number_generations):
        print('trial',trial)
        temp_demes_mean = [0] * number_of_demes
        temp_demes_std = [0] * number_of_demes
                # simulate draws foreward in time
        
        temp_draws = [0.0] * number_of_demes
        
        demes_mean = []
        demes_std = []
        print(demes[trial])
        for i, deme_traits in enumerate(demes[trial]):
            print('i',i)
            demes_mean.append(deme_traits.mean)
            demes_std.append(deme_traits.std)
            
        for i, deme_traits in enumerate(demes[trial]):
        
            
            # migration arrary
            migration_array = [migration_rate] * number_of_demes
            migration_array[i] = 1 - (migration_rate * (number_of_demes - 1))
            print(migration_array)
            number_draw = int(migration_array[i] * number_of_chromosomes)
            temp_draws[i] = np.random.normal(loc=deme_traits.mean,
                                             scale=deme_traits.std,
                                             size=number_draw)
            migration_array = migration_array[:i]+migration_array[i+1:] 
            print(migration_array)
            #if sum(migration array    
            temp_std = demes_std[:i]+demes_std[i+1:]       
            for j, meanj in enumerate(demes_mean[:i]+demes_mean[i+1:]):
                print('j',j)

                number_draw = int(migration_array[j] * number_of_chromosomes)
                np.append(temp_draws[i],
                          np.random.normal(loc=meanj,
                                           scale=temp_std[j],
                                           size=number_draw))
            #else pass
            temp_demes_mean[i] = round(np.mean(temp_draws[i]),2)
            temp_demes_std[i] = round(np.std(temp_draws[i]),2)

        demes_mean = temp_demes_mean
        demes_std = temp_demes_std
        demes.append([])
        for i in range(number_of_demes): 
            print('next trial',trial+1)
            mean = demes_mean[i]
            std = demes_std[i]
            label = alphabet[i]
            demes[trial+1].append(Continuous_trait_deme(mean,std,label))
        
        csv_list.append(str(trial)+','+''.join(str(demes_mean))[1:-1])
    print(csv_list)
    
    
    # write deme data for analysis downstream
    # open file in write mode
    date = datetime.datetime.now()
    date = str(date).split(" ")[0]
    FILENAME = date+"_"+'continuous_'+Geometry+"_"+str(number_of_chromosomes)+"_"
    FILENAME+= str(number_of_demes)+"_"+str(migration_rate)+"_"+str(number_generations)
    CSVFILENAME = "./completegraph_simulation_continuous_trait_evolution.csv"
    CSVFILENAME = FILENAME + "_" + str(trial) + ".csv"
    with open(CSVFILENAME, "w") as fp:
        for item in csv_list:
            # write each item on a new line
            fp.write("%s\n" % item)
        print("Done", trial)

        CSVFILENAMES.append(CSVFILENAME)

    return CSVFILENAMES
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def continuous_trait_evoluion_old(
    Geometry,
    number_of_chromosomes,
    number_of_demes,
    migration_rate,
    number_generations
    ):
    FASTAFILENAMES = []  # csv means
    CSVFILENAMES = []  # pop data
    demes_mean = [[]] * number_of_demes
    demes_std = [[]] * number_of_demes
    labels = []
    for i, mean in enumerate(demes_mean):
        labels.append(alphabet[i])

    csv_list = [',generation,A,B,C,D']#add labels


    ###
    for i, mean in enumerate(demes_mean):
        demes_mean[i] = [0] 
        demes_std[i] = [1] 

    
    

    for trial in range(number_generations):
        temp_demes_mean = [0] * number_of_demes
        temp_demes_std = [0] * number_of_demes
                # simulate draws foreward in time
        
        temp_draws = [0.0] * number_of_demes
        for i, mean in enumerate(demes_mean):
        
            
            # migration arrary
            migration_array = [migration_rate] * number_of_demes
            migration_array[i] = 1 - (migration_rate * (number_of_demes - 1))
            number_draw = int(migration_array[i] * number_of_chromosomes)
            temp_draws[i] = np.random.normal(loc=mean,
                                             scale=demes_std[i],
                                             size=number_draw)
            migration_array = migration_array[:i]+migration_array[i+1:]            
            for j, meanj in enumerate(demes_mean[:i]+demes_mean[i+1:]):

                number_draw = int(migration_array[j] * number_of_chromosomes)
                np.append(temp_draws[i],
                          np.random.normal(loc=mean,
                                           scale=demes_std[i],
                                           size=number_draw))
        
            temp_demes_mean[i] = round(np.mean(temp_draws[i]),2)
            temp_demes_std[i] = round(np.std(temp_draws[i]),2)

        demes_mean = temp_demes_mean
        demes_std = temp_demes_std
        
        csv_list.append(str(trial)+','+''.join(str(demes_mean))[1:-1])
    print(csv_list)
    # write deme data for analysis downstream
    # open file in write mode
    CSVFILENAME = "./completegraph_simulation_continuous_trait_evolution.csv"
    with open(CSVFILENAME, "w") as fp:
        for item in csv_list:
            # write each item on a new line
            fp.write("%s\n" % item)
        print("Done", trial)

        CSVFILENAMES.append(CSVFILENAME)

    return CSVFILENAMES
