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

#
# Simulation
#
###########################################################################
## Reader Implementation

# Libraries
# import modules front end parser
from nyemtaay.parse.parser import read_fasta_files, read_metadata, to_dataframe
from nyemtaay.calculate import populationgeneticstats, informationtheory
from nyemtaay.mathlib import sterling
from nyemtaay.tests.nuetrality import tajimas_d
import logging


def parse(filenames):
    # use parser modules
    # pass list of fasta files to fasta parser
    fastafiles = [filenames[0]]
    sequence_matrix = read_fasta_files(fastafiles)

    # pass metadata to its parser
    metadata = filenames[1]
    data_matrix = read_metadata(metadata, 0)

    # convert matrix to dataframe with indexes matching metadata
    sequence_dataframe = to_dataframe(sequence_matrix, data_matrix)
    print("Done parsing")

    return (sequence_dataframe, data_matrix)


def nei_fst(sequences, metadata, deme_identifier,geo):
    populationgeneticstats.nei_fst(sequences, metadata, deme_identifier,geo)
    logging.info("Analyzing nei fst for gen")
    return None


def by_deme_pairwise_fst(sequences, metadata, deme_identifier,geo):
    populationgeneticstats.by_deme_pairwise_fst(sequences, metadata, deme_identifier,geo)

    return None


def weir_cockerman_fst():
    populationgeneticstats.weir_cockerman_fst(sequences, metadata, deme_identifier,geo)

    return None


def weir_goudet_population_specific_fst(sequences, metadata, deme_identifier,geo):
    populationgeneticstats.weir_goudet_population_specific_fst(
        sequences, metadata, deme_identifier,geo
    )

    return None
    #


def wright_fis(sequences, metadata, deme_identifier,geo):
    populationgeneticstats.wright_fis(sequences, metadata, deme_identifier,geo)

    return None
    
def information_theroy_pipeline(sequences, metadata, deme_identifier,geo):
    (gamete_probabilities, population_dict,pos_allele_probablity_deme_dict) = informationtheory.sequences_to_gamete_prob(sequences, metadata, deme_identifier,geo)
    #print(gamete_probabilities)
    
    #by gamete
    shannon_entropy_by_deme = informationtheory.demes_shannon_entropy(gamete_probabilities)
    
    jsd_by_comp = informationtheory.demes_jsd(gamete_probabilities, population_dict)
#    
#    norm_jsd_by_comp = informationtheory.demes_norm_jsd(gamete_probabilities, population_dict)
#    I_by_comp = informationtheory.demes_info_flow_direction(gamete_probabilities, population_dict)
#    
#    (gamete_probabilities, population_dict) = informationtheory.sequences_to_random_deme_combinations(sequences, metadata, deme_identifier,geo)
#    
#    IR_dict_by_comparison = informationtheory.randomized_information_flow_directionality(gamete_probabilities, I_by_comp)
#    
    
    #by posisitional allele
    informationtheory.demes_allele_shannon_entropy(pos_allele_probablity_deme_dict)
    informationtheory.demes_allele_jsd(pos_allele_probablity_deme_dict, population_dict)
    
    #look at tp 20 geneflowgeometries -c ../tests/config/sequenceconfig.ini -simK 1 -G "chain graph" -R 0.5 -m 0.25 -k 3 -mut 0.0001 -N 20 -L 10 -simT 200
    
    return None


def plot_fst(fst_network_dataframe):
    #networkx target, ..., (un)directed, weight

    return None

def plot_ancenstal_snapshots():

    return None


def plot_():

    return None


# WRIGHT Fst

# Pairwise Fst

# Weir and Cockerham Fst Population specific

# anc bar graph

# network x plot of calculation
# topo
# completegraph()
# build chain int until same size
