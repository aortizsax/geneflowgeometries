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

import datetime
import math
import logging
import json
#from geneflowgeometries.log.multi_handlers_logger import setup_logger


def write_log_analysis():
    logging.info("Analyzing")
    return 


def write_log_continuous(csv_list, outfile_prefix):
    # write deme data for analysis downstream
    # open file in write mode
    date = datetime.datetime.now()
    date = str(date).split(" ")[0]

    CSVFILENAME = outfile_prefix + ".csv"
    with open(CSVFILENAME, "w") as fp:
        for item in csv_list:
            # write each item on a new line
            fp.write("%s\n" % item)
        print("Done", CSVFILENAME)

        CSVFILENAME

    return CSVFILENAME


def write_log_sequences(demes, outfile_prefix, generation, labels, config_dict):
    print("snapshot", generation)
    taxon_labels = []
    csv_list = [",IDV,SUBPOP,ANC"]
    fasta_list = []
    for deme_count, deme in enumerate(demes):
        for chromose_count, chromosome in enumerate(deme):
            sequence_bases = chromosome.sequence

            fasta_list.append(
                ">" + labels[deme_count] + ".chromsome." + str(chromose_count)
            )

            chromosome_count = labels[deme_count] + ".chromsome." + str(chromose_count)

            ind_count_label = (
                labels[deme_count] + "|" + str(math.floor(chromose_count / 2))
            )
            # change to number ploidy

            csv_to_append = ",".join(
                [
                    chromosome_count,
                    ind_count_label,
                    labels[deme_count],
                    chromosome.ancenstral_deme,
                ]
            )
            csv_list.append(csv_to_append)

            fasta_list.append(sequence_bases)
        print()

    # write deme data for analysis downstream
    # open file in write mode
    date = datetime.datetime.now()  #pull up in scripts
    date = str(date).split(" ")[0]

    FILENAME = outfile_prefix
    CSVFILENAME = date + FILENAME + "_" + str(generation) + ".csv"

    logging.info(json.dumps(outfile_prefix))
    logging.info(json.dumps(generation))
    logging.info(json.dumps(config_dict))

    with open(CSVFILENAME, "w") as fp:
        for item in csv_list:
            # write each item on a new line
            fp.write("%s\n" % item)
        print("Done", generation)

    # Write deme sequences for analyis downstream
    # open file in write mode
    FASTAFILENAME = date + FILENAME + "_" + str(generation) + ".fasta"
    with open(FASTAFILENAME, "w") as fp:
        for item in fasta_list:
            # write each item on a new line
            fp.write("%s\n" % item)
        print("Done", generation)

    Resolution = csv_list[0].split(",")[2]

    return (FASTAFILENAME, CSVFILENAME, Resolution)
