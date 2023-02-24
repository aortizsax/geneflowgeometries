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

import os
import pathlib
import sys
import argparse

from geneflowgeometries import simulate
from geneflowgeometries import calculate

def main():
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument(
            "src_paths",
            action="store",
            nargs="+",
            metavar="FILE",
            help="Path to source file(s).")

    parser.add_argument(
        "-S",
        "--simulate-sequences-or-continous",
        action="store",
        default='sequences',
        help="Choose to simulate sequences or continous variable [default=%(default)s].",
    )
    parser.add_argument(
        "-G",
        "--Geometries",
        action="store",
        default='complete graph',
        help="Network/graph to simulate. Options: Complete, chain, or FILE [default=%(default)s].",
    )
    parser.add_argument(
        "-N",
        "--Number-of-chromosomes-per-deme",
        action="store",
        default=1000,
        help="Migration rate across demes [default=%(default)s].",
    )
   parser.add_argument(
        "-P",
        "--Number-of-chromosomes-per-invdividual",
        action="store",
        default=2,
        help="Migration rate across demes [default=%(default)s].",
    )
    parser.add_argument(
        "-K",
        "--Number-of-demes",
        action="store",
        default=9,
        help="Migration rate across demes [default=%(default)s].",
    ) 
    parser.add_argument(
        "-m",
        "--migration-rate",
        action="store",
        default=1,
        help="Migration rate across demes [default=%(default)s].",
    )
    parser.add_argument(
        "-R",
        "--restriction-rate",
        action="store",
        default=0,
        help="Restriction rate across demes [default=%(default)s].",
    )
    parser.add_argument(
        "-mut",
        "--mutation-rate",
        action="store",
        default=10^-4,
        help="Mutation rate for gene simulaiton [default=%(default)s].",
    )
   parser.add_argument(
        "-L",
        "--sequence-length",
        action="store",
        default=10^-4,
        help="Mutation rate for gene simulaiton [default=%(default)s].",
    )
   parser.add_argument(
        "-mean",
        "--mean",
        action="store",
        default=1,
        help="Starting mean for simulating continous varible[default=%(default)s].",
    )
   parser.add_argument(
        "-std",
        "--standard-deviation",
        action="store",
        default=1,
        help="Starting standard deviation for simulating continous varible[default=%(default)s].",
    )
    parser.add_argument(
        "-simT",
        "--simulation-time",
        action="store",
        default=1000,
        help="Number of trials/generations to run simulation for [default=%(default)s].",
    )
    parser.add_argument(
        "-id",
        "--log-id",
        action="store",
        default="log-DATE",
        help="Filename for log files [default=%(default)s].",
    )
    parser.add_argument(
            "-o", "--output-prefix",
            action="store",
            default="output",
            help="Prefix for output files [default=%(default)s].")
    args = parser.parse_args()
    print("Hello, world.")
    
    #make config class dictionary 
    
    simulate.ancestral_deme_sequences()
    
    

if __name__ == '__main__':
    main()
