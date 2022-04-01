#!/bin/bash

################################################################################

# Goal: Use cath-resolve-hits to resolve overlapping hmmscan hits

# Input: $1 to n - 1: All orthogroup architectures
#        $n: cath-resolve-hits executable

################################################################################

args=( "$@" )
crh=${args[-1]}

# iterate up until last argument, which is the crh executable
length=$(expr ${#args[@]} - 1)

for (( i=0; i<${length}; i++ ));
do
    # create directory for resolved hits by cath-resolve-hits
    mkdir ${args[$i]}/crh

    # run CRH on hmmscan output for each ortholog in the orthogroup
    for seq in ${args[$i]}/hmmscan/*
    do
        seq_name=$(basename $seq)  # get ortholog name from file path

        # run cath-resolve-hits with lenient bit score cutoff of 0.1, to allow for
        # low scoring Pfam model hits that passed gathering thresholds
        ${crh} --input-format hmmscan_out --worst-permissible-bitscore 0.1 \
        ${seq} > ${args[$i]}/crh/${seq_name}
    done

done