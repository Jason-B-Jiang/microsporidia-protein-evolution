#!/bin/bash

################################################################################

# Goal: Format single-copy orthogroup sequences for hmmscan search on a
#       particular set of OrthoFinder outputs, and produces hmmscan outputs to
#       run on the OrthoFinder outputs

# Input: Directory containing OrthoFinder results for some microsporidia species
#        and yeast

################################################################################

# Parse in command line args + define global variables
OF_results=$1
results_dir=$2  # where hmmscan results for species will eventually be saved

sp_pair=$(echo $OF_results | grep -E -o '[[:upper:]]_[[:lower:]]{2,4}\.?_S_cere')
src_dir=$(pwd)  # base directory for scripts

function main {
    # Copy folder with single copy orthogroup sequences, so we can split each
    # orthogroup fasta file into individual fasta files for each ortholog
    cp -r ${OF_results}/Single_Copy_Ortho* ${OF_results}/SCO_seqs_split

    # Split all orthogroup fasta files and generate hmmscan commands for all
    # orthogroup orthologs
    prepare_hmmscan_files ${OF_results}/SCO_seqs_split
}

function prepare_hmmscan_files {
    for f in ${1}/OG*
    do
        # Move each orthogroup fasta into its own folder
        local og=$(basename $f)
        mkdir ${1}/${og%.fa}
        mv $f ${1}/${og%.fa}

        # Set up folder in results for outputting hmmscan results for orthogroup
        mkdir -p ${results_dir}/${sp_pair}/${og%.fa}/hmmscan

        # Split orthogroup fasta into separate fasta files for each ortholog
        cd ${1}/${og%.fa}
        awk '/^>/{s=++d".fa"} {print > s}' *.fa
        rm ./${og}  # remove original orthogroup fasta file

        # rename split files by awk to the names of their ortholog sequences
        mv 1.fa $(head -1 1.fa | sed -r 's/>//').fa
		mv 2.fa $(head -1 2.fa | sed -r 's/>//').fa

        cd $src_dir  # change back to script directory
    done
}

main