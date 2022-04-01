#!/bin/bash

################################################################################

# Goal: Using files set up by prep_hmmscan_cmds, create a text file of hmmscan
#       commands to run for all orthogroups in a species pair

################################################################################

# Parse in command line arguments
OF_results=$1
results_dir=$2
cmd_dir=$3
pfam_lib=$4

function main {
    # TODO - split into helper functions once you're done
    local sp_pair_regex='[[:upper:]]_[[:lower:]]{2,4}\.?_S_cere'
    local sp_pair=$(echo $OF_results | grep -E -o ${sp_pair_regex})
    echo $sp_pair

    # create necessary files and directories
    mkdir -p ${cmd_dir}  # create directory if not exist
    touch ${cmd_dir}/${sp_pair}

    for og in ${OF_results}/SCO_seqs_split/OG*
    do
        for seq in ${og}/*.fa
        do
            local og_name=$(basename ${og})
            local seq_name=$(basename ${seq})
            local output_f=${results_dir}/${og_name}/hmmscan/${seq_name%.fa}
            local hmmscan_cmd="hmmscan -o ${output_f} --cut_ga ${pfam_lib} ${seq}"
            echo ${hmmscan_cmd} >> ${cmd_dir}/${sp_pair}
        done
    done
}

main