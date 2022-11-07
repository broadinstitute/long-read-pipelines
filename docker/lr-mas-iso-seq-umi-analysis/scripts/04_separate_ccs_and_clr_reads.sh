#!/usr/bin/env bash

in_file=$1

if [[ $# -ne 1 ]] ;  then
   echo "ERROR: Must provide an input file and an output prefix." 1>&2 
   exit 1
fi

out_prefix=$( echo $1 | sed 's@.bam$@@' )

MIN_CCS_RQ="0.0"

set -e

t_start=$(date +%s.%N)
echo "Starting CCS / CLR filtration."
bamtools filter -tag "rq":">=${MIN_CCS_RQ}" -in ${in_file} -out ${out_prefix}.ccs.bam &
bamtools filter -tag "rq":"<${MIN_CCS_RQ}" -in ${in_file} -out ${out_prefix}.clr.bam &

echo "Waiting for completion."
wait
t_end=$(date +%s.%N)
t_elapsed=$( echo "scale=4;${t_end} - ${t_start}" | bc )

echo 'Done!'
echo "Elapsed time: ${t_elapsed}"

