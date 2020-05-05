#!/usr/bin/env bash

# Designed to be run within the docker container.

alignment_dir="/alignment_plots"
input_file="/cartographer/test_read.fasta"

alignment_dir="/development/lrma/tools/Cartographer/alignment_plots/"
input_file="/development/lrma/tools/Cartographer/test_data/first_50_reads.fasta"
input_file="/development/lrma/tools/Cartographer/test_data/first_1000_reads.fasta"

echo "Running test, please wait."

/usr/bin/time -v /cartographer.py -v --MOSAIC -r ${input_file} -s /cartographer/KNOWN_SEQUENCES.fasta -l /cartographer/layout_order.txt --barcode_length 18 --barcode_file /development/lrma/tools/Cartographer/test_data/BarcodeLookupTable.csv &> /cartographer.$(date +%Y%m%dT%H%M%S).log
#/usr/bin/time -v /cartographer.py -v --MOSAIC -r ${input_file} -s /cartographer/KNOWN_SEQUENCES.fasta -l /cartographer/layout_order.txt --barcode_length 18 &> /cartographer.$(date +%Y%m%dT%H%M%S).log
#/usr/bin/time -v /cartographer.py -v --MOSAIC -r ${input_file} -s /cartographer/KNOWN_SEQUENCES.fasta -l /cartographer/layout_order.txt &> /cartographer.$(date +%Y%m%dT%H%M%S).log
#/usr/bin/time -v /cartographer.py -v --MOSAIC -r ${input_file} -s /cartographer/KNOWN_SEQUENCES.fasta -l /cartographer/layout_order.txt --plot --plot_folder ${alignment_dir} &> /cartographer.$(date +%Y%m%dT%H%M%S).log
#/usr/bin/time -v /cartographer.py -v -r /cartographer/test_read.fasta -s /cartographer/KNOWN_SEQUENCES.fasta -l /cartographer/layout_order.txt &> /cartographer.$(date +%Y%m%dT%H%M%S).log
#/usr/bin/time -v /cartographer.py -v -r /cartographer/test_read2.fasta -s /cartographer/KNOWN_SEQUENCES_rc.fasta -l /cartographer/layout_order_rc.txt &> tee /cartographer.$(date +%Y%m%dT%H%M%S).log

echo "Test complete."
