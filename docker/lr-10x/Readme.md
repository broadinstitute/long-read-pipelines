# 10x Barcode Extraction Tool

This tool can be used for identifying 10x adapter sequences, extracting and correcting barcode and UMI information, and annotating reads accordingly.

## Installation
This tool requires the packages specified in the conda `environment.yml` file found in the `docker` subdirectory, as well as the shared library [Complete-Striped-Smith-Waterman-Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) and a [starcode](https://github.com/gui11aume/starcode) binary executable. Please refer to the linked github repositories for installing these requirements.

## Usage
```
usage: tool.py [-h] -b BAM -a ADAPTER -r REVERSE_ADAPTER -n NAME
               [--whitelist-10x WHITELIST_10X]
               [--whitelist-illumina WHITELIST_ILLUMINA]
               [--max-reads MAX_READS] [--contig CONTIG]
               [--read-end-length READ_END_LENGTH] [--record-umis]
               [--ssw-path SSW_PATH] [--starcode-path STARCODE_PATH]

Reads an input BAM file and tries to extract adapter and barcode sequences.
When found, annotated reads are written to an output file, if an output
filename is provided

optional arguments:
  -h, --help            show this help message and exit
  --whitelist-10x WHITELIST_10X
                        10x whitelist filename. This may be GZIP compressed
                        (has to have extension .gz in that case)
  --whitelist-illumina WHITELIST_ILLUMINA
                        Illumina whitelist filename. This may be GZIP
                        compressed (has to have extension .gz in that case)
  --max-reads MAX_READS
                        Number of reads after which the processing should be
                        terminated
  --contig CONTIG       Perform analysis only on this contig
  --read-end-length READ_END_LENGTH
                        Interval from both ends of the read in which to search
                        for the adapter sequence
  --record-umis         If enabled, all barcodes and UMIs will be written to
                        file. This increases memory usage.
  --ssw-path SSW_PATH   Path to the Striped Smith-Waterman library
  --starcode-path STARCODE_PATH
                        Path to the starcode executable

required named arguments:
  -b BAM, --bam BAM     BAM filename
  -a ADAPTER, --adapter ADAPTER
                        Adapter FASTA filename. BWA index must be present.
  -r REVERSE_ADAPTER, --reverse-adapter REVERSE_ADAPTER
                        Reverse adapter FASTA filename. BWA index must be
                        present.
  -n NAME, --name NAME  Analysis name (output prefix)

Two barcode whitelists can be provided: a 10x whitelist and an Illumina
whitelist. If both whitelists are provided, reads will only be annotated if
the barcode matches the Illumina list. If only a 10x whitelist is provided,
reads will be annotated if the barcode matches the 10x list. In any case,
reads will be annotated with the "raw barcode" tag (CR). If no whitelist is
provided, all reads will be annotated with the barcode tage (CB) after
correction.
```