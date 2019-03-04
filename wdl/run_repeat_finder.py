#!/usr/bin/env python3

import os
os.chdir(os.path.abspath(os.path.dirname(__file__)))

import argparse
import logging

from run_cromwell_job_utils import execute

logging.basicConfig(level=logging.INFO)

FILE_PATH_PREFIX = {
    "local": "/Users/weisburd/project1/",
    "remote": "gs://broad-dsp-spec-ops/scratch/weisburd/",
}

WORKFLOWS = {
    "repeat-finder": { "name": "RepeatFinderWorkflow", "wdl": "repeat-finder/repeat-finder.wdl" },
}


DATA = {
    "hg19_chr22": {
        "ref_fasta_gz_urls": [
            "https://raw.githubusercontent.com/bw2/str-callers/master/wdl/repeat-finder/hg19_chr22.fa.gz",
        ],
    },
    # URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq/; for i in $(curl -s $URL | cut -c 57- | grep ref | grep fa.gz | grep -v mfa | grep -v unloc | grep -v unplace | grep -v alts | sort ) ; do echo \"${URL}$i\", | sed s/ftp:/http:/; done
    "hg19": {
        "ref_fasta_gz_urls": [
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr1.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr10.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr11.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr12.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr13.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr14.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr15.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr16.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr17.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr18.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr19.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr2.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr20.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr21.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr22.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr3.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr4.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr5.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr6.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr7.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr8.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chr9.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chrMT.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chrX.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq//hs_ref_GRCh37.p5_chrY.fa.gz",
        ],
    },
    # URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq/; for i in $(curl -s $URL | cut -c 57- | grep ref | grep fa.gz | grep -v mfa | grep -v unloc | grep -v unplace | grep -v alts | sort ) ; do echo \"${URL}$i\", | sed s/ftp:/http:/; done
    "hg38": {
        "ref_fasta_gz_urls": [
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr1.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr10.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr11.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr12.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr13.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr14.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr15.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr16.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr17.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr18.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr19.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr2.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr20.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr21.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr22.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr3.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr4.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr5.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr6.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr7.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr8.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chr9.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chrMT.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chrX.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq//hs_ref_GRCh38.p12_chrY.fa.gz",
        ],
    },
    # URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/; for i in $(curl -s $URL | cut -c 57- | grep ref | grep fa.gz | grep -v mfa | grep -v unloc | grep -v unplace | grep -v alts | sort ) ; do echo \"${URL}$i\", | sed s/ftp:/http:/; done
    "mm38": {
        "ref_fasta_gz_urls": [
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr1.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr10.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr11.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr12.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr13.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr14.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr15.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr16.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr17.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr18.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr19.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr2.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr3.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr4.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr5.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr6.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr7.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr8.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr9.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chrMT.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chrX.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chrY.fa.gz",
        ],
    },
    # URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/; for i in $(curl -s $URL | cut -c 57- | grep ref | grep fa.gz | grep -v mfa | grep -v unloc | grep -v unplace | grep -v alts | sort ) ; do echo \"${URL}$i\", | sed s/ftp:/http:/; done
    "pan_troglodytes": {
        "ref_fasta_gz_urls": [
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr1.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr10.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr11.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr12.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr13.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr14.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr15.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr16.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr17.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr18.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr19.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr20.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr21.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr22.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr2A.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr2B.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr3.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr4.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr5.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr6.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr7.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr8.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chr9.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chrMT.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chrX.fa.gz",
            "http://ftp.ncbi.nlm.nih.gov/genomes/Pan_troglodytes/Assembled_chromosomes/seq/ptr_ref_Clint_PTRv2_chrY.fa.gz",
        ],
    },
}

OTHER_PARAMS = {
    "default":  { "match_score": 2, "mismatch_score": 7, "indel_score": 7, "minscore": 50, },
    "gangstr1": { "match_score": 2, "mismatch_score": 5, "indel_score": 17, "minscore": 24, },
    "gangstr2": { "match_score": 2, "mismatch_score": 3, "indel_score": 5, "minscore": 24, },  # from https://github.com/gymreklab/GangSTR/blob/master/reference/chr_make_reference.sh#L51-L55
    #"lenient1": { "match_score": 2, "mismatch_score": 3, "indel_score": 5, "minscore": 8, },
    "custom1": { "match_score": 2, "mismatch_score": 5, "indel_score": 17, "minscore": 8, },
    "custom2": { "match_score": 2, "mismatch_score": 19, "indel_score": 21, "minscore": 8, },
    "custom3": { "match_score": 2, "mismatch_score": 31, "indel_score": 31, "minscore": 8, },
    "custom4": { "match_score": 2, "mismatch_score": 5, "indel_score": 17, "minscore": 12, },
    #"custom5": { "match_score": 2, "mismatch_score": 5, "indel_score": 17, "minscore": 24, },  # same as gangstr1
    "custom6": { "match_score": 2, "mismatch_score": 1000000, "indel_score": 1000000, "minscore": 8, },
    "custom7": { "match_score": 2, "mismatch_score": 7, "indel_score": 7, "minscore": 24, },
    "custom8": { "match_score": 2, "mismatch_score": 7, "indel_score": 7, "minscore": 22, },  # lowest minscore cutoff in willems 2014
    #"lenient5": { "match_score": 2, "mismatch_score": 51, "indel_score": 51, "minscore": 8, },
}

#LOCAL_WORKFLOW_OPTIONS = {}
GCS_WORKFLOW_OPTIONS = {
    "write_to_cache": True,
    "read_from_cache": True
}

def parse_args():
    p = argparse.ArgumentParser()

    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--local", action="store_true")
    group.add_argument("--remote", action="store_true")

    p.add_argument("--skip-input-validation", action="store_true", help="Skip input validation.")
    p.add_argument("--test", action="store_true", help="Don't actually submit the job.")

    p.add_argument("tool", choices=WORKFLOWS.keys())
    p.add_argument("data", nargs="+", choices=DATA.keys())
    p.add_argument("params", choices=OTHER_PARAMS.keys())
    #p.add_argument("label", nargs="?", help="optional label describing this run")

    args = p.parse_args()

    return args



def main():
    args = parse_args()

    tool_name = WORKFLOWS[args.tool]['name']
    wdl_path = WORKFLOWS[args.tool]['wdl']
    workflow_options = GCS_WORKFLOW_OPTIONS
    data_root_path = FILE_PATH_PREFIX["local"] if args.local else FILE_PATH_PREFIX["remote"]

    # generate input json
    for data in args.data:
        input_json = {}
        input_json.update(DATA[data])
        input_json.update(OTHER_PARAMS[args.params])
        input_json["output_prefix"] = f"{args.tool}_{data}__{args.params}__" + '__'.join([f"{key.replace('_score', '')}{value}" for key, value in OTHER_PARAMS[args.params].items()])

        execute(args, tool_name, data, wdl_path, workflow_options, input_json, data_root_path)


if __name__ == "__main__":
    main()
