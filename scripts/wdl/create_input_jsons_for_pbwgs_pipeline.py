from __future__ import print_function
import argparse
from google.cloud import storage
import xml.etree.ElementTree as etree
import re
import sys
import os


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


parser = argparse.ArgumentParser(description='Create an input JSON for the LRWGS pipeline')
parser.add_argument('gcs_roots', metavar='G', type=str, nargs='+', help='GCS path(s)')
args = parser.parse_args()

delimiter = '/'
quoted_paths = []

for gcs_root in args.gcs_roots:
    quoted_paths.append(f'"{gcs_root}"')

    gcs_path = gcs_root.replace("gs://", "").split(delimiter)
    bucket_name = gcs_path[0]
    prefix = delimiter.join(gcs_path[1:])

    storage_client = storage.Client()
    bam_blobs = storage_client.list_blobs(bucket_name, prefix=prefix)

    for bam_blob in bam_blobs:
        if bam_blob.name.endswith(".subreads.bam") and not bam_blob.name.endswith(".scraps.bam"):
            gsdir = os.path.dirname(bam_blob.name)
            fc_name = re.sub(".subreads.bam", "", os.path.basename(bam_blob.name))

            md_blobs = storage_client.list_blobs(bucket_name, prefix=gsdir)
            for md_blob in md_blobs:
                if md_blob.name.endswith(".metadata.xml") and not md_blob.name.endswith(".run.metadata.xml"):
                    md_blob.download_to_filename("metadata.xml")
                    a = etree.parse("metadata.xml")
                    info = {}

                    for e in a.iter():
                        if bool(e.attrib):
                            tag = re.sub(r"{.+}", "", e.tag)
                            for k, v in e.attrib.items():
                                key = tag + "." + k
                                info[key] = "".join(v)

                    json_path = "data/PBCCSGenomeSingleSample/{fc_name}.json"
                    sample_name = fc_name
                    if 'BioSample.Name' in info:
                        sample_name = info['BioSample.Name']
                    elif 'WellSample.Name' in info:
                        sample_name = info['WellSample.Name']

                    sample_name = re.sub("[ #]", "_", sample_name)

                    print(info)

                    is_ccs = "unknown"
                    if "IsCCS" in info:
                        is_ccs = info['WellSample.IsCCS']

                    print(f'{bam_blob.name} {sample_name} {fc_name} {is_ccs}')

                    s = f"""{{
    "PBCCSWholeGenomeSingleFlowcell.gcs_input_dir":     "gs://{bucket_name}/{gsdir}",
    "PBCCSWholeGenomeSingleFlowcell.sample_name":       "{sample_name}",

    "PBCCSWholeGenomeSingleFlowcell.ref_fasta":         "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa",
    "PBCCSWholeGenomeSingleFlowcell.ref_fasta_fai":     "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai",
    "PBCCSWholeGenomeSingleFlowcell.ref_dict":          "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict",

    "PBCCSWholeGenomeSingleFlowcell.tandem_repeat_bed": "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/human_GRCh38_no_alt_analysis_set.trf.bed",
    "PBCCSWholeGenomeSingleFlowcell.ref_flat":          "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/refFlat.txt",
    "PBCCSWholeGenomeSingleFlowcell.dbsnp_vcf":         "gs://broad-dsde-methods-long-reads/resources/var_db/dbsnp/v151_GRCh38p7/common_all_20180418.chr.prefix.vcf.gz",
    "PBCCSWholeGenomeSingleFlowcell.dbsnp_tbi":         "gs://broad-dsde-methods-long-reads/resources/var_db/dbsnp/v151_GRCh38p7/common_all_20180418.chr.prefix.vcf.gz.tbi",

    "PBCCSWholeGenomeSingleFlowcell.mt_chr_name":       "chrM",
    "PBCCSWholeGenomeSingleFlowcell.metrics_locus":     "gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/metric_intervals.interval_list",

    "PBCCSWholeGenomeSingleFlowcell.gcs_out_root_dir":  "gs://broad-dsde-methods-long-reads-outgoing/PBCCSWholeGenomeSingleFlowcell"
}}"""

                    #json_path = f'data/PBCCSWholeGenomeSingleFlowcell/{sample_name}.{fc_name}.json'
                    #json_file = open(json_path, 'w')
                    #print(s, file=json_file)
                    #json_file.close()
