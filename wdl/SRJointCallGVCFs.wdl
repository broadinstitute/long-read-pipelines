version 1.0

#############################################################################################################
## A workflow that performs joint calling on single-sample gVCFs from GATK4 HaplotypeCaller using GenomicsDB.
#############################################################################################################

import "SRJointGenotyping.wdl" as SRJOINT
import "tasks/Finalize.wdl" as FF

workflow SRJointCallGVCFs {
    input {
        Array[File] gvcfs
        Array[File] tbis

        File ref_map_file

        File interval_list

        String prefix

        String gcs_out_root_dir
    }

    parameter_meta {
        gvcfs:            "GCS paths to gVCF files"
        tbis:             "GCS paths to gVCF tbi files"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        prefix:           "prefix for output joint-called gVCF and tabix index"
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/SRJointCallGVCFs/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    # Create sample-name map:
    call SRJOINT.CreateSampleNameMap as CreateSampleNameMap {
        input:
            gvcfs = gvcfs,
            prefix = prefix
    }

    # Import our data into GenomicsDB:
    call SRJOINT.ImportGVCFs as ImportGVCFsIntoGenomicsDB {
        input:
            sample_name_map = CreateSampleNameMap.sample_name_map,
            interval_list = interval_list,
            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            ref_dict          = ref_map['dict'],
            prefix = prefix,
            batch_size = 50,
    }

    # Joint call
    call SRJOINT.GenotypeGVCFs as JointCallGVCFs {
        input:
            input_gvcf_data = ImportGVCFsIntoGenomicsDB.output_genomicsdb,
            interval_list = interval_list,
            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            ref_dict          = ref_map['dict'],
            dbsnp_vcf         = ref_map["known_sites_vcf"],
            prefix = prefix,
    }

    # Finalize
    call FF.FinalizeToFile as FinalizeGVCF { input: outdir = outdir, file = JointCallGVCFs.output_vcf }
    call FF.FinalizeToFile as FinalizeTBI { input: outdir = outdir, file = JointCallGVCFs.output_vcf_index }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File joint_gvcf = FinalizeGVCF.gcs_path
        File joint_gvcf_tbi = FinalizeTBI.gcs_path
    }
}


