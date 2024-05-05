version 1.0

import "../QC/Fingerprinting.wdl" as FPUtils
import "../Utility/VariantUtils.wdl"

workflow FPCheckAoU {

    meta {
        description:
        "Check correctness of metadata on a (demultiplexed) alignmed BAM, by genotyping it's BAM generated with its metadata, against a fingerprint VCF. Practically assumes human GRCh38 reference."
    }
    parameter_meta {
        aligned_bam:        "GCS path to aligned BAM file, supposed to be of the same sample as from the fingerprinting (FP) VCF"
        tech:               "The technology used to generate this BAM. Currently, the following values are accepted: [ONT, Sequel, Revio]."
        force:              "If true, will force run the fingerprinting program, and the workflow may fail for various reasons; otherwise, if the BAM is too small, it will automatically fail this QC check."

        fp_vcf_store:       "Name of the bucket and prefix holding the fingerprint VCFs."
        fp_sample_id:       "UUID of the sample at the fingerprint store, used to fetch the fingerprinting VCF"

        ref_specific_haplotype_map: "Happlotype map file for the reference build used. See https://bit.ly/3QyZbwt "

        lod_expected_sample: "An LOD score assuming the BAM is the same sample as the FP VCF, i.e. BAM sourced from the 'expected sample'."

        lod_pass_threshold: "A numeric threshold for LOD above which the sample will be considered passing the FP check."
        lod_fail_threshold: "A numeric threshold for LOD below which the sample will be considered failing the FP check."

        FP_status:           "A single word summary on the result of FP check; one of [PASS, FAIL, BORDERLINE]."
        fingerprint_summary: "A file holding the summaries of LOD (a bit more detail than pass/fail)."
        fingerprint_details: "A file holding the detailed LOD at each FP site."
        chain_file: "Chain file for the GATK LiftoverVcf process, mapping between source and target genome builds."
        target_reference_sequence_fasta_file: "Reference sequence fasta file for the target genome build used in the liftover."
        target_reference_sequence_fasta_file_index: "Index file for the target reference sequence fasta, required for GATK processes."
        target_reference_sequence_fasta_file_dict: "Sequence dictionary for the target reference fasta, required by GATK for reference sequence validation."

        mem: "Memory allocated for the GATK LiftoverVcf process, specified in GB."
        preemptible_attempts: "Number of preemptible attempts for the task. Preemptible instances are cheaper but may be terminated unexpectedly."
        disk_space_gb: "Amount of disk space allocated for the task, specified in GB."
        cpu: "Number of CPU cores allocated for the task."
        boot_disk_size_gb: "Size of the boot disk for the virtual machine running the task, specified in GB."
    }

    input {
        File aligned_bam
        File aligned_bai
        String tech

        String fp_vcf_store
        String fp_sample_id

        File ref_specific_haplotype_map

        # Input args liftedoverVCFGATK generation:
        File chain_file
        File target_reference_sequence_fasta_file
        File target_reference_sequence_fasta_file_index
        File target_reference_sequence_fasta_file_dict
        # Runtime args:
        Int? mem
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb 

        Boolean force = false

        Float lod_pass_threshold =  6.0
        Float lod_fail_threshold = -3.0
    }
    # Generate output file names using fp_sample_id
    String lifted_over_vcf = fp_sample_id + ".lifted_over.vcf"
    String rejectedVcf = fp_sample_id + ".rejected.vcf"

    output {
        Float lod_expected_sample = fingerprint_check_LOD
        String FP_status = fingerprint_check_status

        File fingerprint_summary = select_first([CheckFingerprint.summary_metrics, "None"])
        File fingerprint_details = select_first([CheckFingerprint.detail_metrics, "None"])

    }

    # 1X coverage ~= 1.5GiB (Revio); 2.3GiB (Sequel); 2.8GiB (ONT)
    # below 0.5X, we define the bam as too small to draw a conclusion on, unless we're forced to run the program
    Map[String, Int] teeny_bam_def_sz = {"ONT":    1400,
                                         "Sequel": 1650,
                                         "Revio":   750}

    Int bam_sz_mb = ceil(size(aligned_bam, "MiB"))
    if (bam_sz_mb >= teeny_bam_def_sz[tech] || force) {
        ##### Prep work
        call ResolveFPVCFPath {input: fp_vcf_store = fp_vcf_store, fp_sample_id = fp_sample_id}
        call ReheaderFullGRCh38VCFtoNoAlt {input: full_GRCh38_vcf = ResolveFPVCFPath.fp_vcf}
        # liftover the vcfs
        call LiftoverVcfGATK {
            input:
                input_vcf_file                              = ReheaderFullGRCh38VCFtoNoAlt.reheadered_vcf,
                chain_file                                  = chain_file,
                target_reference_sequence_fasta_file        = target_reference_sequence_fasta_file,
                target_reference_sequence_fasta_file_index  = target_reference_sequence_fasta_file_index,
                target_reference_sequence_fasta_file_dict   = target_reference_sequence_fasta_file_dict,
                lifted_over_vcf_name                        = lifted_over_vcf,
                lifted_over_rejects_vcf_name                = rejectedVcf,
                mem                                         = mem,
                preemptible_attempts                        = preemptible_attempts,
                disk_space_gb                               = disk_space_gb,
                cpu                                         = cpu,
                boot_disk_size_gb                           = boot_disk_size_gb
        }
    
        call VariantUtils.GetVCFSampleName {
            input:
                fingerprint_vcf = LiftoverVcfGATK.lifted_over_vcf
        }
        call FPUtils.FilterGenotypesVCF {
            input:
                fingerprint_vcf = LiftoverVcfGATK.lifted_over_vcf
        }
        call FPUtils.ExtractGenotypingSites {
            input:
                fingerprint_vcf = FilterGenotypesVCF.ready_to_use_vcf
        }
        call FPUtils.ExtractRelevantGenotypingReads {
            input:
                aligned_bam     = aligned_bam,
                aligned_bai     = aligned_bai,
                genotyping_sites_bed = ExtractGenotypingSites.sites,
        }

        ##### check
        call FPUtils.CheckFingerprint {
            input:
                aligned_bam     = ExtractRelevantGenotypingReads.relevant_reads,
                aligned_bai     = ExtractRelevantGenotypingReads.relevant_reads_bai,
                fingerprint_vcf = FilterGenotypesVCF.ready_to_use_vcf,
                vcf_sample_name = GetVCFSampleName.sample_name,
                haplotype_map   = ref_specific_haplotype_map
        }

        ##### wrapup
        Float lod_expected_sample_t = CheckFingerprint.metrics_map['LOD_EXPECTED_SAMPLE']

        String status = if(lod_expected_sample_t < lod_fail_threshold) then "FAIL" else if (lod_expected_sample_t > lod_pass_threshold) then "PASS" else "BORDERLINE"
    }
    # FAIL the bam if coverage is below a certain threshold (not much useful data anyway)
    Float teeny_bam_lod = 0.0
    String fingerprint_check_status = select_first([status, "FAIL"])
    Float fingerprint_check_LOD = select_first([lod_expected_sample_t, teeny_bam_lod])
}

task ResolveFPVCFPath {
    meta {
        desciption:
        "Find the fingerprint VCF at the fingerprint store; project specific."
    }

    input {
        String fp_vcf_store
        String fp_sample_id
        RuntimeAttr? runtime_attr_override
    }

    String fp_vcf_store_formatted = sub(fp_vcf_store, "/$", "")

    command <<<
        set -eux

        # note the addition of the wildcard character *
        FP_SEARCH="~{fp_vcf_store_formatted}/~{fp_sample_id}*.fingerprint.liftedover.vcf"
        # this will error if no paths match, i.e. no FP file exists with this fp_sample_id
        FP_PATH=$(gsutil ls "${FP_SEARCH}" | head -n 1)
        FP_INDEX_PATH="${FP_PATH}.idx"

        echo "${FP_PATH}" > "vcf.gspath"
        echo "${FP_INDEX_PATH}" > "index.gspath"
    >>>

    output {
        String fp_vcf     = read_string("vcf.gspath")
        String fp_vcf_idx = read_string("index.gspath")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ReheaderFullGRCh38VCFtoNoAlt {
    meta {
        desciption:
        "Reheader the fingperint VCF that's generated with full GRCh38 reference to the no_alt header; project specific."
    }

    input {
        File full_GRCh38_vcf
    }

    command <<<
        set -eux

        grep -vF "_decoy,length=" ~{full_GRCh38_vcf} | \
            grep -vF "_alt,length=" | \
            grep -v "^##contig=<ID=HLA-" \
            > "reheadered.fp.vcf"
    >>>

    output {
        File reheadered_vcf = "reheadered.fp.vcf"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
task LiftoverVcfGATK {
    meta {
        description: "Use GATK's LiftoverVcf tool to lift over a VCF from one reference build to another."
    }
    input {
        # Input args:
        File input_vcf_file
        File chain_file
        File target_reference_sequence_fasta_file
        File target_reference_sequence_fasta_file_index
        File target_reference_sequence_fasta_file_dict
    
        # Output Names:
        String lifted_over_vcf_name
        String lifted_over_rejects_vcf_name

        # Runtime args:
        Int? mem
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb 
    
    }
    # Get machine settings:
    #Boolean use_ssd = true
    Int base_disk_space_gb = 100

    Int base_boot_disk_size_gb = 15

    # Timing output file
    #String timing_output_file = "liftover_timing.txt"

    command <<<
        set -euxo pipefail

        # startTime=`date +%s.%N`
        # echo "StartTime: $startTime" > ${timing_output_file}

        time \
        gatk LiftoverVcf \
            -I ~{input_vcf_file} \
            -O ~{lifted_over_vcf_name} \
            -CHAIN ~{chain_file} \
            -REJECT ~{lifted_over_rejects_vcf_name} \
            -R ~{target_reference_sequence_fasta_file} \
            -recover_swapped true
        # endTime=`date +%s.%N`
        # echo "EndTime: $endTime" >> ${timing_output_file}
        # elapsedTime=`python -c "print( $endTime - $startTime )"`
        # echo "Elapsed Time: $elapsedTime" >> ${timing_output_file}        
    >>>
    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.4.0.0"
        memory: select_first([mem, 8]) + " GB"
        cpu: select_first([cpu, 1])
        disks: "local-disk " + select_first([disk_space_gb, base_disk_space_gb]) + " SSD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, base_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 3])
    }
    
    # Outputs:
    output {
        File lifted_over_vcf         = "${lifted_over_vcf_name}"
        File lifted_over_rejects_vcf = "${lifted_over_rejects_vcf_name}"
        # File timing_info             = timing_output_file
    }
}

