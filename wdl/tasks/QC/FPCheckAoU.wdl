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

        fp_store:           "Name of the bucket and prefix holding the fingerprint VCFs."
        sample_id_at_store: "UUID of the sample at the fingerprint store, used to fetch the fingerprinting VCF"

        ref_specific_haplotype_map: "Happlotype map file for the reference build used. See https://bit.ly/3QyZbwt "

        lod_expected_sample: "An LOD score assuming the BAM is the same sample as the FP VCF, i.e. BAM sourced from the 'expected sample'."

        lod_pass_threshold: "A numeric threshold for LOD above which the sample will be considered passing the FP check."
        lod_fail_threshold: "A numeric threshold for LOD below which the sample will be considered failing the FP check."

        FP_status:           "A single word summary on the result of FP check; one of [PASS, FAIL, BORDERLINE]."
        fingerprint_summary: "A file holding the summaries of LOD (a bit more detail than pass/fail)."
        fingerprint_details: "A file holding the detailed LOD at each FP site."
    }

    input {
        File aligned_bam
        File aligned_bai

        String fp_store
        String sample_id_at_store

        File ref_specific_haplotype_map

        Float lod_pass_threshold =  6.0
        Float lod_fail_threshold = -3.0
    }

    ##### Prep work
    call ResolveFPVCFPath {input: fp_store = fp_store, sample_id_at_store = sample_id_at_store}
    call ReheaderFullGRCh38VCFtoNoAlt {input: full_GRCh38_vcf = ResolveFPVCFPath.fp_vcf}

    call VariantUtils.GetVCFSampleName {
        input:
            fingerprint_vcf = ReheaderFullGRCh38VCFtoNoAlt.reheadered_vcf
    }
    call FPUtils.FilterGenotypesVCF {
        input:
            fingerprint_vcf = ReheaderFullGRCh38VCFtoNoAlt.reheadered_vcf
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

    output {
        Float lod_expected_sample = lod_expected_sample_t
        String FP_status = status

        File fingerprint_summary = CheckFingerprint.summary_metrics
        File fingerprint_details = CheckFingerprint.detail_metrics
    }
}

task ResolveFPVCFPath {
    meta {
        desciption:
        "Find the fingerprint VCF at the fingerprint store; project specific."
    }

    input {
        String fp_store
        String sample_id_at_store
        RuntimeAttr? runtime_attr_override
    }

    String fp_store_formatted = sub(fp_store, "/$", "")

    command <<<
        set -eux

        # note the addition of the wildcard character *
        FP_SEARCH="~{fp_store_formatted}/~{sample_id_at_store}*.fingerprint.liftedover.vcf"
        # this will error if no paths match, i.e. no FP file exists with this sample_id_at_store
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
        boot_disk_gb:       25,
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
