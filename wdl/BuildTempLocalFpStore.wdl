version 1.0

import "tasks/Finalize.wdl" as Transfer

import "tasks/utils/qc/Fingerprinting.wdl" as FPUtils
import "tasks/utils/qc/FPCheckAoU.wdl"

workflow BuildTempLocalFpStore {
    meta {
        desciption: "1-off. Throw-away ware."
    }
    input {
        String fingerprint_store
        String local_path

        String name_for_all_fp_sites_bed
        String name_for_samples_without_fp

        File sample_ids_pool
    }

    call FilterDown { input: fingerprint_store = fingerprint_store, sample_ids_pool = sample_ids_pool }
    scatter (vcf in read_lines(FilterDown.available_paths)) {

        call FPCheckAoU.ReheaderFullGRCh38VCFtoNoAlt as reheader { input: full_GRCh38_vcf = vcf }

        call FPUtils.FilterGenotypesVCF { input: fingerprint_vcf = vcf }
        call FPUtils.ExtractGenotypingSites { input: fingerprint_vcf = FilterGenotypesVCF.ready_to_use_vcf }

        String out_name = basename(vcf, ".vcf") + "grch38_noalt.site_filtered.vcf"
        call Transfer.FinalizeToFile {input: file = reheader.reheadered_vcf, outdir = local_path, name = out_name}
    }

    call FPUtils.MergeGenotypingSites {input: all_sites = ExtractGenotypingSites.sites}
    call Transfer.FinalizeToFile as ffbed {input: file = MergeGenotypingSites.merged_sites, outdir = local_path, name = name_for_all_fp_sites_bed}

    call Transfer.FinalizeToFile as ffmask {input: file = FilterDown.unavailable_samples, outdir = local_path, name = name_for_samples_without_fp}

    output {
        File all_current_genotyping_sites = ffbed.gcs_path
        File samples_without_fp = ffmask.gcs_path
    }
}

task FilterDown {
    meta {
        volatile: true
    }
    input {
        String fingerprint_store
        File sample_ids_pool
    }

    String fp_store_formatted = sub(fingerprint_store, "/$", "")
    Array[String] sample_ids = read_lines(sample_ids_pool)

    command <<<
        set -eux

        touch fpvcf.txt
        touch no_fp_samples.txt

        gsutil ls ~{fingerprint_store} > everything.txt

        for sm in ~{sep=' ' sample_ids}; do
            if grep -qF "~{fp_store_formatted}/${sm}" everything.txt; then
                grep -F "~{fp_store_formatted}/${sm}" everything.txt >> fpvcf.txt
            else
                echo "${sm}" >> no_fp_samples.txt
            fi
        done
    >>>

    output {
        File available_paths = "fpvcf.txt"
        File unavailable_samples = "no_fp_samples.txt"
    }

    ###################
    runtime {
        cpu: 2
        memory:  "4 GiB"
        disks: "local-disk 200 HDD"
        docker:"us.gcr.io/broad-dsp-lrma/lr-basic:latest"
    }
}
