version 1.0

import "tasks/Utils.wdl"
import "tasks/VariantUtils.wdl"

import "tasks/utils/qc/Fingerprinting.wdl" as FPUtils

workflow ECIdentity {

    meta {
        description: "Another one-off"
    }

    input {
        File aligned_bam
        File aligned_bai

        String all_vcf_tar_gz

        File gt_sites_union_bed

        File ref_map_file
    }

    parameter_meta {
        aligned_bam:        "GCS path to aligned BAM file of the flowcell"

        all_vcf_tar_gz:     "All FP vcf for suspect pool"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    call FPUtils.ExtractRelevantGenotypingReads {
        input:
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            genotyping_sites_bed = gt_sites_union_bed
    }

    call Hack {
        input:
            bam = ExtractRelevantGenotypingReads.relevant_reads,
            bai = ExtractRelevantGenotypingReads.relevant_reads_bai,
            haplotype_map = ref_map['haplotype_map'],
            all_vcf_tar_gz = all_vcf_tar_gz
    }

    # call FPUtils.ListGenotypedVCFs { input: fingerprint_store = fingerprint_store }
    # Array[String] all_vcfs = read_lines(ListGenotypedVCFs.vcf_gs_paths)
    # scatter (vcf in all_vcfs) {

    #     call VariantUtils.GetVCFSampleName {
    #         input:
    #             fingerprint_vcf = vcf
    #     }

    #     call FPUtils.CheckFingerprint {
    #         input:
    #             aligned_bam     = ExtractRelevantGenotypingReads.relevant_reads,
    #             aligned_bai     = ExtractRelevantGenotypingReads.relevant_reads_bai,
    #             fingerprint_vcf = vcf,
    #             vcf_sample_name = GetVCFSampleName.sample_name,
    #             haplotype_map   = ref_map['haplotype_map']
    #     }
    #     Float non_clr_lod = CheckFingerprint.metrics_map['LOD_EXPECTED_SAMPLE']
    # }

    # call FindMaxLOD { input: lod_scores = non_clr_lod }

    # if (FindMaxLOD.max_lod < 6) {
    if (Hack.best_lod < -3.0 ) {
        String bbid_na = "NotFound"
    }
    if (Hack.best_lod >= -3.0 ) {
        String bbid = Hack.best_bbid
    }
    # if (FindMaxLOD.max_lod >= 6) {
    #     String matching_vcf = all_vcfs[FindMaxLOD.idx]
    #     call GetIdentityInfo { input: vcf = matching_vcf }
    #     String bbid = GetIdentityInfo.resolved_bbid
    # }

    output {
        # Float highest_lod_from_exhaustive = FindMaxLOD.max_lod
        Float highest_lod_from_exhaustive = Hack.best_lod

        String true_bbid = select_first([bbid_na, bbid])
    }
}

# task FindMaxLOD {
#     input {
#         Array[String] lod_scores
#     }

#     Int n = length(lod_scores) - 1

#     command <<<

#         set -eux

#         seq 0 ~{n} > indices.txt
#         echo "~{sep='\n' lod_scores}" > scores.txt
#         paste -d' ' indices.txt scores.txt > to.sort.txt

#         sort -k2 -n to.sort.txt > sorted.txt
#         cat sorted.txt

#         tail -n 1 sorted.txt | awk '{print $1}' > "idx.txt"
#         tail -n 1 sorted.txt | awk '{print $2}' > "max_lod.txt"
#     >>>

#     output {
#         Int idx = read_int("idx.txt")
#         Float max_lod = read_float("max_lod.txt")
#     }

#     ###################
#     runtime {
#         cpu: 2
#         memory:  "4 GiB"
#         disks: "local-disk 50 HDD"
#         bootDiskSizeGb: 10
#         preemptible_tries:     3
#         max_retries:           2
#         docker:"gcr.io/cloud-marketplace/google/ubuntu2004:latest"
#     }
# }

# task GetIdentityInfo {
#     meta {
#         description: "Get collaborator participant, sample id, and SMID from a fingerprint VCF that follows the naming convension of smid__collabSmId_collabPartId.vcf"
#     }
#     input {
#         String vcf
#     }

#     String vcf_name = basename(basename(vcf, ".gz"), ".vcf")

#     command <<<
#         set -eux
#         echo ~{vcf_name} | grep -Eo "A[0-9]{9}" > "bbid.txt"
#     >>>
#     output {
#         String resolved_bbid = read_string("bbid.txt")
#     }

#     ###################
#     runtime {
#         cpu: 2
#         memory:  "4 GiB"
#         disks: "local-disk 50 HDD"
#         bootDiskSizeGb: 10
#         docker:"gcr.io/cloud-marketplace/google/ubuntu2004:latest"
#     }
# }

task Hack {
    input {
        File bam
        File bai
        File haplotype_map
        File all_vcf_tar_gz
    }

    String bam_base = basename(bam, ".bam")

    command <<<
        set -eux

        # unpack all
        mkdir -p all_vcfs
        tar -xf ~{all_vcf_tar_gz} -C /cromwell_root/all_vcfs/

        # one-by-one
        mkdir -p results
        touch bbid_lod.tsv
        for vcf in all_vcfs/*.vcf; do
            postfix=$(echo "${vcf}" | grep -Eo "A[0-9]+(\.[0-9]+)?")
            prefix="~{bam_base}__${postfix}"

            bash /opt/1vall.sh ~{bam} ${vcf} ~{haplotype_map} ${prefix} &> "${prefix}.log"

            lod=$(grep -F 'LOD_EXPECTED_SAMPLE' "${prefix}.metrics_map.txt" | awk '{print $2}')

            sm=$(echo "${postfix}" | grep -Eo "A[0-9]+")
            echo -e "${sm}\t${lod}" >> bbid_lod.tsv
            mv *.txt *.log results/
        done

        ls results | wc -l
        wc -l bbid_lod.tsv

        tar -zcf results.tar.gz -C results/ .

        # find max
        sort -r -n -k2,2 bbid_lod.tsv | head -n 1 > best.txt
        awk '{print $1}' best.txt > best_bbid.txt
        awk '{print $2}' best.txt > best_lod.txt
    >>>

    output {
        File detailed_results = "results.tar.gz"
        File result = "best.txt"
        String best_bbid = read_string("best_bbid.txt")
        Float best_lod = read_float("best_lod.txt")
    }

    runtime {
        cpu: 8
        memory: "32 GiB"
        disks: "local-disk 375 LOCAL"
        docker: "us.gcr.io/broad-dsp-lrma/lr-gatk:throwaway"
    }
}
