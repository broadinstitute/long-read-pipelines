version 1.0

import "../../../tasks/VariantCalling/GLNexus.wdl" as GLNexus
import "../../../tasks/Utility/Hail.wdl" as Hail
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow LRJointCallGVCFs {

    meta {
        description: "A workflow that performs joint calling on gVCFs (usually from DeepVariant) using GLNexus. For R9/PEPPER-Margin-DeepVariant inputs it first drops FILTER=NoCall records, which carry a malformed FORMAT/PL vector that makes GLnexus abort (see DropNoCallRecords)."
    }
    parameter_meta {
        gvcfs:                "GCS paths to gVCF files"
        tbis:                 "GCS paths to gVCF tbi files"
        ref_map_file:         "table indicating reference sequence and auxillary file locations"
        prefix:               "prefix for output joint-called gVCF and tabix index"
        gcs_out_root_dir:     "GCS bucket to store the reads, variants, and metrics files"
        drop_nocall_records:  "if true (default), strip FILTER=NoCall records from each input gVCF before joint calling. Required for R9/PEPPER-Margin-DeepVariant (DeepVariant 1.3.0) gVCFs, whose NoCall <*> sites have a one-short PL vector that trips GLnexus' `PL vector length N, expected M` check. Lossless (those records are GT=./., AD=0,0,0; the sample falls back to its reference band)."
    }

    input {
        Array[File] gvcfs
        Array[File] tbis
        File ref_map_file

        String prefix

        String gcs_out_root_dir

        Boolean drop_nocall_records = true
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/JointCallGVCFs/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    # For R9/PEPPER-Margin-DeepVariant inputs, drop FILTER=NoCall records whose
    # malformed PL vector would abort GLnexus. See DropNoCallRecords.
    if (drop_nocall_records) {
        scatter (i in range(length(gvcfs))) {
            call DropNoCallRecords { input: gvcf = gvcfs[i], tbi = tbis[i], idx = i }
        }
    }
    Array[File] jc_gvcfs = select_first([DropNoCallRecords.clean_gvcf, gvcfs])
    Array[File] jc_tbis  = select_first([DropNoCallRecords.clean_tbi,  tbis])

    # Gather across multiple input gVCFs
    call GLNexus.JointCall {
        input:
            gvcfs = jc_gvcfs,
            tbis = jc_tbis,
            dict = ref_map['dict'],
            prefix = prefix
    }

    call Hail.ConvertToHailMT {
        input:
            gvcf = JointCall.joint_gvcf,
            tbi = JointCall.joint_gvcf_tbi,
            prefix = prefix,
            outdir = outdir
    }

    # Finalize
    call FF.FinalizeToFile as FinalizeGVCF { input: outdir = outdir, file = JointCall.joint_gvcf }
    call FF.FinalizeToFile as FinalizeTBI { input: outdir = outdir, file = JointCall.joint_gvcf_tbi }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File joint_gvcf = FinalizeGVCF.gcs_path
        File joint_gvcf_tbi = FinalizeTBI.gcs_path
        String joint_mt = ConvertToHailMT.gcs_path
    }
}

task DropNoCallRecords {

    meta {
        description: "Remove FILTER=NoCall records from a (PEPPER-Margin-)DeepVariant gVCF. At <*> multiallelic NoCall sites, DeepVariant 1.3.0 emits a FORMAT/PL vector one entry short of the allele count (e.g. length 5 where 3 alleles need 6), which makes GLnexus abort with 'genotyper: unexpected result when fetching record FORMAT field ... PL vector length N, expected M' once such a site is unified across the cohort. These records are uncalled (GT=./., AD=0,0,0), so dropping them is lossless: the sample falls back to its underlying reference band at those positions. Ref: GLnexus issue #287. Verified end-to-end against a reproduced failure."
    }

    parameter_meta {
        gvcf: "input gVCF to filter"
        tbi:  "index for the input gVCF"
        idx:  "scatter index; used only to build a short, unique output basename (see out_prefix)"
    }

    input {
        File gvcf
        File tbi
        Int idx

        Int cpu = 2
        Int mem_gb = 4
    }

    # GLnexus builds a per-input "dataset name" from the sharded file name
    # (ShardVCFByRanges: "<contig-index>.<this-basename>.locus_<range>") and rejects
    # any that exceeds 100 chars (GLnexus regex_id, src/types.cc). Emit a short,
    # unique basename ("s<scatter-index>") so the dataset name stays under the cap.
    # Per-sample identity is unaffected: GLnexus reads sample names from the gVCF
    # header (untouched here), not from the file/dataset name.
    String out_prefix = "s~{idx}"
    Int disk_size = 10 + 2*ceil(size(gvcf, "GB"))

    command <<<
        set -euxo pipefail

        # Exclude only records whose FILTER is exactly NoCall (empirically a perfect
        # 1:1 with the malformed-PL records; GT=./. alone would also drop well-formed
        # RefCall bands, so filter on FILTER, not GT).
        bcftools view -e 'FILTER="NoCall"' -O z -o "~{out_prefix}.g.vcf.gz" ~{gvcf}
        bcftools index --tbi --force "~{out_prefix}.g.vcf.gz"
    >>>

    output {
        File clean_gvcf = "~{out_prefix}.g.vcf.gz"
        File clean_tbi  = "~{out_prefix}.g.vcf.gz.tbi"
    }

    runtime {
        cpu:                    cpu
        memory:                 mem_gb + " GiB"
        disks: "local-disk " +  disk_size + " SSD"
        bootDiskSizeGb:         10
        preemptible:            1
        maxRetries:             0
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
