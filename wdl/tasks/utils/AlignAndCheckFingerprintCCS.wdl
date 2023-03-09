version 1.0

import "qc/FPCheckAoU.wdl" as FPCheck
import "../metrics/CollectPacBioAlignedMetrics.wdl" as AlnMetrics
import "../PBUtils.wdl" as PB
import "../Utils.wdl"
import "GeneralUtils.wdl"

workflow AlignAndCheckFingerprintCCS {
    meta {
        desciption:
        "Given an unaligned CCS/HiFi BAM for a sample, align and verify fingerprint."
    }

    input {
        File uBAM
        File uPBI
        String bam_sample_name
        String? library

        Boolean turn_off_fingperprint_check
        String fp_store
        String sample_id_at_store

        File ref_map_file
    }

    parameter_meta {
        turn_off_fingperprint_check: "Please turn of fingerprint check if the reference is not GRCh38."
        fp_store:           "Bucket name and prefix (gs://...) storing the fingerprinting VCFs"
        sample_id_at_store: "Name of the sample at the fingerprint store."

        # outputs
        alignment_metrics_tar_gz : "A tar.gz file holding the custom alignment metrics."

        fp_status : "A summary string on the fingerprint check result, value is one of [PASS, FAIL, BORDERLINE]."
        fp_lod_expected_sample : "An LOD score assuming the BAM is the same sample as the FP VCF, i.e. BAM sourced from the 'expected sample'."
        fingerprint_detail_tar_gz : "A tar.gz file holding the fingerprinting details."
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    Int shard_and_align_threshold = 50

    Map[String, String] map_presets = {
        'CLR':    'SUBREAD',
        'CCS':    'CCS',
        'ISOSEQ': 'ISOSEQ',
        'MASSEQ': 'SUBREAD',
    }
    ###################################################################################
    if (ceil(size(uBAM, "GB")) > shard_and_align_threshold) {# shard & align, but practically never true

        call Utils.ComputeAllowedLocalSSD as Guess {input: intended_gb = 3*ceil(size(uBAM, "GB") + size(uPBI, "GB"))}
        call Utils.RandomZoneSpewer as arbitrary {input: num_of_zones = 3}
        Int num_shards = 10

        call PB.ShardLongReads {
            input:
                unaligned_bam = uBAM, unaligned_pbi = uPBI,
                num_shards = num_shards, num_ssds = Guess.numb_of_local_ssd, zones = arbitrary.zones
        }

        scatter (unaligned_bam in ShardLongReads.unmapped_shards) {
            # sometimes we see the sharded bams mising EOF marker, which then fails record counts, use this as a checkpoint
            call Utils.CountBamRecords as ValidateShard {input: bam = unaligned_bam}

            call PB.Align as AlignReads {
                input:
                    bam         = unaligned_bam,
                    ref_fasta   = ref_map['fasta'],
                    sample_name = bam_sample_name,
                    library     = library,
                    map_preset  = map_presets['CCS'],
                    drop_per_base_N_pulse_tags = true
            }
            # call Utils.BamToFastq as ShardedFQ { input: bam = unaligned_bam, prefix = basename(unaligned_bam, ".bam") }
        }

        call Utils.MergeBams as MergeAlignedReads { input: bams = AlignReads.aligned_bam, prefix = basename(uBAM, ".bam") }
        # call Utils.MergeFastqs as MergeAllFastqs { input: fastqs = ShardedFQ.reads_fq }
    }
    if (! (ceil(size(uBAM, "GB")) > shard_and_align_threshold)) {
        call PB.Align as AlignReadsTogether {
            input:
                bam         = uBAM,
                ref_fasta   = ref_map['fasta'],
                sample_name = bam_sample_name,
                library     = library,
                map_preset  = map_presets['CCS'],
                drop_per_base_N_pulse_tags = true
        }
        # call Utils.BamToFastq { input: bam = uBAM, prefix = basename(uBAM, ".bam") }
    }

    File aBAM = select_first([MergeAlignedReads.merged_bam, AlignReadsTogether.aligned_bam])
    File aBAI = select_first([MergeAlignedReads.merged_bai, AlignReadsTogether.aligned_bai])
    # call PB.PBIndex as IndexAlignedReads { input: bam = aBAM }
    # File FastQ = select_first([MergeAllFastqs.merged_fastq, BamToFastq.reads_fq])

    ###################################################################################
    # alignment metrics and fingerprint check
    call AlnMetrics.CollectPacBioAlignedMetrics {
        input:
            aligned_bam = aBAM,
            aligned_bai = aBAI,
            # aligned_pbi = IndexAlignedReads.pbi
    }

    call GeneralUtils.TarGZFiles as saveAlnMetrics {
        input:
            files = flatten([[CollectPacBioAlignedMetrics.custom_aln_metrics_summary, CollectPacBioAlignedMetrics.nanoplot_stats], CollectPacBioAlignedMetrics.nanoplot_pngs]),
            name = "alignment.metrics"
    }

    if (!turn_off_fingperprint_check){
        call FPCheck.FPCheckAoU {
            input:
                aligned_bam = aBAM,
                aligned_bai = aBAI,
                fp_store = fp_store,
                sample_id_at_store = sample_id_at_store,
                ref_specific_haplotype_map = ref_map['haplotype_map']
        }
        call GeneralUtils.TarGZFiles as saveFPRes {input: files = [FPCheckAoU.fingerprint_summary, FPCheckAoU.fingerprint_details], name = 'fingerprint_check.summary_and_details'}
    }

    output {
        File aligned_bam = aBAM
        File aligned_bai = aBAI
        # File aligned_pbi = IndexAlignedReads.pbi

        File alignment_metrics_tar_gz = saveAlnMetrics.you_got_it

        # File fastq = FastQ

        Float? fp_lod_expected_sample = FPCheckAoU.lod_expected_sample
        String? fp_status = FPCheckAoU.FP_status
        File? fingerprint_detail_tar_gz = saveFPRes.you_got_it
    }
}
