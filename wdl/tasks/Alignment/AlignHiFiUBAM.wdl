version 1.0

import "../Utility/Utils.wdl"
import "../Utility/PBUtils.wdl" as PB
import "../Utility/BAMutils.wdl" as BU

workflow AlignHiFiUBAM {
    meta {
        desciption:
        "Align an CCS/HiFi unaligned BAM from a single readgroup."
    }

    parameter_meta {
        application: "Application of the data; accepted values are: [CCS, HIFI, ISOSEQ, MASSEQ]"
        library: "If provided, the output aligned BAM will use this value for the LB field in its readgroup header line (@RG)"
        drop_per_base_N_pulse_tags: "If false, the per-base and pulse tags will not be stripped away from the input uBAM, if they are available (please think through the storage implications)"
    }

    input {
        File  uBAM
        File? uPBI
        String bam_sample_name
        String? library

        String application

        Boolean drop_per_base_N_pulse_tags = true
        File ref_map_file

        String disk_type
    }

    output {
        File aligned_bam = aBAM
        File aligned_bai = aBAI
        File aligned_pbi = IndexAlignedReads.pbi

        String movie = movie_name
        String total_runtime = select_first([sumRuntimes.total_runtime, AlignReadsTogether.wallclocktime])
    }

    # todo: verify if this is still necessary
    Map[String, String] map_presets = {
        # 'CLR':    'SUBREAD',  # don't support any more
        'CCS':    'CCS',
        'HIFI':   'HIFI',
        'ISOSEQ': 'ISOSEQ',
        'MASSEQ': 'SUBREAD',
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    call BU.GetReadGroupInfo as RG {input: bam = uBAM, keys = ['SM', 'LB', 'PU']}
    String LB = select_first([library, RG.read_group_info['LB']])
    String movie_name = RG.read_group_info['PU']

    ###################################################################################
    Int shard_threshold = 50

    if (ceil(size(uBAM, "GiB")) > shard_threshold) {# shard & align

        if (! defined(uPBI) ) {
            call PB.PBIndex as PBIndex {input: bam = uBAM}
        }

        Int num_shards = ceil(size(uBAM, "GiB") / shard_threshold )
        call Utils.ComputeAllowedLocalSSD as Guess {input: intended_gb = 3*ceil(size(uBAM, "GiB") + size(uPBI, "GiB"))}
        call PB.ShardLongReads { input:
            unaligned_bam = uBAM, unaligned_pbi = select_first([uPBI, PBIndex.pbi]),
            num_shards = num_shards,
            num_ssds = Guess.numb_of_local_ssd
        }

        scatter (unaligned_bam in ShardLongReads.unmapped_shards) {
            # sometimes we see the sharded bams mising EOF marker, which then fails record counts, use this as a checkpoint
            call Utils.CountBamRecords as ValidateShard {input: bam = unaligned_bam}

            call PB.Align as AlignReads { input:
                bam         = unaligned_bam,
                sample_name = bam_sample_name,
                library     = LB,
                ref_fasta   = ref_map['fasta'],
                map_preset  = map_presets[application],
                drop_per_base_N_pulse_tags = drop_per_base_N_pulse_tags,
                disk_type = disk_type
            }
        }

        call Utils.MergeBams as MergeAlignedReads { input: bams = AlignReads.aligned_bam, prefix = basename(uBAM, ".bam") }
        call sumRuntimes { input: runtimes = AlignReads.wallclocktime }
    }
    if (! (ceil(size(uBAM, "GiB")) > shard_threshold)) {
        call PB.Align as AlignReadsTogether { input:
            bam         = uBAM,
            sample_name = bam_sample_name,
            library     = LB,
            ref_fasta   = ref_map['fasta'],
            map_preset  = map_presets[application],
            drop_per_base_N_pulse_tags = drop_per_base_N_pulse_tags,
            disk_type = disk_type
        }
    }

    File aBAM = select_first([MergeAlignedReads.merged_bam, AlignReadsTogether.aligned_bam])
    File aBAI = select_first([MergeAlignedReads.merged_bai, AlignReadsTogether.aligned_bai])
    call PB.PBIndex as IndexAlignedReads { input: bam = aBAM }
}

task sumRuntimes {
    input {
        Array[String] runtimes
    }

    output {
        String total_runtime = read_string("total.txt")
    }
    command <<<
        total_minutes=0

        # Read each runtime and sum up
        while IFS= read -r line; do
            # Extract hours and minutes from format "XX hours, YY minutes"
            hours=$(echo "$line" | grep -oP '\d+(?=H)')
            minutes=$(echo "$line" | grep -oP '\d+(?=M)')

            # Convert to total minutes and add
            total_minutes=$((total_minutes + 10#$hours * 60 + 10#$minutes))
        done < ~{write_lines(runtimes)}

        # Convert back to hours and minutes
        final_hours=$((total_minutes / 60))
        final_minutes=$((total_minutes % 60))

        # Output in same format with padding
        printf "%02dH%02dM\n" $final_hours $final_minutes > total.txt
    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
