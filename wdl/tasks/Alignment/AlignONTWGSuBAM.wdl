version 1.0

import "../Utility/Utils.wdl"
import "../Utility/BAMutils.wdl" as BU
import "AlignReads.wdl" as AR

workflow AlignONTWGSuBAM {

    meta {
        desciption:
        "Given an unaligned ONT BAM, align to a reference."
    }

    parameter_meta {
        bam_sample_name:       "Name of the sample to use in the aligned BAM."
        uBAM_tags_to_preserve: "Please be careful with this: methylation tags MM, ML are usually what you want to keep at minimum"
        flowcell:              "Name of the flowcell."

        aln_disk_type: "An optimization specifying which type of disk to use for the minimap2 alignment step."
    }

    input {
        File uBAM
        Array[String] uBAM_tags_to_preserve # please be careful with this: methylation tags MM, ML are usually what you want to keep
        String bam_sample_name
        String flowcell

        String aln_disk_type = 'SSD'

        File ref_map_file
    }

    output {
        File aligned_bam = aBAM
        File aligned_bai = aBAI
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String output_prefix = bam_sample_name + "." + flowcell

    ###################################################################################
    # massage the RG line to prepare for minimap2

    # for transferring the read group line
    call BU.GetReadGroupLines { input: bam = uBAM }
    if (1!=length(GetReadGroupLines.read_group_ids)) {
        call Utils.StopWorkflow { input: reason = "This workflow only works for single read-group inputs." }
    }

    # if the intended sample id isn't the same as that in the input, change the read group line
    call BU.GetReadGroupInfo { input: bam = uBAM, keys = ['SM']}
    if (bam_sample_name != GetReadGroupInfo.read_group_info['SM']) {
        call OneOffReplaceSampleID { input:
            orig_rg_line = GetReadGroupLines.read_group_lines[0],
            old_sample_name = GetReadGroupInfo.read_group_info['SM'],
            new_sample_name = bam_sample_name
        }
    }

    String fix_this_rg_line = select_first([OneOffReplaceSampleID.fixed_rg_line, GetReadGroupLines.read_group_lines[0]])
    call OneOffHandleSpacesInRGLine { input: orig_rg_line = fix_this_rg_line }
    call OneOffHandleTabInRGLine { input: orig_rg_line = OneOffHandleSpacesInRGLine.fixed_rg_line } # this MUST be the last fix
    String ready_to_use_rg_line = OneOffHandleTabInRGLine.fixed_rg_line

    ###################################################################################
    # (fundamental assumption: PacBio reads are handled by pbmm2 in a separate route)
    # shard-align when BAM is bigger than a certain threshold
    # here our enemy is GCP preemption:
    # emperically, 576K reads ~ 8GiB BAM for ONT reads
    # which takes roughly 2 hours for minimap2+sort task on SSD-type PD with 16 cores, and that rarely gets preemptied
    # we want to balance between wall-clock time and size of VM
    # both higher-count-CPU and longer wall-clock time lead to higher chance of preemption
    # also, since LOCAL SSD comes in 375G/pieces, we can compsensate the reads number with IO speed
    Pair[Int, Int] emperical_bamsz_and_readcnt = (8, 576000)
    Int ssd_pd_scale_factor = 3
    Int local_ssd_scale_factor = 4
    Int emperical_bam_sz_threshold = (if('LOCAL'==aln_disk_type) then local_ssd_scale_factor else ssd_pd_scale_factor) * emperical_bamsz_and_readcnt.left
    Int emperical_n_reads_per_shard = (if('LOCAL'==aln_disk_type) then local_ssd_scale_factor else ssd_pd_scale_factor) * emperical_bamsz_and_readcnt.right

    if (emperical_bam_sz_threshold < ceil(size(uBAM, "GiB"))) {
        call BU.SplitNameSortedUbam {
            input: uBAM = uBAM, n_reads = emperical_n_reads_per_shard
        }
        Int num_shards = length(SplitNameSortedUbam.split)
        scatter (i in range(num_shards)) {
            File shard_ubam = SplitNameSortedUbam.split[i]
            String shard_out_prefix = output_prefix + "." + "~{i}"
            call AR.Minimap2 as MapShard {
                input:
                    reads = [shard_ubam],
                    reads_file_basenames = [basename(shard_ubam)],
                    ref_fasta = ref_map['fasta'],
                    map_preset = 'map-ont',

                    RG = ready_to_use_rg_line,

                    tags_to_preserve = uBAM_tags_to_preserve,

                    prefix = shard_out_prefix,
                    disk_type = aln_disk_type
            }
        }
        call BU.MergeBamsWithSamtools as MergeBams { input: bams = MapShard.aligned_bam, out_prefix = output_prefix, disk_type = if ('LOCAL' == aln_disk_type) then 'SSD' else 'LOCAL'}
    }
    if (emperical_bam_sz_threshold >= ceil(size(uBAM, "GiB"))) {
        call AR.Minimap2 {
            input:
                reads = [uBAM],
                reads_file_basenames = [basename(uBAM)],
                ref_fasta = ref_map['fasta'],
                map_preset = 'map-ont',

                RG = ready_to_use_rg_line,

                tags_to_preserve = uBAM_tags_to_preserve,

                prefix = output_prefix,
                disk_type = aln_disk_type
        }
    }
    File aBAM = select_first([Minimap2.aligned_bam, MergeBams.merged_bam])
    File aBAI = select_first([Minimap2.aligned_bai, MergeBams.merged_bai])
}

###################################################################################
# one-off fix that is project specific
###################################################################################

task OneOffReplaceSampleID {
    meta {
        desciption: "Fixing sample ID in a BAM's read group line (because what's desired in the final BAM and what's in the input aren't the same, and reheader is expensive.)"
    }
    input {
        String orig_rg_line
        String old_sample_name
        String new_sample_name
    }

    output {
        String fixed_rg_line = read_lines("result.txt")[0]
    }
    command <<<
        set -eux

        echo -e "~{orig_rg_line}" | \
            sed "s/SM:~{old_sample_name}/SM:~{new_sample_name}/g" \
        > "result.txt"
    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task OneOffHandleTabInRGLine {
    meta {
        desciption: "RG input to minimap needs \\t, not just \t"
    }
    input {
        String orig_rg_line
    }

    output {
        String fixed_rg_line = read_lines("result.txt")[0]
    }
    command <<<
        set -eux

        echo -e "~{orig_rg_line}" | \
            sed 's#\t#\\t#g' \
        > "result.txt"
    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task OneOffHandleSpacesInRGLine {
    meta {
        desciption: "RG input sometimes contrain spaces (in DS field), re-format as ';' here. Otherwise this causes troubles in downstream WDLs because spaces may be used as indicating new args to commands."
    }
    input {
        String orig_rg_line
    }

    output {
        String fixed_rg_line = read_lines("result.txt")[0]
    }
    command <<<
        set -eux

        # replace any sequence of one or more spaces with ";"
        SPACE=" "
        echo -e "~{orig_rg_line}" | \
            sed "s/${SPACE}${SPACE}*/;/g" \
        > "result.txt"
    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}