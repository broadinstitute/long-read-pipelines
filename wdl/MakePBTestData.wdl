version 1.0

##########################################################################################
## A workflow that performs CCS correction and variant calling on PacBio HiFi reads from a
## single flow cell. The workflow shards the subreads into clusters and performs CCS in
## parallel on each cluster.  Error-corrected reads are then variant-called.  A number of
## metrics and figures are produced along the way.
##########################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/CallSVs.wdl" as SV
import "tasks/Figures.wdl" as FIG
import "tasks/Finalize.wdl" as FF
import "tasks/CallSmallVariants.wdl" as SMV

workflow MakePBTestData {
    input {
        File subreads_bam
        File subreads_pbi

        File ref_map_file
        String participant_name

        Array[String] loci

        String? gcs_out_root_dir
    }

    parameter_meta {
        subreads_bam:     "unaligned subreads BAM"
        subreads_pbi:     "index for unaligned subreads BAM"

        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        participant_name: "name of the participant from whom these samples were obtained"

        loci:             "loci to select"

        gcs_out_root_dir: "[optional] GCS bucket to store data"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    call Utils.GetDefaultDir { input: workflow_name = "MakePBTestData" }
    String outdir = sub(select_first([gcs_out_root_dir, GetDefaultDir.path]), "/$", "") + "/" + participant_name

#    # break one raw BAM into fixed number of shards
#    call PB.ShardLongReads { input: unaligned_bam = subreads_bam, unaligned_pbi = subreads_pbi, num_shards = 300 }
#
#    # scatter across subreads
#    scatter (subreads in ShardLongReads.unmapped_shards) {
#        call PB.Align as AlignUncorrected {
#            input:
#                bam         = subreads,
#                ref_fasta   = ref_map['fasta'],
#                sample_name = participant_name,
#                map_preset  = "SUBREAD",
#                runtime_attr_override = { 'mem_gb': 64 }
#        }
#
#        scatter (locus in loci) {
#            call Utils.SubsetBam {
#                input:
#                    bam   = AlignUncorrected.aligned_bam,
#                    bai   = AlignUncorrected.aligned_bai,
#                    locus = locus
#            }
#
#            call SelectReadsByName {
#                input:
#                    subreads_bam = subreads,
#                    subset_bam = SubsetBam.subset_bam
#            }
#        }
#
#        call Utils.MergeBams as MergeSelectedSubsets { input: bams = SelectReadsByName.selected_bam, prefix = "~{participant_name}.selected" }
#    }
#
#    # merge the subsetted BAMs into a single BAM
#    call Utils.MergeBams as MergeAllSelectedSubsets { input: bams = MergeSelectedSubsets.merged_bam, prefix = "~{participant_name}.selected" }
#    call PB.PBIndex { input: bam = MergeAllSelectedSubsets.merged_bam }

    call ListPBInputFiles { input: bam = subreads_bam }

#    call FixXMLSampleName {
#        input:
#            xml = ListPBInputFiles.files['subreadset.xml'],
#            sample_name = participant_name
#    }

    ##########
    # store the results into designated bucket
    ##########

#    call FF.FinalizeToDir as FinalizeUnmodifiedInputs {
#        input:
#            files = [
#                ListPBInputFiles.files['adapters.fasta'],
#                ListPBInputFiles.files['baz2bam_1.log'],
#                ListPBInputFiles.files['controls.adapters.fasta'],
#                ListPBInputFiles.files['controls.fasta'],
#                ListPBInputFiles.files['sts.xml'],
#                ListPBInputFiles.files['transferdone'],
#            ],
#            outdir = outdir
#    }

#    call FF.FinalizeToDir as FinalizeSubreadsetXML {
#        input:
#            files = [ FixXMLSampleName.new_xml ],
#            outdir = outdir
#    }
#
#    call FF.FinalizeToFile as FinalizeSubreadsBam {
#        input:
#            file = MergeAllSelectedSubsets.merged_bam,
#            outfile = outdir + "/" + basename(subreads_bam)
#    }
#
#    call FF.FinalizeToFile as FinalizeSubreadsPbi {
#        input:
#            file = PBIndex.pbi,
#            outfile = outdir + "/" + basename(subreads_bam) + ".pbi"
#    }
}

# A utility to subset a BAM by name
task SelectReadsByName {
    input {
        File subreads_bam
        File subset_bam
        String prefix = "selected"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        subreads_bam: "original subreads"
        subset_bam:   "subsetted bam"
        prefix: "prefix for output bam and bai file names"
    }

    Int disk_size = 4*ceil(size([subreads_bam, subset_bam], "GB"))

    command <<<
        set -euxo pipefail

        samtools view ~{subset_bam} | awk '{ print $1 }' | sed 's/_/ /g' | sort -t'/' -k2 -k3 -n | sed 's/ /_/g' | uniq > readnames.txt

        bri index -v ~{subreads_bam}
        samtools view -H ~{subreads_bam} > header.txt
        ((cat header.txt) && (cat readnames.txt | bri get ~{subreads_bam} -)) | samtools view -b > ~{prefix}.bam
    >>>

    output {
        File selected_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bri:0.1.22"
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

task FixXMLSampleName {
    input {
        File xml
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        xml: "subreadset.xml file path"
        sample_name: "new sample name"
    }

    Int disk_size = 4*ceil(size(xml, "GB"))

    command <<<
        set -euxo pipefail

        cat ~{xml} | perl -pe 's/(Well|Bio)Sample Name=".*?"/$1Sample Name="~{sample_name}"/' > ~{basename(xml)}
    >>>

    output {
        File new_xml = "~{basename(xml)}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task ListPBInputFiles {
    input {
        String bam

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: "bam file in the directory to list"
    }

    String dir = sub(sub(bam, basename(bam), ""), "/+$", "")

    command <<<
        set -euxo pipefail

        gsutil ls "~{dir}" | awk '{ a=$1; sub(".*/m[0-9]+_[0-9]+_[0-9]+.", "", a); print a "\t" $1 }' > files_map.txt
    >>>

    output {
        Map[String, File] files = read_map("files_map.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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
