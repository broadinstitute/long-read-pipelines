version 1.0

############################################################################################
## A workflow that preprocess a population-scale dataset for easy visualization.
############################################################################################

import "tasks/Utils.wdl"
import "tasks/VariantUtils.wdl"
import "tasks/Truvari.wdl"
import "tasks/SVTK.wdl"
import "tasks/Finalize.wdl" as FF

workflow LRVisualizeVariants {
    input {
        Array[File]+ aligned_bams
        Array[File]+ aligned_bais

        Array[File] callset1_vcfs
        Array[File] callset2_vcfs
        Array[File] callset3_vcfs

        File ref_map_file

        String prefix

        Int num_simultaneous_jobs = 100

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/LRVisualizeVariants/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    scatter (vcf in flatten([callset1_vcfs, callset2_vcfs, callset3_vcfs])) {
        call VariantUtils.VCFToBed { input: vcf = vcf }
    }

    call Utils.MergeBedFiles { input: beds = VCFToBed.bed, ref_fai = ref_map['fai'], dist = 50000, prefix = "union" }

    # Parallelize downloads over num_simultaneous_downloads jobs.
    scatter (n in range(num_simultaneous_jobs)) {
        call ExtractLoci {
            input:
                aligned_bams     = aligned_bams,
                aligned_bais     = aligned_bais,
                bed              = MergeBedFiles.merged_bed,
                nth              = n,
                num_jobs         = num_simultaneous_jobs,
                gcs_out_root_dir = outdir
        }
    }

    ##########
    # store the results into designated bucket
    ##########

    output {
        String viz_path = outdir
    }
}

task ExtractLoci {
    input {
        Array[String] aligned_bams
        Array[String] aligned_bais

        File bed
        Int nth
        Int num_jobs

        String gcs_out_root_dir

        RuntimeAttr? runtime_attr_override
    }

    # estimate required disk size in GB
    Int disk_size = 1

    command <<<
        # Authorize streaming access to GCS data
        source /opt/re-auth.sh
        set -euxo pipefail

        RET=0

        # Iterate over every nth line in the bed file
        awk 'NR % ~{num_jobs} == ~{nth} { print $0 }' ~{bed} | while read line; do
            # Extract locus pieces
            chr=$(echo $line | awk '{ print $1 }')
            start=$(echo $line | awk '{ print $2 }')
            stop=$(echo $line | awk '{ print $3 }')

            # Set up local and remote BAM file paths
            local_bam_file="viz.${chr}_${start}_${stop}.bam"
            local_bai_file="viz.${chr}_${start}_${stop}.bam.bai"
            gcs_bam_file="~{gcs_out_root_dir}/$local_bam_file"
            gcs_bai_file="~{gcs_out_root_dir}/$local_bai_file"

            if gsutil -q stat "$gcs_bam_file" ; then
                # Skip processing BAMs that we've already written to GCS
                echo "$gcs_bam_file already exists."
            else
                # Remove any lingering subset BAMs
                rm -f *.subset.bam*

                # Loop over all sample BAMs
                while read sample_bam; do
                    # Reauthorize streaming access to GCS with each iteration
                    source /opt/re-auth.sh

                    sample_bai=$(echo $sample_bam | sed 's/.bam/.bam.bai/')
                    name=$(basename $sample_bam .bam)

                    # Extract locus from sample BAM
                    samtools view \
                        -bhX \
                        -M \
                        -@ 1 \
                        --verbosity=8 \
                        --write-index \
                        -o "${name}.subset.bam##idx##${name}.subset.bai" \
                        $sample_bam $sample_bai \
                        ${chr}:${start}-${stop}
                done <~{write_lines(aligned_bams)}

                # Merge all pieces into a single BAM file
                samtools merge -o $local_bam_file *.subset.bam
                samtools index $local_bam_file

                # Upload BAM file to GCS
                gsutil -m cp $local_bam_file $gcs_bam_file
                gsutil -m cp $local_bai_file $gcs_bai_file
            fi
        done

        exit $RET
    >>>

    output {
        String out = read_string(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  10,
        max_retries:        0,
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
