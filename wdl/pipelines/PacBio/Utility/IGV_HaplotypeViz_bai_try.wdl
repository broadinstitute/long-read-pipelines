version 1.0

workflow IGVScreenshotWorkflow {

    input {
        File aligned_bam_hap1
        File aligned_bam_hap1_bai
        File aligned_bam_hap2
        File aligned_bam_hap2_bai
        File alignments
        File alignments_bai

        Array[File]? orthogonal_alignments
        Array[File]? orthogonal_alignments_bai
        File? vcf
        File? vcf_tbi

        File bed_file
        File fasta_file
        File fasta_file_fai   # Include the .fai file
        String sample_name
        Int image_height = 500
        Int memory_mb = 4000
        Int disk_gb = 100           # Disk size in GB, default to 100 GB
        String docker_image_tag = "v982024"  # The Docker image to use
    }

    call RunIGVScreenshot {
        input:
            aligned_bam_hap1 = aligned_bam_hap1,
            aligned_bam_hap1_bai = aligned_bam_hap1_bai,
            aligned_bam_hap2 = aligned_bam_hap2,
            aligned_bam_hap2_bai = aligned_bam_hap2_bai,
            alignments = alignments,
            alignments_bai = alignments_bai,
            orthogonal_alignments = orthogonal_alignments,
            orthogonal_alignments_bai = orthogonal_alignments_bai,
            vcf = vcf,
            vcf_tbi = vcf_tbi,
            bed_file = bed_file,
            fasta_file = fasta_file,
            fasta_file_fai = fasta_file_fai,
            sample_name = sample_name,
            image_height = image_height,
            memory_mb = memory_mb,
            disk_gb = disk_gb,
            docker_image_tag = docker_image_tag
    }

    output {
        Array[File] snapshots = RunIGVScreenshot.snapshots
    }
}

task RunIGVScreenshot {

    input {
        File aligned_bam_hap1
        File aligned_bam_hap1_bai
        File aligned_bam_hap2
        File aligned_bam_hap2_bai
        File alignments
        File alignments_bai

        Array[File]? orthogonal_alignments
        Array[File]? orthogonal_alignments_bai
        File? vcf
        File? vcf_tbi

        File bed_file
        File fasta_file
        File fasta_file_fai
        String sample_name
        Int image_height
        Int memory_mb
        Int disk_gb
        String docker_image_tag
    }

    Boolean fail_early_a = (defined(vcf) == defined(vcf_tbi))
    Boolean fail_early_b = (defined(orthogonal_alignments) == defined(orthogonal_alignments_bai))
    Boolean fail_early_c = if defined(orthogonal_alignments) then (length(select_first([orthogonal_alignments])) == length(select_first([orthogonal_alignments_bai]))) else false
    Boolean fail_early = (fail_early_a && fail_early_b && fail_early_c)

    # this is when you get to hate WDL
    Array[File] all_alignments = if defined(orthogonal_alignments) then flatten([[aligned_bam_hap1, aligned_bam_hap2, alignments],
                                                                                 select_first([orthogonal_alignments])])
                                                                   else [aligned_bam_hap1, aligned_bam_hap2, alignments]

    Array[File] all_tracks = if defined(vcf) then flatten([select_all([vcf]), all_alignments])
                                             else all_alignments

    Int scaled_image_height = image_height + 100*(length(all_tracks) - 3)

    command <<<
        set -euo pipefail

        if ~{fail_early} ; then
            echo "ERROR: VCF and VCF index must be provided together. Same for orthogonal alignments and their indexes, plus they must of the same length."
            exit 1
        fi

        # Ensure the snapshots directory exists under the mounted disk path
        mkdir -p 'output/IGV_Snapshots/'

        # Start a virtual frame buffer to allow IGV to render
        Xvfb :1 -screen 0 1024x768x16 &> xvfb.log &
        export DISPLAY=:1

        # Run the IGV screenshot script with the provided inputs, no --snapshot-dir
        python3 /opt/IGV_Linux_2.18.2/make_igv_screenshot.py \
          ~{sep=' ' all_tracks} \
          -r ~{bed_file} \
          -ht ~{scaled_image_height} \
          -bin /opt/IGV_Linux_2.18.2/igv.sh \
          -mem ~{memory_mb} \
          --fasta_file ~{fasta_file} \
          --sample_name ~{sample_name}

        # Move the screenshots to the output directory
        # mv -- *.png 'output/IGV_Snapshots/'
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/igv_screenshot_docker:~{docker_image_tag}"
        memory: "~{memory_mb} MB"
        cpu: 2
        disks: "local-disk ~{disk_gb} SSD"
    }

    output {
        Array[File] snapshots = glob("output/IGV_Snapshots/*.png")
    }
}
