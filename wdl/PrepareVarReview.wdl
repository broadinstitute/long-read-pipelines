version 1.0

import "tasks/Structs.wdl"
import "tasks/BWA.wdl" as BWA
import "tasks/AlignReads.wdl" as AlignReads
import "tasks/WGA.wdl" as WGA
import "tasks/Finalize.wdl" as FF

# This workflow aids in building a variant calling truth set, by comparing two (high quality) assemblies
# with tools like mummer and minimap2. This workflow then generates several plots and reports for each variant
# identified which can be used to manually validate each variant.

workflow PrepareVarReview {
    input {
        String sample1
        String sample2

        File ref1
        File ref2

        File? ref1_gff

        Array[File]+ s1_illumina_reads
        File? s1_nanopore_reads
        File? s1_pacbio_reads

        Array[File]+ s2_illumina_reads
        File? s2_nanopore_reads
        File? s2_pacbio_reads

        String output_dir
    }

    String gcs_output_dir = sub(output_dir, "/+$", "")

    parameter_meta {
        sample1: "Name of the first sample (e.g. strain name)"
        sample2: "Name of the second sample (e.g. strain name)"

        ref1: "FASTA file containing the assembled genome of sample 1"
        ref2: "FASTA file containing the assembled genome of sample 2"

        ref1_gff: "Gene annotations for the first reference, to be displayed in variant review reports"

        s1_illumina_reads: "Illumina read data from sample 1."
        s1_nanopore_reads: "Oxford nanopore reads from sample 1."
        s1_pacbio_reads: "Pacbio reads from sample 1."

        s2_illumina_reads: "Illumina read data from sample 2."
        s2_nanopore_reads: "Oxford nanopore reads from sample 2."
        s2_pacbio_reads: "Pacbio reads from sample 2."
    }

    # Read alignments to both references
    # ----------------------------------
    #
    
    call BWA.BWAMem2Align as S1_SelfAlignmentIllumina {
        input: ref=ref1, reads=s1_illumina_reads, output_prefix="illumina_self_~{sample1}"
    }

    call BWA.BWAMem2Align as S2_SelfAlignmentIllumina {
        input: ref=ref2, reads=s2_illumina_reads, output_prefix="illumina_self_~{sample2}"
    }

    call BWA.BWAMem2Align as S2_to_S1_Illumina {
        input: ref=ref1, reads=s2_illumina_reads, output_prefix="illumina_~{sample2}_to_~{sample1}"
    }

    call BWA.BWAMem2Align as S1_to_S2_Illumina {
        input: ref=ref2, reads=s1_illumina_reads, output_prefix="illumina_~{sample1}_to_~{sample2}"
    }

    if(defined(s1_nanopore_reads)) {
        File np_reads_s1 = select_first([s1_nanopore_reads])

        # Self alignment
        call AlignReads.Minimap2 as S1_SelfAlignmentNanopore {
            input: ref_fasta=ref1, reads=[np_reads_s1], map_preset="map-ont",
                prefix="ont_self_~{sample1}"
        }

        # Cross alignment
        call AlignReads.Minimap2 as S1_to_S2_Nanopore {
            input: ref_fasta=ref2, reads=[np_reads_s1], map_preset="map-ont",
                prefix="ont_~{sample1}_to_~{sample2}"
        }
    }

    if(defined(s1_pacbio_reads)) {
        File pb_reads_s1 = select_first([s1_pacbio_reads])

        # Self alignment
        call AlignReads.Minimap2 as S1_SelfAlignmentPacBio {
            input: ref_fasta=ref1, reads=[pb_reads_s1], map_preset="map-pb",
            prefix="pb_self_~{sample1}"
        }

        # Cross alignment
        call AlignReads.Minimap2 as S1_to_S2_PacBio {
            input: ref_fasta=ref2, reads=[pb_reads_s1], map_preset="map-pb",
                prefix="pb_~{sample1}_to_~{sample2}"
        }
    }

    if(defined(s2_nanopore_reads)) {
        File np_reads_s2 = select_first([s2_nanopore_reads])

        # Self-alignment
        call AlignReads.Minimap2 as S2_SelfAlignmentNanopore {
            input: ref_fasta=ref2, reads=[np_reads_s2], map_preset="map-ont",
                prefix="ont_self_~{sample2}"
        }

        # Cross alignment
        call AlignReads.Minimap2 as S2_to_S1_Nanopore {
            input: ref_fasta=ref1, reads=[np_reads_s2], map_preset="map-ont",
                prefix="ont_~{sample2}_to_~{sample1}"
        }
    }

    if(defined(s2_pacbio_reads)) {
        File pb_reads_s2 = select_first([s2_pacbio_reads])

        # Self-alignment
        call AlignReads.Minimap2 as S2_SelfAlignmentPacBio {
            input: ref_fasta=ref2, reads=[pb_reads_s2], map_preset="map-pb",
            prefix="pb_self_~{sample2}"
        }

        # Cross alignment
        call AlignReads.Minimap2 as S2_to_S1_PacBio {
            input: ref_fasta=ref1, reads=[pb_reads_s2], map_preset="map-pb",
                prefix="pb_~{sample2}_to_~{sample1}"
        }
    }

    # Whole genome alignments with Nucmer and Minimap2
    # ------------------------------------------------
    # 
    call WGA.Minimap2 as S2_to_S1_MinimapWGA {
        input: ref1=ref1, ref2=ref2, prefix="wga_~{sample2}_to_~{sample1}_minimap2"
    }

    call WGA.Nucmer as S2_to_S1_NucmerWGA {
        input: ref1=ref1, ref2=ref2, prefix="wga_~{sample2}_to_~{sample1}_nucmer"
    }

    call CallVariantsNucmer {
        input: reference=ref1, query=ref2, delta=S2_to_S1_NucmerWGA.filtered_delta, 
            prefix="~{sample2}_to_~{sample1}.nucmer"
    }

    call CallVariantsMinimap2 {
        input: reference=ref1, query=ref2, prefix="~{sample2}_to_~{sample1}.minimap2"
    }

    call MergeVcfs {
        input: vcf1=CallVariantsMinimap2.vcf, vcf2=CallVariantsNucmer.vcf,
            prefix="~{sample2}_to_~{sample1}.merged"
    }

    Array[File] bams = select_all([
        S2_to_S1_Illumina.aligned_bam,
        S2_to_S1_Nanopore.aligned_bam,
        S2_to_S1_PacBio.aligned_bam
    ])
    Array[String] names = select_all([
        "Illumina",
        if(defined(s2_nanopore_reads)) then "Nanopore" else "",
        if(defined(s2_pacbio_reads)) then "PacBio" else "", 
    ])

    scatter(bam in bams) {
        File bais = bam + ".bai"
    }

    call SVPlots {
        input: reference=ref1, query=ref2, delta=S2_to_S1_NucmerWGA.filtered_delta, vcf=MergeVcfs.bcf,
            bams=bams, bais=bais, names=names
    }

    call FF.FinalizeToDir as FinalizeIllumina {
        input: files=[
            S1_SelfAlignmentIllumina.aligned_bam,
            S1_SelfAlignmentIllumina.aligned_bai,
            S2_SelfAlignmentIllumina.aligned_bam,
            S2_SelfAlignmentIllumina.aligned_bai,

            S2_to_S1_Illumina.aligned_bam,
            S2_to_S1_Illumina.aligned_bai,
            S1_to_S2_Illumina.aligned_bam,
            S1_to_S2_Illumina.aligned_bai
        ], outdir=gcs_output_dir + "/illumina"
    }

    Array[File] s1_ont_files = select_all([
        S1_SelfAlignmentNanopore.aligned_bam,
        S1_SelfAlignmentNanopore.aligned_bai,
        S1_to_S2_Nanopore.aligned_bam,
        S1_to_S2_Nanopore.aligned_bai
    ])

    Array[File] s2_ont_files = select_all([
        S2_SelfAlignmentNanopore.aligned_bam,
        S2_SelfAlignmentNanopore.aligned_bai,
        S2_to_S1_Nanopore.aligned_bam,
        S2_to_S1_Nanopore.aligned_bai
    ])

    Array[File] s1_pb_files = select_all([
        S1_SelfAlignmentNanopore.aligned_bam,
        S1_SelfAlignmentNanopore.aligned_bai,
        S1_to_S2_Nanopore.aligned_bam,
        S1_to_S2_Nanopore.aligned_bai
    ])

    Array[File] s2_pb_files = select_all([
        S2_SelfAlignmentPacBio.aligned_bam,
        S2_SelfAlignmentPacBio.aligned_bai,
        S2_to_S1_PacBio.aligned_bam,
        S2_to_S1_PacBio.aligned_bai
    ])

    Array[File] ont_files = flatten([s1_ont_files, s2_ont_files])
    Array[File] pb_files = flatten([s2_pb_files, s2_pb_files])

    if(length(ont_files) > 0) {
        call FF.FinalizeToDir as FinalizeONT {
            input: files=ont_files, outdir=gcs_output_dir + "/oxford_nanopore"
        }
    }

    if(length(pb_files) > 0) {
        call FF.FinalizeToDir as FinalizePB {
            input: files=pb_files, outdir=gcs_output_dir + "/pacbio"
        }
    }

    call FF.FinalizeToDir as FinalizeWGA {
        input: files=[
            S2_to_S1_NucmerWGA.delta,
            S2_to_S1_NucmerWGA.filtered_delta,
            S2_to_S1_MinimapWGA.wga_bam,
            S2_to_S1_MinimapWGA.wga_bai
        ], outdir=gcs_output_dir + "/wga"
    }
            
    call FF.FinalizeToDir as FinalizeVariants {
        input: files=[
            CallVariantsNucmer.vcf,
            CallVariantsNucmer.insertions_fasta,
            CallVariantsMinimap2.vcf,
            MergeVcfs.bcf,
            MergeVcfs.tbi
        ], outdir=gcs_output_dir
    }

    call FF.FinalizeToDir as FinalizePlots {
        input: files=SVPlots.plots, outdir=gcs_output_dir + "/review/plots"
    }

}

task CallVariantsNucmer {
    input {
        File reference
        File query
        File delta

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    output {
        File vcf = "~{prefix}.vcf"
        File insertions_fasta = "~{prefix}.insertions.fasta"
    }

    command <<<
        set -euxo pipefail

        # Ensure correct file paths in delta
        patch_nucmer_delta -d "~{delta}" "~{reference}" "~{query}" > patched.delta

        nucmer_to_vcf "~{reference}" "~{query}" "patched.delta" \
            | bcftools sort \
            | bcftools norm -f "~{reference}" - > "~{prefix}.full.vcf"

        # Move long ALT sequences of large insertions to separate FASTA
        move_large_insertions -a "~{prefix}.insertions.fasta" "~{prefix}.full.vcf" > "~{prefix}.vcf"
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-variant-review:0.1.3"
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


task CallVariantsMinimap2 {
    input {
        File reference
        File query

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    output {
        File vcf = "~{prefix}.vcf"
    }
    
    command <<<
        set -euxo pipefail

        minimap2 -cx asm5 --cs "~{reference}" "~{query}" \
            | sort -k6,6 -k8,8n \
            | paftools.js call -f "~{reference}" - > "~{prefix}.vcf"
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-variant-review:0.1.3"
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

task MergeVcfs {
    input {
        File vcf1
        File vcf2

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    output {
        File bcf = "~{prefix}.vcf.gz"
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    command <<<
        set -euxo pipefail

        bgzip "~{vcf1}"
        bgzip "~{vcf2}"

        tabix -p vcf "~{vcf1}.gz"
        tabix -p vcf "~{vcf2}.gz"

        # Somehow -O b doesn't work with tabix?
        bcftools merge "~{vcf1}.gz" "~{vcf2}.gz" -O v > "~{prefix}.vcf"
        bgzip "~{prefix}.vcf"
        tabix "~{prefix}.vcf.gz"
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-variant-review:0.1.3"
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


task SVPlots {
    input {
        File reference
        File query
        File delta
        File vcf

        Array[File] bams
        Array[File] bais
        Array[String] names

        RuntimeAttr? runtime_attr_override
    }


    output {
        Array[File] plots = glob("*.png")
    }

    command <<<
        patch_nucmer_delta -d "~{delta}" "~{reference}" "~{query}" > patched.delta

        plot_sv_loci -r "~{reference}" -q "~{query}" -d "patched.delta" \
            -V "~{vcf}" -b ~{sep=" " bams} -n ~{sep=" " names} -o .
    >>>
        
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-variant-review:0.1.3"
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

