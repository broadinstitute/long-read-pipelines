version 1.0

import "../../structs/Structs.wdl"

workflow Hifiasm {

    meta {
        description: "We run two HiFiasm jobs, one for getting alternative contigs and one for getting the haplotigs. And we take the primary assembly from the first job."
    }
    parameter_meta {
        reads:    "reads (in fasta or fastq format, compressed or uncompressed)"
        prefix:   "prefix to apply to assembly output filenames"
    }

    input {
        File reads
        File? ont_reads_fastq
        String prefix

        Int kmer_size = 51
        Int bloom_filter_bits = 37
        Int minimizer_window_size = 51

        Boolean haploid = false

        String extra_args = ""

        String? telomere_5_prime_sequence

        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        RuntimeAttr? runtime_attr_override
    }

    call AssembleForAltContigs {
        input:
            reads  = reads,
            ont_reads_fastq = ont_reads_fastq,
            prefix = prefix,
            telomere_5_prime_sequence = telomere_5_prime_sequence,
            haploid = haploid,
            kmer_size = kmer_size,
            bloom_filter_bits = bloom_filter_bits,
            minimizer_window_size = minimizer_window_size,
            extra_args = extra_args,
            zones = zones,
            runtime_attr_override = runtime_attr_override
    }

    call AssembleForHaplotigs {
        input:
            reads  = reads,
            ont_reads_fastq = ont_reads_fastq,
            prefix = prefix,
            telomere_5_prime_sequence = telomere_5_prime_sequence,
            kmer_size = kmer_size,
            bloom_filter_bits = bloom_filter_bits,
            minimizer_window_size = minimizer_window_size,
            haploid = haploid,
            extra_args = extra_args,
            zones = zones,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File primary_gfa  = AssembleForAltContigs.primary_gfa
        File primary_tigs = AssembleForAltContigs.primary_tigs

        File alternate_gfa  = AssembleForAltContigs.alternate_gfa
        File alternate_tigs = AssembleForAltContigs.alternate_tigs

        File alt_contig_assembly_all_outputs_gz = AssembleForAltContigs.all_outputs_gz

        File log_in_pVSa_mode = AssembleForAltContigs.log

        ###########
        Array[File]? phased_gfas = AssembleForHaplotigs.phased_gfas
        Array[File]? phased_tigs = AssembleForHaplotigs.phased_tigs

        File haplotig_assembly_all_outputs_gz = AssembleForHaplotigs.all_outputs_gz

        File log_in_hap_mode = AssembleForHaplotigs.log

        # these two are saved, but the one generated in the primary VS alternate mode are preferred
        File primary_gfa_in_hap_mode  = AssembleForHaplotigs.primary_gfa
        File primary_tigs_in_hap_mode = AssembleForHaplotigs.primary_fa
    }
}

task AssembleForHaplotigs {
    input {
        File reads
        File? ont_reads_fastq
        String prefix = "out"
        String? telomere_5_prime_sequence
        
        Boolean haploid = false

        Int kmer_size = 51
        Int bloom_filter_bits = 37
        Int minimizer_window_size = 51

        String extra_args = ""

        String zones

        RuntimeAttr? runtime_attr_override
    }

    String telomere_5_prime_sequence_arg = if defined(telomere_5_prime_sequence) then "--telo-m " + select_first([telomere_5_prime_sequence]) else ""

    Int proposed_memory = 4 * ceil(size(reads, "GB"))
    Int memory = if proposed_memory < 96 then 96 else proposed_memory  # this 96 magic number is purely empirical
    Int n = memory / 4  # this might be an odd number
    Int num_cpus_proposal = if (n/2)*2 == n then n else n+1  # a hack because WDL doesn't have modulus operator
    Int num_cpus = if num_cpus_proposal > 96 then 96 else num_cpus_proposal

    Int disk_size = 50 + (10 * ceil(size(reads, "GB")) + 10 * ceil(size(ont_reads_fastq, "GB")))

    command <<<
        set -euxo pipefail

        echo "--------------------------------" | tee hifiasm.log
        echo "Proposed memory: ~{proposed_memory} GiB" | tee -a hifiasm.log
        echo "Memory requested: ~{memory} GiB" | tee -a hifiasm.log
        echo "Memory available: $(free -m | grep '^Mem' | awk '{print $2}') MB" | tee -a hifiasm.log
        echo "n: ~{n}" | tee -a hifiasm.log
        echo "num_cpus_proposal: ~{num_cpus_proposal}" | tee -a hifiasm.log
        echo "num_cpus: ~{num_cpus}" | tee -a hifiasm.log
        echo "default_attr: ~{default_attr}" | tee -a hifiasm.log
        echo "runtime_attr: ~{runtime_attr}" | tee -a hifiasm.log
        echo "--------------------------------" | tee -a hifiasm.log

        time hifiasm \
            -o ~{prefix} \
            -t~{num_cpus} \
            -k ~{kmer_size} \
            -f ~{bloom_filter_bits} \
            -m ~{minimizer_window_size} \
            ~{true="-l0" false="" haploid} \
            ~{true="--n-hap 1" false="" haploid} \
            ~{telomere_5_prime_sequence_arg} \
            ~{extra_args} \
            ~{reads} \
            ~{ont_reads_fastq} \
            2>&1 | tee -a hifiasm.log

        tree -h .

        # GFA graph to contigs, primary
        # outputs generated this way has "bp" in their names
        awk '/^S/{print ">"$2; print $3}' \
            ~{prefix}.bp.p_ctg.gfa \
            > ~{prefix}.bp.p_ctg.fa

        # GFA graph to contigs, for each haplotig set
        for gfa in *.gfa; do
            filename=$(basename -- "${gfa}")
            haplotype="${filename%.*}"
            awk '/^S/{print ">"$2; print $3}' \
                "${gfa}" \
                > "${haplotype}".fa
        done

        du -hs 

        # Save everything in a tar.gz file:
        tar -zcf ~{prefix}_hifiasm_assembleforhaplotigs.tar.gz ~{prefix}.*
    >>>

    output {
        # these are saved, but the one with alt contigs genearted will be preferred for now
        File primary_gfa = "~{prefix}.bp.p_ctg.gfa"
        File primary_fa = "~{prefix}.bp.p_ctg.fa"

        Array[File]? phased_gfas = glob("~{prefix}.bp.hap*.p_ctg.gfa")
        Array[File]? phased_tigs = glob("~{prefix}.bp.hap*.p_ctg.fa")

        File all_outputs_gz = "~{prefix}_hifiasm_assembleforhaplotigs.tar.gz"

        File log = "hifiasm.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.24.0"
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
        zones:                  zones
    }
}

task AssembleForAltContigs {
    input {
        File reads
        File? ont_reads_fastq
        String prefix = "out"
        String? telomere_5_prime_sequence

        Int kmer_size = 51
        Int bloom_filter_bits = 37
        Int minimizer_window_size = 51
        Boolean haploid = false

        String extra_args = ""

        String zones

        RuntimeAttr? runtime_attr_override
    }

    String telomere_5_prime_sequence_arg = if defined(telomere_5_prime_sequence) then "--telo-m " + select_first([telomere_5_prime_sequence]) else ""

    Int proposed_memory = 4 * ceil(size(reads, "GB"))
    Int memory = if proposed_memory < 96 then 96 else if proposed_memory > 512 then 512 else proposed_memory # this 96 magic number is purely empirical
    Int n = memory / 4  # this might be an odd number
    Int num_cpus_proposal = if (n/2)*2 == n then n else n+1  # a hack because WDL doesn't have modulus operator
    Int num_cpus = if num_cpus_proposal > 96 then 96 else num_cpus_proposal

    Int disk_size = 50 + (10 * ceil(size(reads, "GB")) + 10 * ceil(size(ont_reads_fastq, "GB")))

    command <<<
        set -euxo pipefail

        echo "--------------------------------" | tee hifiasm.log
        echo "Proposed memory: ~{proposed_memory} GiB" | tee -a hifiasm.log
        echo "Memory requested: ~{memory} GiB" | tee -a hifiasm.log
        echo "Memory available: $(free -m | grep '^Mem' | awk '{print $2}') MB" | tee -a hifiasm.log
        echo "n: ~{n}" | tee -a hifiasm.log
        echo "num_cpus_proposal: ~{num_cpus_proposal}" | tee -a hifiasm.log
        echo "num_cpus: ~{num_cpus}" | tee -a hifiasm.log
        echo "default_attr: ~{default_attr}" | tee -a hifiasm.log
        echo "runtime_attr: ~{runtime_attr}" | tee -a hifiasm.log
        echo "--------------------------------" | tee -a hifiasm.log

        time hifiasm \
            -o ~{prefix} \
            -t~{num_cpus} \
            --primary \
            -k ~{kmer_size} \
            -f ~{bloom_filter_bits} \
            -m ~{minimizer_window_size} \
            ~{true="-l0" false="" haploid} \
            ~{true="--n-hap 1" false="" haploid} \
            ~{telomere_5_prime_sequence_arg} \
            ~{extra_args} \
            ~{reads} \
            ~{ont_reads_fastq} \
            2>&1 | tee -a hifiasm.log

        tree -h .

        # tricky, outputs generated this way has no "bp" in their file names
        # GFA graph to contigs, primary
        awk '/^S/{print ">"$2; print $3}' \
            ~{prefix}.p_ctg.gfa \
            > ~{prefix}.p_ctg.fa

        # GFA graph to contigs, alternate
        awk '/^S/{print ">"$2; print $3}' \
            ~{prefix}.a_ctg.gfa \
            > ~{prefix}.a_ctg.fa

        du -hs 

        # Save everything in a tar.gz file:
        tar -zcf ~{prefix}_hifiasm_assembleforaltcontigs.tar.gz ~{prefix}.*
    >>>

    output {
        File primary_gfa  = "~{prefix}.p_ctg.gfa"
        File primary_tigs = "~{prefix}.p_ctg.fa"

        File alternate_gfa  = "~{prefix}.a_ctg.gfa"
        File alternate_tigs = "~{prefix}.a_ctg.fa"

        File all_outputs_gz = "~{prefix}_hifiasm_assembleforaltcontigs.tar.gz"

        File log = "hifiasm.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.24.0"
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
        zones:                  zones
    }
}
