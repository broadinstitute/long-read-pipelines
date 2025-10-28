version 1.0

import "../../structs/Structs.wdl"

workflow HifiasmTrio {

    meta {
        description: "We run two HiFiasm jobs, one for getting alternative contigs and one for getting the haplotigs. And we take the primary assembly from the first job."
    }
    parameter_meta {
        reads:    "reads (in fasta or fastq format, compressed or uncompressed)"
        ont_ultralong_reads: "ONT reads (ideally ultralong, in fastq format, compressed or uncompressed)"
        prefix:   "prefix to apply to assembly output filenames"
    }

    input {
        File reads
        File? ont_ultralong_reads

        File? maternal_fastq_1
        File? maternal_fastq_2
        File? paternal_fastq_1
        File? paternal_fastq_2

        String? telomere_5_prime_sequence

        String prefix

        Boolean haploid = false

        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }

    call AssembleTrioData as t001_AssembleTrioData {
        input:
            reads  = reads,
            ont_ultralong_reads = ont_ultralong_reads,
            maternal_fastq_1 = maternal_fastq_1,
            maternal_fastq_2 = maternal_fastq_2,
            paternal_fastq_1 = paternal_fastq_1,
            paternal_fastq_2 = paternal_fastq_2,
            telomere_5_prime_sequence = telomere_5_prime_sequence,
            haploid = haploid,
            prefix = prefix,
            zones = zones
    }

    output {
        File maternal_phased_gfa  = t001_AssembleTrioData.maternal_phased_gfa
        File maternal_phased_fa = t001_AssembleTrioData.maternal_phased_fa

        File paternal_phased_gfa  = t001_AssembleTrioData.paternal_phased_gfa
        File paternal_phased_fa = t001_AssembleTrioData.paternal_phased_fa

        File merged_unitigs_gfa = t001_AssembleTrioData.merged_unitigs_gfa
        File merged_unitigs_fa = t001_AssembleTrioData.merged_unitigs_fa

        File all_outputs_gz = t001_AssembleTrioData.all_outputs_gz

        File log = t001_AssembleTrioData.log
    }
}

task AssembleTrioData {
    input {
        File reads
        String prefix = "out"
        String zones

        File? ont_ultralong_reads

        File? maternal_fastq_1
        File? maternal_fastq_2
        File? paternal_fastq_1
        File? paternal_fastq_2

        String? telomere_5_prime_sequence

        Boolean haploid = false

        Int kmer_size = 51
        Int homopolymer_kmer_size = 37
        Int yak_bloom_filter_size = 37

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int proposed_memory = 4 * ceil(size(reads, "GB"))
    Int memory = if proposed_memory < 96 then 96 else proposed_memory  # this 96 magic number is purely empirical
    Int n = memory / 4  # this might be an odd number
    Int num_cpus_proposal = if (n/2)*2 == n then n else n+1  # a hack because WDL doesn't have modulus operator
    Int num_cpus = if num_cpus_proposal > 96 then 96 else num_cpus_proposal

    Int disk_size = 10 + 10 * ceil(size(reads, "GB")) + 10 * ceil(size(ont_ultralong_reads, "GB"))

    String ont_ultralong_reads_arg = if defined(ont_ultralong_reads) then "--ul " + select_first([ont_ultralong_reads]) else ""
    String telomere_5_prime_sequence_arg = if defined(telomere_5_prime_sequence) then "--telo-m " + select_first([telomere_5_prime_sequence]) else ""

    command <<<
        set -uxo pipefail

        # Check if we have the parental fastq files for trio assembly.
        # They must be either all defined or undefined:
        defined_count=0
        if [[ -n "~{maternal_fastq_1}" ]]; then ((defined_count++)); fi
        if [[ -n "~{maternal_fastq_2}" ]]; then ((defined_count++)); fi
        if [[ -n "~{paternal_fastq_1}" ]]; then ((defined_count++)); fi
        if [[ -n "~{paternal_fastq_2}" ]]; then ((defined_count++)); fi
        
        if [[ $defined_count -ne 4 ]] && [[ $defined_count -ne 2 ]] && [[ $defined_count -ne 0 ]]; then
            echo "Error: Either 2, 4, or none of the parental fastq files must be defined." 1>&2
            echo "       Maternal fastq 1: ~{maternal_fastq_1}" 1>&2
            echo "       Maternal fastq 2: ~{maternal_fastq_2}" 1>&2
            echo "       Paternal fastq 1: ~{paternal_fastq_1}" 1>&2
            echo "       Paternal fastq 2: ~{paternal_fastq_2}" 1>&2
            exit 1
        fi

        # Make sure we use all our processors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            let np=${np}-1
        fi

        set -euxo pipefail

        # Generate parental kmer databases if provided:
        trio_arg=""
        if [[ $defined_count -eq 4 ]] ; then
            ./yak count -b~{yak_bloom_filter_size} -t${np} -o maternal.yak <(zcat ~{maternal_fastq_1}) <(zcat ~{maternal_fastq_2})
            ./yak count -b~{yak_bloom_filter_size} -t${np} -o paternal.yak <(zcat ~{paternal_fastq_1}) <(zcat ~{paternal_fastq_2})
            trio_arg="-1 maternal.yak -2 paternal.yak"
        elif [[ $defined_count -eq 2 ]] ; then
            yak count -b~{yak_bloom_filter_size} -t${np} -o maternal.yak ~{maternal_fastq_1}
            yak count -b~{yak_bloom_filter_size} -t${np} -o paternal.yak ~{paternal_fastq_1}
            trio_arg="-1 maternal.yak -2 paternal.yak"
        fi

        time hifiasm \
            -o ~{prefix} \
            -t$((np-1)) \
            -k ~{kmer_size} \
            -f ~{homopolymer_kmer_size} \
            ~{ont_ultralong_reads_arg} \
            ${trio_arg} \
            ~{true="-l0" false="" haploid} \
            ~{true="--n-hap 1" false="" haploid} \
            ~{telomere_5_prime_sequence_arg} \
            ~{reads} \
            ~{extra_args} \
            2>&1 | tee hifiasm.log

        tree -h .

        # Convert all GFAs to fasta files:
        for gfa in ~{prefix}.*.gfa ; do
            bn=$(basename -- "${gfa}" | sed 's#.gfa$##')

            awk '/^S/{print ">"$2; print $3}' \
                "${gfa}" \
                > "${bn}.fa"
        done

        ls -l 

        tar -zcf ~{prefix}.hifiasm.assembleforhaplotigs.tar.gz ~{prefix}.*
    >>>

    output {
        File maternal_phased_gfa = "~{prefix}.dip.hap1.p_ctg.gfa"
        File maternal_phased_fa = "~{prefix}.dip.hap1.p_ctg.fa"
        File paternal_phased_gfa = "~{prefix}.dip.hap2.p_ctg.gfa"
        File paternal_phased_fa = "~{prefix}.dip.hap2.p_ctg.fa"

        File merged_unitigs_gfa = "~{prefix}.dip.p_utg.gfa"
        File merged_unitigs_fa = "~{prefix}.dip.p_utg.fa"

        File all_outputs_gz = "~{prefix}.hifiasm.assembleforhaplotigs.tar.gz"
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
