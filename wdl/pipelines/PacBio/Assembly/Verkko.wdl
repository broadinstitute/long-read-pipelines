version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/Utility/Utils.wdl"

workflow Verkko {

    meta {
        description: "Perform a genome assembly using Verkko2."
    }
    parameter_meta {
        pacbio_hifi_reads: "PacBio HiFi reads"
        sample_name: "Sample name"
        is_haploid: "Whether the sample is haploid (default: false)"
        
        nanopore_scaffolding_read_basecall_dir: "Nanopore scaffolding reads basecall directory"

        maternal_fastq_files: "Maternal fastq files for trio assembly from which to generate hapmer databases (optional)"
        paternal_fastq_files: "Paternal fastq files for trio assembly from which to generate hapmer databases (optional)"
        maternal_hapmer_database_tar_gz: "Maternal hapmer database tar.gz (optional)"
        paternal_hapmer_database_tar_gz: "Paternal hapmer database tar.gz (optional)"
        hap_kmers_type: "Hapmer database type for trio assembly (optional)"
    }
    input {
        Array[File] pacbio_hifi_reads
        String sample_name
        Boolean is_haploid = false
        
        String? nanopore_scaffolding_read_basecall_dir

        Array[File]? maternal_fastq_files
        Array[File]? paternal_fastq_files
        File? maternal_hapmer_database_tar_gz
        File? paternal_hapmer_database_tar_gz
        String? hap_kmers_type
    }

    # Validate that either all hapmer database options are given or none of them are given:
    Int hapmer_database_defined_count = (if defined(maternal_hapmer_database_tar_gz) then 1 else 0) + 
                        (if defined(paternal_hapmer_database_tar_gz) then 1 else 0) +
                        (if defined(hap_kmers_type) then 1 else 0)


    # Validate that either all hapmer fastq file set options are given or none of them are given:
    Int hapmer_fastq_defined_count = (if defined(maternal_fastq_files) then 1 else 0) + 
                        (if defined(paternal_fastq_files) then 1 else 0) +
                        (if defined(hap_kmers_type) then 1 else 0)

    if ((hapmer_fastq_defined_count != 3) && (hapmer_fastq_defined_count != 0)) {
        call Utils.StopWorkflow as t_000_StopWorkflow_InvalidHapmeFastqOptions {
            input:
                reason = "Either all or none of the parental hapmer fastq file sets must be defined.  If both are defined, the hap_kmers_type must be defined as well."
        }
    }
    if ((hapmer_database_defined_count != 3) && (hapmer_database_defined_count != 0)) {
        call Utils.StopWorkflow as t_000_StopWorkflow_InvalidHapmerDatabaseOptions {
            input:
                reason = "Either all or none of the parental hapmer database file sets must be defined.  If both are defined, the hap_kmers_type must be defined as well."
        }
    }

    if ((hapmer_database_defined_count == 3) && (hapmer_fastq_defined_count == 3)) {
        call Utils.StopWorkflow as t_000_StopWorkflow_InvalidHapmerDatabaseAndFastqOptions {
            input:
                reason = "You must either provide hapmer fastq file sets or hapmer database files, but not both."
        }
    }

    if ((hapmer_fastq_defined_count == 3) && (hapmer_database_defined_count != 3)) {
        call CreateParentalHapmerDatabases as t_001_CreateParentalHapmerDatabases {
            input:
                maternal_fastq_files = select_first([maternal_fastq_files]),
                paternal_fastq_files = select_first([paternal_fastq_files]),
                prefix = sample_name
        }
    }

    if ((hapmer_database_defined_count == 3) || (hapmer_fastq_defined_count == 3)) {
        File final_maternal_hapmer_database_tar_gz = select_first([maternal_hapmer_database_tar_gz, t_001_CreateParentalHapmerDatabases.maternal_hapmer_database_tar_gz])
        File final_paternal_hapmer_database_tar_gz = select_first([paternal_hapmer_database_tar_gz, t_001_CreateParentalHapmerDatabases.paternal_hapmer_database_tar_gz])
    }

    call VerkoAssemble as t_002_VerkkoAssemble {
        input:
            pacbio_hifi_reads = pacbio_hifi_reads,
            nanopore_scaffolding_read_basecall_dir = nanopore_scaffolding_read_basecall_dir,
            prefix = sample_name,
            maternal_hapmer_database_tar_gz = final_maternal_hapmer_database_tar_gz,
            paternal_hapmer_database_tar_gz = final_paternal_hapmer_database_tar_gz,
            hap_kmers_type = hap_kmers_type
    }

    output {
        File assembly_fasta = t_002_VerkkoAssemble.assembly_fasta
        File output_dir_tar_gz = t_002_VerkkoAssemble.output_dir_tar_gz    
    }
}

task VerkoAssemble {
    input {
        Array[File] pacbio_hifi_reads
        String? nanopore_scaffolding_read_basecall_dir
        String prefix

        Boolean is_haploid = true

        File? maternal_hapmer_database_tar_gz
        File? paternal_hapmer_database_tar_gz
        String? hap_kmers_type

        String extra_args = ""

        Int default_preemptible_tries = 0
        Int default_max_retries = 0
        RuntimeAttr? runtime_attr_override
    }

    String out_folder_name = "~{prefix}_assembly_verkko"

    String hap_kmers_arg = if (defined(maternal_hapmer_database_tar_gz) && defined(paternal_hapmer_database_tar_gz) && defined(hap_kmers_type)) then "--hap-kmers " + basename(select_first([maternal_hapmer_database_tar_gz]), ".tar.gz") + " " + basename(select_first([paternal_hapmer_database_tar_gz]), ".tar.gz") + " " + hap_kmers_type else ""

    # Nanopore scaffolding is about 30GB, but we don't know how much space it will take up.
    Int disk_size = 10 + 2*(
        2 * ceil(size(pacbio_hifi_reads, "GB")) + 
        30 + 
        11 * ceil(size(maternal_hapmer_database_tar_gz, "GB")) + 
        11 * ceil(size(paternal_hapmer_database_tar_gz, "GB"))
    )

    command <<<
        # Check if we have the parental hapmer databases.
        # They are either all defined or all undefined:
        defined_count=0
        if [[ -n "~{maternal_hapmer_database_tar_gz}" ]]; then ((defined_count++)); fi
        if [[ -n "~{paternal_hapmer_database_tar_gz}" ]]; then ((defined_count++)); fi
        if [[ -n "~{hap_kmers_type}" ]]; then ((defined_count++)); fi
        if [[ $defined_count -eq 3 ]]; then
            tar -xzf ~{maternal_hapmer_database_tar_gz}
            tar -xzf ~{paternal_hapmer_database_tar_gz}
        elif [[ $defined_count -ne 0 ]]; then
            echo "Error: Either both or none of the parental hapmer databases must be defined.  If both are defined, the hap_kmers_type must be defined as well." 1>&2
            echo "       Maternal hapmer database: ~{maternal_hapmer_database_tar_gz}" 1>&2
            echo "       Paternal hapmer database: ~{paternal_hapmer_database_tar_gz}" 1>&2
            echo "       Hap_kmers type: ~{hap_kmers_type}" 1>&2
            exit 1
        fi

        # Setup nanopore scaffolding if we have any:
        if [[ -n "~{nanopore_scaffolding_read_basecall_dir}" ]]; then
            mkdir nanopore_scaffolding
            gsutil -m cp -r "~{nanopore_scaffolding_read_basecall_dir}/*" nanopore_scaffolding/
            fastq_file_exemplar=$(find nanopore_scaffolding -type f -name "*fastq*" | head -n1)
            nanopore_scaffolding_fastq_folder=$(dirname ${fastq_file_exemplar})
            nanopore_scaffolding_arg="--nano ${nanopore_scaffolding_fastq_folder}"
        else
            nanopore_scaffolding_arg=""
        fi

        ############################################################
        set -euxo pipefail

        time verkko \
            -d ~{out_folder_name} \
            --hifi ~{sep=' ' pacbio_hifi_reads} \
            ${nanopore_scaffolding_arg} \
            ~{true="--haploid" false="" is_haploid} \
            ~{hap_kmers_arg}

        tar -czf ~{out_folder_name}.tar.gz ~{out_folder_name}
    >>>

    output {
        File assembly_fasta = "~{out_folder_name}/~{prefix}.fasta"
        File output_dir_tar_gz = "~{out_folder_name}.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          12,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  default_preemptible_tries,
        max_retries:        default_max_retries,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-verkko:2.2.1"
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


task CreateParentalHapmerDatabases {
    input {
        Array[File] maternal_fastq_files
        Array[File] paternal_fastq_files
        String prefix
        
        Int kmer_size = 30

        Boolean use_compressed_space = true

        String extra_args = ""

        Int default_preemptible_tries = 0
        Int default_max_retries = 0
        RuntimeAttr? runtime_attr_override
    }

    String out_base_name = "hapmer_database" + if (use_compressed_space) then "_compress" else "" + ".k~{kmer_size}"

    Int disk_size = 10 + 2*(2 * ceil(size(maternal_fastq_files, "GB")) + 2 * ceil(size(paternal_fastq_files, "GB")))

    command <<<
        # Make sure we use all our processors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            let np=${np}-1
        fi

        $ Make sure we use all our memory, too:
        memory_gb=$(free -h | awk '{print $2}' | tail -n+2 | head -1 | sed 's@\([0-9]*\)[A-z]*@\1@')
        memory_gb=$((memory_gb-1))

        set -euxo pipefail

        time meryl count compress k=30 threads=$((np-1)) memory=${memory_gb} ~{sep=' ' maternal_fastq_files} output ~{prefix}_maternal_~{out_base_name}.meryl
        time meryl count compress k=30 threads=$((np-1)) memory=${memory_gb} ~{sep=' ' paternal_fastq_files} output ~{prefix}_paternal_~{out_base_name}.meryl

        $MERQURY/trio/hapmers.sh ~{prefix}_maternal_~{out_base_name}.meryl ~{prefix}_paternal_~{out_base_name}.meryl

        tar -czf ~{prefix}_maternal_~{out_base_name}.meryl.tar.gz ~{prefix}_maternal_~{out_base_name}.meryl
        tar -czf ~{prefix}_paternal_~{out_base_name}.meryl.tar.gz ~{prefix}_paternal_~{out_base_name}.meryl
    >>>

    output {
        File maternal_hapmer_database_tar_gz = "~{prefix}_maternal_~{out_base_name}.meryl.tar.gz"
        File paternal_hapmer_database_tar_gz = "~{prefix}_paternal_~{out_base_name}.meryl.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          12,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  default_preemptible_tries,
        max_retries:        default_max_retries,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-verkko:2.2.1"
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
