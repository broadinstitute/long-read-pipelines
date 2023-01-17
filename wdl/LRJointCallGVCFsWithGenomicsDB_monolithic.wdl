version 1.0

#############################################################################################################
## A workflow that performs joint calling on single-sample gVCFs from GATK4 HaplotypeCaller using GenomicsDB.
#############################################################################################################

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

struct DataTypeParameters {
    Int num_shards
    String map_preset
}

workflow LRJointCallGVCFsWithGenomicsDB {
    input {
        Array[File] gvcfs
        Array[File] gvcf_indices

        File ref_map_file

        File interval_list

        Float snp_filter_level = 99.7
        Array[String] snp_recalibration_annotation_values = ["QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum"]
        Array[Float] snp_recalibration_tranche_values = [100.0, 99.95, 99.9, 99.8, 99.6, 99.5, 99.4, 99.3, 99.0, 98.0, 97.0, 90.0 ]

        Float indel_filter_level = 99.0
        Array[String] indel_recalibration_annotation_values = ["QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum"]
        Array[Float] indel_recalibration_tranche_values = [100.0, 99.95, 99.9, 99.5, 99.0, 97.0, 96.0, 95.0, 94.0, 93.5, 93.0, 92.0, 91.0, 90.0]

        Array[File]?   annotation_bed_files
        Array[File]?   annotation_bed_file_indexes
        Array[String]? annotation_bed_file_annotation_names

        String prefix

        String gcs_out_root_dir
    }

    parameter_meta {
        gvcfs:            "GCS paths to gVCF files"
        gvcf_indices:     "GCS paths to gVCF tbi files"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        prefix:           "prefix for output joint-called gVCF and tabix index"
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/LRJointCallGVCFsWithGenomicsDB/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    # From WARP:
    # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
    # We allow overriding this default behavior for testing / special requests.
    Boolean is_small_callset = length(gvcfs) <= 1000

    # Create sample-name map:
    call CreateSampleNameMap as CreateSampleNameMap {
        input:
            gvcfs = gvcfs,
            prefix = prefix
    }

    # Import our data into GenomicsDB:
    call ImportGVCFs as ImportGVCFsIntoGenomicsDB {
        input:
            sample_name_map = CreateSampleNameMap.sample_name_map,
            interval_list   = interval_list,
            ref_fasta       = ref_map['fasta'],
            ref_fasta_fai   = ref_map['fai'],
            ref_dict        = ref_map['dict'],
            prefix          = prefix,
            batch_size      = 50,
    }

    # Joint call
    call GenotypeGVCFs as JointCallGVCFs {
        input:
            input_gvcf_data = ImportGVCFsIntoGenomicsDB.output_genomicsdb,
            interval_list   = interval_list,
            ref_fasta       = ref_map['fasta'],
            ref_fasta_fai   = ref_map['fai'],
            ref_dict        = ref_map['dict'],
            dbsnp_vcf       = ref_map["known_sites_vcf"],
            prefix          = prefix,
    }

    # First make a sites-only VCF for recal (smaller file, easier to work with):
    call MakeSitesOnlyVcf as MakeSitesOnlyGVCF {
        input:
            vcf = JointCallGVCFs.output_vcf,
            vcf_index = JointCallGVCFs.output_vcf_index,
            prefix = prefix
    }

    # Now we run VariantRecalibrator for indels and snps:
    call IndelsVariantRecalibrator as TrainVQSROnHCIndelVariants {
        input:
            vcf = MakeSitesOnlyGVCF.sites_only_vcf,
            vcf_index = MakeSitesOnlyGVCF.sites_only_vcf_index,
            prefix = prefix + ".indels",
            recalibration_tranche_values = snp_recalibration_tranche_values,
            recalibration_annotation_values = snp_recalibration_annotation_values,
            known_reference_variants = [ref_map["known_sites_vcf"]],
            known_reference_variants_index = [ref_map["known_sites_index"]],
            known_reference_variants_identifier = ["pfcrosses"],
            is_known = [true],
            is_training = [true],
            is_truth = [true],
            prior = [15],
            use_allele_specific_annotations = true,
            max_gaussians = 8,
    }

    call SNPsVariantRecalibratorCreateModel as TrainVQSROnHCSnpVariants {
        input:
            vcf = MakeSitesOnlyGVCF.sites_only_vcf,
            vcf_index = MakeSitesOnlyGVCF.sites_only_vcf_index,
            prefix = prefix + ".snps",
            recalibration_tranche_values = snp_recalibration_tranche_values,
            recalibration_annotation_values = snp_recalibration_annotation_values,
            known_reference_variants = [ref_map["known_sites_vcf"]],
            known_reference_variants_index = [ref_map["known_sites_index"]],
            known_reference_variants_identifier = ["pfcrosses"],
            is_known = [true],
            is_training = [true],
            is_truth = [true],
            prior = [15],
            use_allele_specific_annotations = true,
            max_gaussians = 8,
    }

    call ApplyVqsr as ApplyVqsr {
        input:
            vcf = JointCallGVCFs.output_vcf,
            vcf_index = JointCallGVCFs.output_vcf_index,

            prefix = prefix + ".vqsr_filtered",

            snps_recalibration = TrainVQSROnHCSnpVariants.recalibration,
            snps_recalibration_index = TrainVQSROnHCSnpVariants.recalibration_index,
            snps_tranches = TrainVQSROnHCSnpVariants.tranches,
            snp_filter_level = snp_filter_level,

            indels_recalibration = TrainVQSROnHCIndelVariants.recalibration,
            indels_recalibration_index = TrainVQSROnHCIndelVariants.recalibration_index,
            indels_tranches = TrainVQSROnHCIndelVariants.tranches,
            indel_filter_level = indel_filter_level,

            use_allele_specific_annotations = true,
    }

    # Now we need to annotate our variants by region:
    if (defined(annotation_bed_files)) {
        call AnnotateVcfWithBedRegions as AnnotateVcfRegions {
            input:
                vcf = ApplyVqsr.recalibrated_vcf,
                vcf_index = ApplyVqsr.recalibrated_vcf_index,
                bed_files = select_first([annotation_bed_files]),
                bed_file_indexes = select_first([annotation_bed_file_indexes]),
                bed_file_annotation_names = select_first([annotation_bed_file_annotation_names]),
                prefix = prefix + ".region_annotated"
        }
    }

    # Finally convert the output to a HAIL Matrix Table:
    call ConvertToHailMT as CreateHailMatrixTable {
        input:
            gvcf = select_first([AnnotateVcfRegions.annotated_vcf, ApplyVqsr.recalibrated_vcf]),
            tbi = select_first([AnnotateVcfRegions.annotated_vcf_index, ApplyVqsr.recalibrated_vcf_index]),
            reference = sub(sub(ref_map["fasta"], "^.*/", ""), "\.[fasta]*$", ""),
            ref_fasta = ref_map["fasta"],
            ref_fai = ref_map["fai"],
            prefix = prefix,
            outdir = outdir
    }

    # Finalize:
    File keyfile = select_first([AnnotateVcfRegions.annotated_vcf_index, ApplyVqsr.recalibrated_vcf_index])

    call FinalizeToFile as FinalizeGenomicsDB { input: outdir = outdir, keyfile = keyfile, file = ImportGVCFsIntoGenomicsDB.output_genomicsdb }

    call FinalizeToFile as FinalizeRawVCF { input: outdir = outdir, keyfile = keyfile, file = JointCallGVCFs.output_vcf }
    call FinalizeToFile as FinalizeRawTBI { input: outdir = outdir, keyfile = keyfile, file = JointCallGVCFs.output_vcf_index }

    call FinalizeToFile as FinalizeIndelRecalFile { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.recalibration }
    call FinalizeToFile as FinalizeIndelRecalIndex { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.recalibration_index }
    call FinalizeToFile as FinalizeIndelRecalTranches { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.tranches }
    call FinalizeToFile as FinalizeIndelRecalModelReport { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCIndelVariants.model_report }

    call FinalizeToFile as FinalizeSnpRecalFile { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.recalibration }
    call FinalizeToFile as FinalizeSnpRecalIndex { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.recalibration_index }
    call FinalizeToFile as FinalizeSnpRecalTranches { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.tranches }
    call FinalizeToFile as FinalizeSnpRecalModelReport { input: outdir = outdir, keyfile = keyfile, file = TrainVQSROnHCSnpVariants.model_report }

    call FinalizeToFile as FinalizeVQSRVCF { input: outdir = outdir, keyfile = keyfile, file = ApplyVqsr.recalibrated_vcf }
    call FinalizeToFile as FinalizeVQSRTBI { input: outdir = outdir, keyfile = keyfile, file = ApplyVqsr.recalibrated_vcf_index }

    if (defined(annotation_bed_files)) {
        call FinalizeToFile as FinalizeRegionAnnotatedVcf { input: outdir = outdir, keyfile = keyfile, file = select_first([AnnotateVcfRegions.annotated_vcf]) }
        call FinalizeToFile as FinalizeRegionAnnotatedVcfIndex { input: outdir = outdir, keyfile = keyfile, file = select_first([AnnotateVcfRegions.annotated_vcf_index]) }
    }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File genomicsDB = FinalizeGenomicsDB.gcs_path

        File raw_joint_vcf     = FinalizeRawVCF.gcs_path
        File raw_joint_vcf_tbi = FinalizeRawTBI.gcs_path

        File? vqsr_indel_recal_file         = FinalizeIndelRecalFile.gcs_path
        File? vqsr_indel_recal_file_index   = FinalizeIndelRecalIndex.gcs_path
        File? vqsr_indel_recal_tranches     = FinalizeIndelRecalTranches.gcs_path
        File? vqsr_indel_recal_model_report = FinalizeIndelRecalModelReport.gcs_path

        File? vqsr_snp_recal_file         = FinalizeSnpRecalFile.gcs_path
        File? vqsr_snp_recal_file_index   = FinalizeSnpRecalIndex.gcs_path
        File? vqsr_snp_recal_tranches     = FinalizeSnpRecalTranches.gcs_path
        File? vqsr_snp_recal_model_report = FinalizeSnpRecalModelReport.gcs_path

        File joint_recalibrated_vcf     = FinalizeVQSRVCF.gcs_path
        File joint_recalibrated_vcf_tbi = FinalizeVQSRTBI.gcs_path

        File? annotated_joint_vcf     = AnnotateVcfRegions.annotated_vcf
        File? annotated_joint_vcf_tbi = AnnotateVcfRegions.annotated_vcf_index

        File joint_mt = CreateHailMatrixTable.gcs_path
    }
}


task CreateSampleNameMap {

    meta {
        description: "Creates the sample / name-map file of the GVCFs for ingest into ImportGVCFs.  NOTE: Some of this functionality is duplicated from Utils.InferSampleName.  This is intentional - we don't want to localize all these files or shard over potentially thousands of input GVCFs."
    }

    input {
        Array[File] gvcfs
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        gvcfs: {
            help: "Array of single-sample GVCF files.",
            localization_optional: true
        }
    }

    Int disk_size_gb = 20

    String outfile_name = "~{prefix}.sample_name_map.tsv"

    # Every so often we should reauthorize so `bcftools` can continue to access our data:
    Int re_auth_interval = 50

    command <<<
        set -euxo pipefail

        # Put our gvcfs into a file we can iterate over:
        gvcf_file_list=~{write_lines(gvcfs)}

        # Initialize a file for the sample names:
        [ -e ~{outfile_name} ] && rm -rf ~{outfile_name}

        # Set our access token:
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        let i=1
        while read file_path ; do

            # Get our sample list from our file:
            bcftools query -l ${file_path} > sample_names.txt

            # Make sure we only have one sample name:
            [[ $(wc -l sample_names.txt | awk '{print $1}') -ne 1 ]] && echo "Incorrect number of sample names found in GVCF (there can be only one!): ${file_path}" && exit 1

            # Make sure the samplename has an actual name:
            [ $(grep -iq "unnamedsample" sample_names.txt) ] && echo "Sample name found to be unnamedsample in GVCF: ${file_path}" && exit 1

            # Add the sample name and GVCF path to the sample name file:
            echo -e "$(cat sample_names.txt)\t${file_path}" >> ~{outfile_name}

            let i=$i+1
            if [[ $i -gt ~{re_auth_interval} ]] ; then
                # Periodically we should update the token so we don't have problems with long file lists:
                export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
                i=0
            fi
        done < ${gvcf_file_list}
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
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

    output {
        File sample_name_map = outfile_name
    }
}

task ImportGVCFs {

    input {
        File sample_name_map

        File interval_list

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix

        Int batch_size = 50

        RuntimeAttr? runtime_attr_override
    }

    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB"))

    Int disk_size = 1 + 4*ref_size

    command <<<
        set -euxo pipefail

        # Make sure that the output directory does not exist:
        [ -e ~{prefix} ] && rm -rf ~{prefix}

        #
        # Notes from WARP Team:
        #
        # We've seen some GenomicsDB performance regressions related to intervals, so we're going to pretend we only have a single interval
        # using the --merge-input-intervals arg
        # There's no data in between since we didn't run HaplotypeCaller over those loci so we're not wasting any compute

        # The memory setting here is very important and must be several GiB lower
        # than the total memory allocated to the VM because this tool uses
        # a significant amount of non-heap memory for native libraries.
        # Also, testing has shown that the multithreaded reader initialization
        # does not scale well beyond 5 threads, so don't increase beyond that.
        gatk --java-options "-Xms8000m -Xmx25000m" \
            GenomicsDBImport \
                --genomicsdb-workspace-path ~{prefix}.genomicsDB \
                --batch-size ~{batch_size} \
                -L ~{interval_list} \
                --sample-name-map ~{sample_name_map} \
                --reader-threads 5 \
                --merge-input-intervals \
                --consolidate

        tar -cf ~{prefix}.genomicsDB.tar ~{prefix}.genomicsDB
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.3.0.0"
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

    output {
        File output_genomicsdb = "~{prefix}.genomicsDB.tar"
    }
}

task GenotypeGVCFs {

    input {
        File input_gvcf_data
        File? input_gvcf_index  # Required if passing a VCF file.

        File interval_list

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String dbsnp_vcf

        String prefix

        Boolean keep_combined_raw_annotations = false
        RuntimeAttr? runtime_attr_override
    }

    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB"))
    Int db_snp_size = ceil(size(dbsnp_vcf, "GB"))

    Int disk_size = 1 + 4*ceil(size(input_gvcf_data, "GB")) + ref_size + db_snp_size

    parameter_meta {
        input_gvcf_data: { help: "Either a single GVCF file or a GenomicsDB Tar file." }
        interval_list: {
            localization_optional: true
        }
    }

    command <<<
        set -euxo pipefail

        # We must determine if our input variants are in a genomicsdb file or in a VCF.
        # The easiest way is to see if the input is a .tar file:

        is_genomics_db=true
        filename=$(basename -- "~{input_gvcf_data}")
        extension="${filename##*.}"
        if [[ "${extension}" != "tar" ]] ; then
            is_genomics_db=false
        fi

        if $is_genomics_db ; then
            tar -xf ~{input_gvcf_data}
            INPUT_FILE="gendb://$(basename ~{input_gvcf_data} .tar)"
        else
            INPUT_FILE=~{input_gvcf_data}
        fi

        gatk --java-options "-Xms8000m -Xmx25000m" \
            GenotypeGVCFs \
                -R ~{ref_fasta} \
                -O ~{prefix}.vcf.gz \
                -D ~{dbsnp_vcf} \
                -G StandardAnnotation -G AS_StandardAnnotation \
                --only-output-calls-starting-in-intervals \
                -V ${INPUT_FILE} \
                -L ~{interval_list} \
                ~{true='--keep-combined-raw-annotations' false='' keep_combined_raw_annotations} \
                --merge-input-intervals
    >>>
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             26,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.3.0.0"
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

    output {
        File output_vcf = "~{prefix}.vcf.gz"
        File output_vcf_index = "~{prefix}.vcf.gz.tbi"
    }
}


task HardFilterVcf {

    input {
        File vcf
        File vcf_index

        String prefix

        # From WARP:
        # ExcessHet is a phred-scaled p-value. We want a cutoof anything more extreme
        # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
        Float excess_het_threshold = 54.69

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size([vcf, vcf_index], "GB"))

    command <<<
        set -euo pipefail

        # Get amount of memory to use:
        mem_available=$(free -m | grep '^Mem' | awk '{print $2}')
        let mem_start=${mem_available}-1000
        let mem_max=${mem_available}-750

        gatk --java-options "-Xms${mem_start}m -Xmx${mem_max}m" \
            VariantFiltration \
            --filter-expression "ExcessHet > ~{excess_het_threshold}" \
            --filter-name ExcessHet \
            -V ~{vcf} \
            -O ~{prefix}.hard_filtered.vcf.gz
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

    output {
        File variant_filtered_vcf = "~{prefix}.hard_filtered.vcf.gz"
        File variant_filtered_vcf_index = "~{prefix}.hard_filtered.vcf.gz.tbi"
    }
}

task MakeSitesOnlyVcf {

    input {
        File vcf
        File vcf_index

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size([vcf, vcf_index], "GB"))

    command <<<
        set -euo pipefail

        # Get amount of memory to use:
        mem_available=$(free -m | grep '^Mem' | awk '{print $2}')
        let mem_start=${mem_available}-1000
        let mem_max=${mem_available}-750

        gatk --java-options "-Xms${mem_start}m -Xmx${mem_max}m" \
            MakeSitesOnlyVcf \
            -I ~{vcf} \
            -O ~{prefix}.sites_only.vcf.gz
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

    output {
        File sites_only_vcf = "~{prefix}.sites_only.vcf.gz"
        File sites_only_vcf_index = "~{prefix}.sites_only.vcf.gz.tbi"
    }
}

task AnnotateVcfWithBedRegions {
    input {
        File vcf
        File vcf_index

        Array[File] bed_files
        Array[File] bed_file_indexes
        Array[String] bed_file_annotation_names

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size([vcf, vcf_index, bed_files, bed_file_indexes], "GB"))

    command <<<
        set -euxo pipefail

        # Get amount of memory to use:
        mem_available=$(free -m | grep '^Mem' | awk '{print $2}')
        let mem_start=${mem_available}-1000
        let mem_max=${mem_available}-750

        # We need to generate argument strings from the input arrays.
        # First we check that the arrays are the same length:
        if [[ ~{length(bed_files)} -ne ~{length(bed_file_indexes)} ]] || \
           [[ ~{length(bed_files)} -ne ~{length(bed_file_annotation_names)} ]] ; then
            echo "ERROR: Not all input arrays for known variants contain the same number of elements: " 1>&2
            echo "       bed_files                 = ~{length(bed_files)}" 1>&2
            echo "       bed_file_indices          = ~{length(bed_file_indexes)}" 1>&2
            echo "       bed_file_annotation_names = ~{length(bed_file_annotation_names)}" 1>&2
            false
        fi

        # Now we can write out the arrays into a TSV file and add them line by line to the execution:
        # Create the TSV:
        options_tsv=~{write_tsv(transpose([bed_files, bed_file_annotation_names]))}

        # Now we have to run `VariantFiltration` multiple times on its own output so that it can
        # annotate each region in the file:
        # NOTE: This is dumb, but must be done because the `--mask` and `--mask-name` inputs are not arrays.

        input_vcf=~{vcf}
        output_vcf=~{prefix}.intermediate.vcf.gz
        while read mask_options ; do

            bed_file=$(echo "${mask_options}" | awk -F'\t' '{print $1}')
            mask_name=$(echo "${mask_options}" | awk -F'\t' '{print $2}')

            echo -e "RUNNING GATK ON NEW MASK: ${mask_name}\t${bed_file}"

            gatk --java-options "-Xms${mem_start}m -Xmx${mem_max}m" \
                VariantFiltration \
                -V ${input_vcf} \
                -O ${output_vcf} \
                --mask ${bed_file} \
                --mask-name ${mask_name}

            mv ${output_vcf} ~{prefix}.new_input.vcf.gz
            mv ${output_vcf}.tbi ~{prefix}.new_input.vcf.gz.tbi
            input_vcf=~{prefix}.new_input.vcf.gz
        done < ${options_tsv}

        # Because of the `mv` at the end of the loop we need to move the "new_input" files here:
        mv ~{prefix}.new_input.vcf.gz ~{prefix}.vcf.gz
        mv ~{prefix}.new_input.vcf.gz.tbi ~{prefix}.vcf.gz.tbi
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_index = "~{prefix}.vcf.gz.tbi"
    }
}

task IndelsVariantRecalibrator {

    input {
        File vcf
        File vcf_index

        String prefix

        Array[String] recalibration_tranche_values
        Array[String] recalibration_annotation_values

        Array[File] known_reference_variants
        Array[File] known_reference_variants_index
        Array[String] known_reference_variants_identifier
        Array[Boolean] is_known
        Array[Boolean] is_training
        Array[Boolean] is_truth
        Array[Int] prior

        Boolean use_allele_specific_annotations
        Int max_gaussians = 4

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        vcf:   "Sites only VCF.  Can be pre-filtered using hard-filters."
        vcf_index: "Tribble Index for sites only VCF."
        known_reference_variants: "Array of known reference VCF files.  For humans, dbSNP is one example."
        known_reference_variants_index: "Array of index files for known reference VCF files."
        known_reference_variants_identifier: "Array of boolean values the identifier / name for the known_reference_variant file at the same array position.  Must be the same length as `known_reference_variants`."
        is_known: "Array of boolean values indicating if the known_reference_variant file at the same array position contains known variants.  Must be the same length as `known_reference_variants`."
        is_training: "Array of boolean values indicating if the known_reference_variant file at the same array position contains training data.  Must be the same length as `known_reference_variants`."
        is_truth: "Array of boolean values indicating if the known_reference_variant file at the same array position contains truth data.  Must be the same length as `known_reference_variants`."
        prior: "Array of integer values indicating the priors for the known_reference_variant file at the same array position.  Must be the same length as `known_reference_variants`."
    }


    Int disk_size = 10 + ceil(size(known_reference_variants, "GB"))
                  + 4*ceil(size(vcf, "GB"))
                  + 2*ceil(size(vcf_index, "GB"))

    command <<<
        set -euxo pipefail

        # We need to generate resource strings from the input arrays.
        # First we check that the arrays are the same length:
        if [[ ~{length(known_reference_variants)} -ne ~{length(known_reference_variants_identifier)} ]] || \
           [[ ~{length(known_reference_variants)} -ne ~{length(known_reference_variants_index)} ]] || \
           [[ ~{length(known_reference_variants)} -ne ~{length(is_known)} ]] || \
           [[ ~{length(known_reference_variants)} -ne ~{length(is_training)} ]] || \
           [[ ~{length(known_reference_variants)} -ne ~{length(is_truth)} ]] || \
           [[ ~{length(known_reference_variants)} -ne ~{length(prior)} ]] ; then
            echo "ERROR: Not all input arrays for known variants contain the same number of elements: " 1>&2
            echo "       known_reference_variants            = ~{length(known_reference_variants)}" 1>&2
            echo "       known_reference_variants            = ~{length(known_reference_variants_index)}" 1>&2
            echo "       known_reference_variants_identifier = ~{length(known_reference_variants_identifier)}" 1>&2
            echo "       is_known                            = ~{length(is_known)}" 1>&2
            echo "       is_training                         = ~{length(is_training)}" 1>&2
            echo "       is_truth                            = ~{length(is_truth)}" 1>&2
            echo "       prior                               = ~{length(prior)}" 1>&2
            false
        fi

        # Now we can write out the arrays into a TSV file and add them line by line to the execution:
        # Create the TSV:
        options_tsv=~{write_tsv(transpose([known_reference_variants_identifier, is_known, is_training, is_truth, prior, known_reference_variants]))}

        # Now read them into a string:
        resource_flags=$(awk '{printf("--resource:%s,known=%s,training=%s,truth=%s,prior=%d %s ", $1, $2, $3, $4, $5, $6)}' ${options_tsv})

        # Get amount of memory to use:
        mem_available=$(free -g | grep '^Mem' | awk '{print $2}')
        let mem_start=${mem_available}-2
        let mem_max=${mem_available}-1

        gatk --java-options "-Xms${mem_start}g -Xmx${mem_max}g" \
            VariantRecalibrator \
                -V ~{vcf} \
                -O ~{prefix}.recal \
                --tranches-file ~{prefix}.tranches \
                --trust-all-polymorphic \
                -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
                -an ~{sep=' -an ' recalibration_annotation_values} \
                ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
                -mode INDEL \
                --output-model ~{prefix}.model.report \
                --max-gaussians ~{max_gaussians} \
                ${resource_flags}
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             26,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

    output {
        File recalibration = "~{prefix}.recal"
        File recalibration_index = "~{prefix}.recal.idx"
        File tranches = "~{prefix}.tranches"
        File model_report = "~{prefix}.model.report"
    }
}

task SNPsVariantRecalibratorCreateModel {

    input {
        File vcf
        File vcf_index

        String prefix

        Array[String] recalibration_tranche_values
        Array[String] recalibration_annotation_values

        Array[File] known_reference_variants
        Array[File] known_reference_variants_index
        Array[String] known_reference_variants_identifier
        Array[Boolean] is_known
        Array[Boolean] is_training
        Array[Boolean] is_truth
        Array[Int] prior

        Int? downsampleFactor

        Boolean use_allele_specific_annotations
        Int max_gaussians = 6

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        vcf:   "Sites only VCF.  Can be pre-filtered using hard-filters."
        vcf_index: "Tribble Index for sites only VCF."
        known_reference_variants: "Array of known reference VCF files.  For humans, dbSNP is one example."
        known_reference_variants_index: "Array of index files for known reference VCF files."
        known_reference_variants_identifier: "Array of boolean values the identifier / name for the known_reference_variant file at the same array position.  Must be the same length as `known_reference_variants`."
        is_known: "Array of boolean values indicating if the known_reference_variant file at the same array position contains known variants.  Must be the same length as `known_reference_variants`."
        is_training: "Array of boolean values indicating if the known_reference_variant file at the same array position contains training data.  Must be the same length as `known_reference_variants`."
        is_truth: "Array of boolean values indicating if the known_reference_variant file at the same array position contains truth data.  Must be the same length as `known_reference_variants`."
        prior: "Array of integer values indicating the priors for the known_reference_variant file at the same array position.  Must be the same length as `known_reference_variants`."
    }

    Int disk_size = 10 + ceil(size(known_reference_variants, "GB"))
              + 4*ceil(size(vcf, "GB"))
              + 2*ceil(size(vcf_index, "GB"))

    String downsample_factor_arg = if defined(downsampleFactor) then " --sample-every-Nth-variant " else ""

    command <<<
        set -euxo pipefail

        # We need to generate resource strings from the input arrays.
        # First we check that the arrays are the same length:
        if [[ ~{length(known_reference_variants)} -ne ~{length(known_reference_variants_identifier)} ]] || \
           [[ ~{length(known_reference_variants)} -ne ~{length(known_reference_variants_index)} ]] || \
           [[ ~{length(known_reference_variants)} -ne ~{length(is_known)} ]] || \
           [[ ~{length(known_reference_variants)} -ne ~{length(is_training)} ]] || \
           [[ ~{length(known_reference_variants)} -ne ~{length(is_truth)} ]] || \
           [[ ~{length(known_reference_variants)} -ne ~{length(prior)} ]] ; then
            echo "ERROR: Not all input arrays for known variants contain the same number of elements: " 1>&2
            echo "       known_reference_variants            = ~{length(known_reference_variants)}" 1>&2
            echo "       known_reference_variants            = ~{length(known_reference_variants_index)}" 1>&2
            echo "       known_reference_variants_identifier = ~{length(known_reference_variants_identifier)}" 1>&2
            echo "       is_known                            = ~{length(is_known)}" 1>&2
            echo "       is_training                         = ~{length(is_training)}" 1>&2
            echo "       is_truth                            = ~{length(is_truth)}" 1>&2
            echo "       prior                               = ~{length(prior)}" 1>&2
            false
        fi

        # Now we can write out the arrays into a TSV file and add them line by line to the execution:
        # Create the TSV:
        options_tsv=~{write_tsv(transpose([known_reference_variants_identifier, is_known, is_training, is_truth, prior, known_reference_variants]))}

        # Now read them into a string:
        resource_flags=$(awk '{printf("--resource:%s,known=%s,training=%s,truth=%s,prior=%d %s ", $1, $2, $3, $4, $5, $6)}' ${options_tsv})

        # Get amount of memory to use:
        mem_available=$(free -g | grep '^Mem' | awk '{print $2}')
        let mem_start=${mem_available}-2
        let mem_max=${mem_available}-1

        gatk --java-options "-Xms${mem_start}g -Xmx${mem_max}g" \
            VariantRecalibrator \
                -V ~{vcf} \
                -O ~{prefix}.recal \
                --tranches-file ~{prefix}.tranches \
                --trust-all-polymorphic \
                -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
                -an ~{sep=' -an ' recalibration_annotation_values} \
                ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
                -mode SNP \
                ~{downsample_factor_arg}~{default="" sep=" --sample-every-Nth-variant " downsampleFactor} \
                --output-model ~{prefix}.model.report \
                --max-gaussians ~{max_gaussians} \
                ${resource_flags}
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

    output {
        File recalibration = "~{prefix}.recal"
        File recalibration_index = "~{prefix}.recal.idx"
        File tranches = "~{prefix}.tranches"
        File model_report = "~{prefix}.model.report"
    }
}

task ApplyVqsr {

    input {
        File vcf
        File vcf_index

        String prefix

        File snps_recalibration
        File snps_recalibration_index
        File snps_tranches
        Float snp_filter_level

        File indels_recalibration
        File indels_recalibration_index
        File indels_tranches
        Float indel_filter_level

        Boolean use_allele_specific_annotations

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(size([vcf, vcf_index], "GB"))
          + 2*ceil(size([snps_recalibration, snps_recalibration_index, snps_tranches], "GB"))
          + 2*ceil(size([indels_recalibration, indels_recalibration_index, indels_tranches], "GB"))

    command <<<
        set -euxo pipefail

        # Get amount of memory to use:
        mem_available=$(free -m | grep '^Mem' | awk '{print $2}')
        let mem_start=${mem_available}-2000
        let mem_max=${mem_available}-500

        gatk --java-options "-Xms${mem_start}m -Xmx${mem_max}m" \
            ApplyVQSR \
                -V ~{vcf} \
                -O tmp.indel.recalibrated.vcf.gz \
                --recal-file ~{indels_recalibration} \
                ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
                --tranches-file ~{indels_tranches} \
                --truth-sensitivity-filter-level ~{indel_filter_level} \
                --create-output-variant-index true \
                -mode INDEL

        gatk --java-options "-Xms${mem_start}m -Xmx${mem_max}m" \
            ApplyVQSR \
                -V tmp.indel.recalibrated.vcf.gz \
                -O ~{prefix}.recalibrated.vcf.gz \
                --recal-file ~{snps_recalibration} \
                ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
                --tranches-file ~{snps_tranches} \
                --truth-sensitivity-filter-level ~{snp_filter_level} \
                --create-output-variant-index true \
                -mode SNP
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             7,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

    output {
        File recalibrated_vcf = "~{prefix}.recalibrated.vcf.gz"
        File recalibrated_vcf_index = "~{prefix}.recalibrated.vcf.gz.tbi"
    }
}


task ConvertToHailMT {
    meta {
        description: "Convert a .vcf.bgz file to a Hail MatrixTable and copy it to a final gs:// URL."
    }

    input {
        File gvcf
        File tbi

        String reference = "GRCh38"
        String? ref_fasta
        String? ref_fai
        String prefix = "out"

        String outdir

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 3*ceil(size(gvcf, "GB"))

    command <<<
        set -x

        python3 <<EOF

        import hail as hl
        hl.init(default_reference='GRCh38')

        if '~{defined(ref_fasta)}' == 'true' and '~{defined(ref_fai)}' == 'true':
            ref = hl.ReferenceGenome.from_fasta_file('~{reference}', '~{ref_fasta}', '~{ref_fai}')

        callset = hl.import_vcf(
            '~{gvcf}',
            array_elements_required=False,
            force_bgz=True,
            reference_genome='~{reference}'
        )

        callset.write('~{prefix}.mt')

        EOF

        gsutil -m rsync -Cr ~{prefix}.mt ~{outdir}/~{prefix}.mt
    >>>

    output {
        String gcs_path = "~{outdir}/~{prefix}.mt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "hailgenetics/hail:0.2.105"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FinalizeToFile {
    input {
        File file
        String outdir
        String? name

        File? keyfile

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        file: {
            description: "file to finalize",
            localization_optional: true
        }
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
        outdir: "directory to which files should be uploaded"
        name:   "name to set for uploaded file"
    }

    String gcs_output_dir = sub(outdir, "/+$", "")
    String gcs_output_file = gcs_output_dir + "/" + select_first([name, basename(file)])

    command <<<
        set -euxo pipefail

        gsutil -m cp "~{file}" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
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

