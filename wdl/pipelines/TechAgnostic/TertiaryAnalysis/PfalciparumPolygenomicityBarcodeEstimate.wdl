version 1.0

import "../../../tasks/Utility/VariantUtils.wdl" as VARUTIL

workflow PfalciparumPolygenomicityBarcodeEstimate {
    meta {
        author: "Jonn Smith"
        description: "This workflow estimates whether a given sample is polygenomic or monogenomic based on a given fingerprint / barcode.  Fingerprint is assumed to be 24 SNP barcode from Daniels et al.  If no fingerprint / barcode is given, the fingerprint will be extracted from the given VCF."
    }

    parameter_meta {
        sample_name: "Name of the sample to be processed."
        barcode_string: "Barcode to use for genomicity estimation.  Assumed to be 24 SNP barcode from Daniels et al."

        ref_map_file:  "Reference map file indicating reference sequence and auxillary file locations."
        vcf: "VCF file from which to extract the barcode."
        vcf_index: "Index for the VCF file."
        fingerprint_haploytpe_db_file: "Haplotype DB file from which to fingerprint the input data."
    }

    input {
        String sample_name
        String? barcode_string
        String? ref_map_file
        File? vcf
        File? vcf_index
        File? fingerprint_haploytpe_db_file
    }

    # Input error checking:
    if (!defined(barcode_string) &&
          (!defined(fingerprint_haploytpe_db_file) &&
            !defined(ref_map_file) && !defined(vcf) && !defined(vcf_index))
        ) {
        call FailWithMessage{input: message="If barcode is not provided, then fingerprint_haploytpe_db_file, ref_map_file, vcf, and vcf_index must be provided."}
    }

    if (!defined(barcode_string)) {
        Map[String, String] ref_map = read_map(select_first([ref_map_file]))

        call VARUTIL.ExtractFingerprintAndBarcode as t_001_FingerprintAndBarcodeVcf {
            input:
                vcf = select_first([vcf]),
                vcf_index = select_first([vcf_index]),
                haplotype_database_file = select_first([fingerprint_haploytpe_db_file]),
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                prefix = sample_name
        }
    }

    call EstimateGenomicityFromBarcode as t_002_EstimateGenomicityFromBarcode {
        input:
            barcode = select_first([barcode_string, t_001_FingerprintAndBarcodeVcf.barcode])
    }

    output {
        ########################################
        # Primary outputs:
        String genomicity = t_002_EstimateGenomicityFromBarcode.genomicity

        File? fingerprint_vcf = t_001_FingerprintAndBarcodeVcf.output_vcf
        String? barcode = t_001_FingerprintAndBarcodeVcf.barcode
    }
}

task EstimateGenomicityFromBarcode {
    input {
        String barcode

        Int hets_for_poly = 2

        RuntimeAttr? runtime_attr_override
    }

    String genomicity_file = "genomicity_estimate.txt"

    command <<<
        set -euxo pipefail
        python3 <<CODE

        with open("~{genomicity_file}", "w") as f:
            if "~{barcode}".count("N") > ~{hets_for_poly}:
                f.write("polygenomic")
            else:
                f.write("monogenomic")
        CODE
    >>>

    output {
        String genomicity = read_string(genomicity_file)
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "python:3.9.19-slim-bullseye"
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

task FailWithMessage {
    input {
        String message
    }

    command <<<
        set -euxo pipefail
        echo "~{message}"
        exit 1
    >>>

    output {}

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "alpine:3.20.2"
    }
    runtime {
        cpu:                    select_first([default_attr.cpu_cores])
        memory:                 select_first([default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([default_attr.boot_disk_gb])
        preemptible:            select_first([default_attr.preemptible_tries])
        maxRetries:             select_first([default_attr.max_retries])
        docker:                 select_first([default_attr.docker])
    }
}