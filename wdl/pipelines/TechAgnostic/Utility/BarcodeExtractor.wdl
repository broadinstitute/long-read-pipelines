version 1.0
import "../../../structs/Structs.wdl"


workflow BarcodeExtractor {
    meta {
        description: "Extract barcode from BAM file and create movie_name.barcode for each sample."
    }
    input {
        File input_file
        String id
        String movie_name
    }

    call ExtractBarcode {
        input: input_bam = input_file
    }

    output {
        String barcode = ExtractBarcode.barcode
        String movie_name_barcode = movie_name + "." + ExtractBarcode.barcode
    }
}

task ExtractBarcode {
    input {
        File input_bam
        RuntimeAttr? runtime_attr_override
    }

    command <<<
    python3 <<CODE
import re
import sys

def extract_barcode(input_bam):
    match = re.search(r'\\.hifi_reads\\.(bc\\d+)\\.bam', input_bam)
    return match.group(1) if match else None

input_bam = "${input_bam}"
barcode = extract_barcode(input_bam)
if barcode:
    print(barcode)
else:
    print("No barcode found", file=sys.stderr)
CODE
    >>>

    output {
        String barcode = read_string(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-compile-sv-stats:0.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks:                  "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker, default_attr.docker])
    }
}
