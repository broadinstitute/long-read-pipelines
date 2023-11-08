version 1.0

import "../../../tasks/VariantCalling/Sniffles2.wdl" as SNF
import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/Utility/Utils.wdl"

workflow Sniffles2JointCallSV {
	input {
        File? list_of_snfs
        Array[File]? snfs

        String out_prefix

        Int cores
        String disk_type = "HDD"

        String gcs_out_root_dir
    }
    output {
    	File vcf = FinalizeVCF.gcs_path
        File tbi = FinalizeTBI.gcs_path
    }

    if (!(defined(list_of_snfs) || defined(snfs))) {
        call Utils.StopWorkflow as MissingInputs {input: reason = "Please provide as least one valid input"}
    }
    if (defined(list_of_snfs) && defined(snfs)) {
        call Utils.StopWorkflow as MutexInputs {input: reason = "Please provide only one of [list_of_snfs, snfs]"}
    }

    String workflow_name = 'Sniffles2JointCallSV'
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}"

    ##################
    if (defined(list_of_snfs) ) {
        call ReadLines { input: one_record_per_line = select_first([list_of_snfs])}
    }
    Array[File] cohort_snfs = select_first([ReadLines.records, snfs])

    call SNF.Joint as SNF2Join { input:
        snfs = cohort_snfs, prefix = out_prefix, cores = cores, disk_type = disk_type
    }
    ##################
    call FF.FinalizeToFile as FinalizeVCF { input: outdir = outdir, file = SNF2Join.vcf }
    call FF.FinalizeToFile as FinalizeTBI { input: outdir = outdir, file = SNF2Join.tbi }
}

task ReadLines {
    meta {
        desciption: "When a file is long enough, WDL stdlib read_lines can get stuck, so we do the work here ourselves"
    }
    input {
        File one_record_per_line
    }
    output {
        Array[String] records = read_lines("unnecessary.txt")
    }

    command <<<
        set -eux
        cat ~{one_record_per_line} > "unnecessary.txt"
    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
