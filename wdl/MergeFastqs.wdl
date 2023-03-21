version 1.0

import "tasks/Utils.wdl"
import "tasks/Finalize.wdl" as FF

workflow MergeFastqs {
    input {
        Array[File] fastq_gzs
        String sample

        String gcs_out_root_dir
    }

    if (1<length(fastq_gzs)) {
        call Utils.MergeFastqGzs {input:
            fastq_gzs = fastq_gzs,
            prefix = sample
        }
    }
    File to_save = select_first([MergeFastqGzs.merged_fastq_gz, fastq_gzs[0]])
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/FASTQs/"
    call FF.FinalizeToFile {input:
        file = to_save, outdir = outdir, name = "~{sample}.fastq.gz"
    }

    output {
        File FASTQ = FinalizeToFile.gcs_path
    }
}
