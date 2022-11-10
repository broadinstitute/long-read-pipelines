version 1.0

import "tasks/Quast.wdl" as Quast

workflow QuastEval {
    input {
        File reference
        File? reference_gff
        Array[File] assemblies

        Boolean icarus = false
    }

    call Quast.Quast as Eval {
        input:
            ref = reference,
            ref_gff = reference_gff,
            icarus = icarus,

            assemblies = assemblies
    }

    output {
        File report_html = Eval.report_html
        File report_txt = Eval.report_txt
        File report_pdf = Eval.report_pdf
        Array[File] plots = Eval.plots
        Array[File] icarus_main = Eval.icarus_main
        Array[File] icarus_viewers = Eval.icarus_viewers

        Map[String, String] metrics = Eval.metrics
    }
}
