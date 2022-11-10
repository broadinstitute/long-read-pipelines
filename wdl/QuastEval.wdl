version 1.0

import "tasks/Quast.wdl" as Quast

workflow QuastEval {
    input {
        File reference

        Array[File] assemblies
    }

    call Quast.Quast as Eval {
        input:
            ref = reference,
            assemblies = assemblies
    }

    output {
        File report_html = Eval.report_html
        File report_txt = Eval.report_txt
        File report_pdf = Eval.report_pdf
        Array[File] plots = Eval.plots

        Map[String, String] metrics = Eval.metrics
    }
}
