version 1.0

import "tasks/ConvertONTQCNotebook.wdl" as Convert


workflow ONTQCReport {
    input {
        String workspace_namespace
        String workspace_name

        String sample_name
        String notebook_path
    }

    call Convert.ConvertONTQCNotebook as Report {
        input:
            workspace_namespace=workspace_namespace,
            workspace_name=workspace_name,
            sample_name=sample_name,
            notebook_path=notebook_path
    }

    output {
        File ont_qc_report = Report.ont_qc_report
    }
}
