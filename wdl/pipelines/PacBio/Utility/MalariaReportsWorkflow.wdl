version 1.0

# Import malaria reports/summary generation as MRS
import "../../../tasks/Visualization/MalariaReports.wdl" as MRS

workflow GenerateMalariaReports {

    input {
        Int mem_gb
        String docker_image
    }
    
    call MRS.RunReportScript { 
        input: 
            mem_gb = mem_gb 
            docker_image = docker_image
    }

    meta {
        author: "Bridget Knight"
        description: "## Report Generation \n This workflow calls the Python script that generates a sequencing report."
    }

}