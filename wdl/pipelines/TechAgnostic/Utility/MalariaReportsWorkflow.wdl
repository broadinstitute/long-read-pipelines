version 1.0

# Import malaria reports/summary generation as MRS
import "../../../tasks/Visualization/MalariaReports.wdl" as MRS

workflow GenerateMalariaReports {

    meta {
        author: "Bridget Knight"
        description: "## Report Generation \n This workflow calls the Python script that generates a sequencing report."
    }

    parameter_meta {

    }

    input {
        
    }
    
    call MRS.RunReportScript { 
        input: 

    }
}