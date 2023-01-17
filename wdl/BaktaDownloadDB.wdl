version 1.0

import "tasks/bacteria/Bakta.wdl" as Bakta

workflow BaktaDownloadDB {
    call Bakta.BaktaDBDownload as Download { }

    output {
        File bakta_db = Download.bakta_db
    }
}
