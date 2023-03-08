version 1.0


workflow ParseRevioFolder {
    meta {
        description:
        "A lightweight workflow to pick apart files of various types in a Revio SMRT cell folder"
    }

    input {
        String smrtcell_folder
    }

    parameter_meta {
        smrtcell_folder: "The (cloud) location of SMRT cell folder"
    }

    String formatted_folder_path = sub(smrtcell_folder, "/$", "")

    call GetRunId { input: smrtcell_folder = formatted_folder_path }

    # barcode to various file paths
    call MapBarcodeToFilePaths { input: smrtcell_folder = formatted_folder_path }

    # barcode to various sample ids
    Array[Pair[String, String]] coerced_pair_from_map = MapBarcodeToFilePaths.bc_2_xml
    scatter (pair in coerced_pair_from_map) {
        call MapBarcodeToSampleIDs { input: barcode = pair.left, xml_for_barcode = pair.right}
    }
    call CoercePairsOfArrayToMap as AssociateToBioSample { input: keys = MapBarcodeToSampleIDs.redundant_barcode, values = MapBarcodeToSampleIDs.biosample }
    call CoercePairsOfArrayToMap as AssociateToAliquot { input: keys = MapBarcodeToSampleIDs.redundant_barcode, values = MapBarcodeToSampleIDs.aliquot }
    call CoercePairsOfArrayToMap as AssociateToLIMS { input: keys = MapBarcodeToSampleIDs.redundant_barcode, values = MapBarcodeToSampleIDs.limsid }

    output {
        String run_id = GetRunId.run_id
        String cell_idex = GetRunId.cell_idx
        String movie = MapBarcodeToFilePaths.movie_name

        Map[String, String] bc_2_bam        = MapBarcodeToFilePaths.bc_2_bam

        Map[String, String] bc_2_aliquot_id   = AssociateToAliquot.res
        Map[String, String] bc_2_biosample_id = AssociateToBioSample.res
        Map[String, String] bc_2_lims_id      = AssociateToLIMS.res

        Map[String, String] bc_2_oninstrument_metrics = MapBarcodeToFilePaths.bc_2_oninstrument_metrics
    }
}

task GetRunId {
    meta {
        desciption: "The the run id of a particular SMRT cell"
    }
    input {
        String smrtcell_folder
    }

    command <<<
        set -eux
        # example path "gs://you_cannot_guess_where_I_belong/r1234567_123456_123456/1_A01"
        echo "~{smrtcell_folder}" | awk -F '/' '{print $(NF-1)}' > run_id.txt
        echo "~{smrtcell_folder}" | awk -F '/' '{print $NF}' > cell_index.txt
    >>>
    output {
        String run_id = read_string("run_id.txt")
        String cell_idx = read_string("cell_index.txt")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task MapBarcodeToFilePaths {
    meta {
        desciption:
        "Separate out barcode and their associated files and metadata."
    }
    input {
        String smrtcell_folder
    }

    parameter_meta {
        smrtcell_folder: "The (cloud) location of SMRT cell folder"
    }

    String metrics_tsv = "metrics.tsv"
    String barcode_2_hifibam = "bc_2_hifi.tsv"
    String barcode_2_failbam = "bc_2_fail.tsv"
    String barcode_2_xml     = "bc_2_xml.tsv"

    command <<<
        set -eux

        movie_name=$(gsutil ls ~{smrtcell_folder}/"metadata"/*.ccs.log | awk -F '/' '{print $NF}' | sed 's#.ccs.log$##')
        echo "${movie_name}" > movie_name.txt

        # to simply parse and extract barcode and various sample IDs
        gsutil ls ~{smrtcell_folder}/"pb_formats"/ | grep -E "${movie_name}.hifi_reads.bc[0-9]{4}.consensusreadset.xml" > xml.files.txt
        while IFS= read -r xml
        do
            barcode=$(echo "$xml" | grep -Eo "bc[0-9]{4}")
            echo -e "${barcode}\t${xml}" >> ~{barcode_2_xml}
        done < xml.files.txt

        # hifi pass BAM
        gsutil ls ~{smrtcell_folder}/"hifi_reads"/*.bam | grep -vF 'unassigned' > hifi_reads.txt
        while IFS= read -r bam
        do
            barcode=$(echo "$bam" | grep -Eo "bc[0-9]{4}")
            echo -e "${barcode}\t${bam}" >> ~{barcode_2_hifibam}
        done < hifi_reads.txt

        # hifi fail BAM
        gsutil ls ~{smrtcell_folder}/"fail_reads"/*bam | grep -vF 'unassigned' > fail_reads.txt
        while IFS= read -r bam
        do
            barcode=$(echo "$bam" | grep -Eo "bc[0-9]{4}")
            echo -e "${barcode}\t${bam}" >> ~{barcode_2_failbam}
        done < fail_reads.txt

        # on-instrument SMRT cell level metric files
        gsutil ls ~{smrtcell_folder}/"statistics"/ | grep -vF 'fail_reads' > metadata.txt
        loc=$(grep -F 'summary.json' metadata.txt)
        echo -e "smrtcell_summary_json\t${loc}" >> ~{metrics_tsv}
        loc=$(grep -F 'ccs_report.txt' metadata.txt)
        echo -e "ccs_report_txt\t${loc}" >> ~{metrics_tsv}
        loc=$(grep -F 'ccs_report.json' metadata.txt)
        echo -e "ccs_report_json\t${loc}" >> ~{metrics_tsv}
        loc=$(grep -F 'hifi_reads.5mc_report.json' metadata.txt)
        echo -e "5mc_report_json\t${loc}" >> ~{metrics_tsv}
        loc=$(grep -F 'hifi_reads.lima_counts.txt' metadata.txt)
        echo -e "lima_counts_txt\t${loc}" >> ~{metrics_tsv}
        loc=$(grep -F 'hifi_reads.lima_summary.txt' metadata.txt)
        echo -e "lima_summary_txt\t${loc}" >> ~{metrics_tsv}

        # for now, we omit exploring file contents in the metadata folder (holding log files)
        gsutil ls ~{smrtcell_folder}/"metadata"/ > statistics.txt
    >>>
    output {
        String movie_name = read_string("movie_name.txt")

        Map[String, String] bc_2_xml = read_map(barcode_2_xml)

        Map[String, String] bc_2_bam          = read_map(barcode_2_hifibam)
        Map[String, String] bc_2_fail_bam     = read_map(barcode_2_failbam)
        Map[String, String] bc_2_oninstrument_metrics = read_map(metrics_tsv)
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:alpine"
    }
}

task MapBarcodeToSampleIDs {
    meta {
        desciption: "Parse a PacBio XML for a Revio SMRT cell, and pick out the various forms of sample IDs. Note this is lab-dependent."
    }
    input {
        String barcode
        File xml_for_barcode
    }
    parameter_meta {
        barcode:         ""
        xml_for_barcode: ""
    }

    command <<<
        set -eux
        mv ~{xml_for_barcode} current.xml
        echo ~{barcode} > 'current.barcode.txt'

        python << CODE
        import xmltodict, json
        with open('current.barcode.txt') as inf:
            current_barcode=inf.read().strip()
        with open('current.xml') as inf:
            doc = xmltodict.parse(inf.read())

        sample_name_meta = doc['pbds:ConsensusReadSet']['pbds:DataSetMetadata']['Collections']['CollectionMetadata']['pbmeta:WellSample']

        # this section is convention-dependent, hence may change from lab to lab
        aliquot_id = sample_name_meta['@Description']
        lims_id    = sample_name_meta['@Name']
        biosample  = sample_name_meta['pbsample:BioSamples']['pbsample:BioSample']['@Name']

        with open('biosample.id', 'w') as outf:
            outf.write(f'{biosample}\n')
        with open('aliq.id', 'w') as outf:
            outf.write(f'{aliquot_id}\n')
        with open('lims.id', 'w') as outf:
            outf.write(f'{lims_id}\n')
        CODE
    >>>
    output {
        String redundant_barcode = barcode
        String biosample = read_string("biosample.id")
        String limsid    = read_string("lims.id")
        String aliquot   = read_string("aliq.id")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-smrtlinkparser:0.0.1"
    }
}

task CoercePairsOfArrayToMap {
    meta {
        desciption: "What can I say."
    }
    input {
        Array[String] keys
        Array[String] values
    }
    parameter_meta {
        keys:   "Array to be used as Map keys"
        values: "Array to be used as Map values"
    }

    command <<<
        set -eux

        paste ~{write_lines(keys)} ~{write_lines(values)} > "whatever.tsv"
    >>>
    output {
        Map[String, String] res = read_map("whatever.tsv")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}










task TaskName {
    meta {
        desciption: ""
    }
    input {

    }
    parameter_meta {

    }

    command <<<
        set -eux

    >>>
    output {

    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
