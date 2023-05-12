version 1.0

import "../../Illumina/Preprocessing/IllumFastP.wdl" as FastP
import "../../../tasks/Assembly/Unicycler.wdl" as Unicycler

workflow AssembleWithUnicycler {
    input {
        String sample_name

        File illumina_fq1
        File illumina_fq2
        File? long_reads
    }

    call FastP.IllumFastP as ProcessIllum {
        input:
            illumina_fq1=illumina_fq1,
            illumina_fq2=illumina_fq2
    }

    call Unicycler.Unicycler as Assemble {
        input:
            sample_name=sample_name,
            illumina_fq1=ProcessIllum.processed_fq1,
            illumina_fq2=ProcessIllum.processed_fq2,
            illumina_unpaired=ProcessIllum.unpaired_fq,
            long_reads=long_reads
    }

    output {
        File processed_fq1 = ProcessIllum.processed_fq1
        File processed_fq2 = ProcessIllum.processed_fq2
        File unpaired_fq = ProcessIllum.unpaired_fq
        File fastp_report = ProcessIllum.fastp_report
        File fastp_json = ProcessIllum.fastp_json

        File assembly_fasta = Assemble.assembly_fasta
        File assembly_gfa = Assemble.assembly_gfa
        File unicycler_log = Assemble.unicycler_log
        Array[File] unicycler_intermediate_graphs = Assemble.unicycler_intermediate_graphs
    }

}