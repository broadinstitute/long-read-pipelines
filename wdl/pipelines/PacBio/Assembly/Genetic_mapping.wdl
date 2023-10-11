workflow genetic_mapping{
    meta{
        description: "a workflow that extract vcfs and bams in a genomic interval"
    }
    input{
        File genetic_mapping
    }
    
    output {
        Map[String, String] mapping = read_map(genetic_mapping)
    }
}

    