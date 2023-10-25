version 1.0


import "../../../tasks/Kmers/Jellyfish.wdl" 


workflow KmerCounting {
    meta{
        description : "..."
    }
    parameter_meta {
    }

    input {
        Array[File] jfs
    }

    call Jellyfish.CountingKmerInFile as kmercounttable { input:
        jfs = jfs
    }

    output{
        File final_kmercount = kmercounttable.final_counttable
    }
}