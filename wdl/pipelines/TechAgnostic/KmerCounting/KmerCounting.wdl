version 1.0


import "../../../tasks/Kmers/Jellyfish.wdl" 


workflow KmerCounting {
    meta{
        description : "..."
    }
    parameter_meta {
    }

    input {
        Array[File] fastas
        Array[String] sampleIDs
        File reference_fasta
        Int k
    }
    
    Int data_length = length(fastas)
    Array[Int] indexes= range(data_length)

    scatter (idx in indexes)  {
        File fa = fastas[idx]
        String sampleid = sampleIDs[idx]

        call Jellyfish.KmerCounts as JF_kmercount { input: 
            fasta = fa,
            kmer_size = k,
            prefix = sampleid
            }
    }
        
        call Jellyfish.KmerCounts as JF_kmercount_ref { input: 

            fasta = reference_fasta,
            kmer_size = k,
            prefix = "NC_012920.1"
            }
    call Jellyfish.CountingKmerInFile as kmercounttable { input:
        jfs = JF_kmercount.kmercount,
        ref_jf = JF_kmercount_ref.kmercount
    }

    output{
        File final_kmercount = kmercounttable.final_counttable
    }
}