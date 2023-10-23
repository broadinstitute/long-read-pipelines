version 1.0


import "../../../tasks/Kmers/Jellyfish.wdl" 


workflow HybridPhase {
    meta{
        description : "..."
    }
    parameter_meta {
    }

    input {
        Array[File] fastas
        Array[File] sampleIDs
        File reference_fasta
        String prefix
        Int k
        Int num_t
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
        

    call Jellyfish.CountingKmerInFile as kmercounttable { input:
        jfs = JF_kmercount.kmercount
    }

    output{
        File final_kmercount = kmercounttable.final_counttable
    }
}