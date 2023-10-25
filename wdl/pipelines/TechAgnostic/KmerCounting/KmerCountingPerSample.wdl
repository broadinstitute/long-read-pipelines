version 1.0


import "../../../tasks/Kmers/Jellyfish.wdl" 


workflow KmerCountingPerSample {
    meta{
        description : "..."
    }
    parameter_meta {
    }

    input {
        File fa
        String sampleid
        Int k
    }
    
    call Jellyfish.KmerCounts as JF_kmercount { input: 
        fasta = fa,
        kmer_size = k,
        prefix = sampleid
        }
    
    output{
        File KmerCounts = JF_kmercount.kmercount
    }
}