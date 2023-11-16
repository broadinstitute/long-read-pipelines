version 1.0


import "../../../tasks/Kmers/Jellyfish.wdl" 
import "../../../tasks/Utility/Finalize.wdl" as FF

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
        String gcs_out_root_dir
    }
    
    call Jellyfish.KmerCounts as JF_kmercount { input: 
        fasta = fa,
        kmer_size = k,
        prefix = sampleid
        }
    
    call FF.FinalizeToFile as Finalizescaffold {
        input: outdir = gcs_out_root_dir, file = JF_kmercount.kmercount
    }
    
    output{
        File KmerCounts = JF_kmercount.kmercount
    }
}