version 1.0

import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/Utility/ReadLengths.wdl" as ReLU

workflow DystPeaker {
    meta {
        description: "Collect read length information from a long reads BAM."
    }
    input {
        File input_file
        Boolean input_is_bam
        String id
        Int short_reads_threshold

        String gcs_out_root_dir
    }
    parameter_meta {
        gcs_out_root_dir: "Cloud storage output directory"
        id: "A distinguishing ID that's going to impact how the files are named and where they are placed in the directories."
        short_reads_threshold: "A threshold below which the reads will be classified as short"

        read_lengths_hist: "Read length histogram"
        peaks: "Estimated peaks in the read length distritbution"
        reverse_yield: "A lenth-9 array of lengths at which a certain fraction of reads are shorter than. The fraction bins are 10% to 90% with 10% increments."
        read_len_summaries: "A summary on some other metrics related to read length"
    }

    String relative_dir = "ReadLengthMetrics"
    String output_dir = sub(gcs_out_root_dir, "/$", "") + "/" + relative_dir + "/" + id

    # collect
    if (input_is_bam) {
        call ReLU.GetLengthsFromBam   { input: bam   = input_file }
    }
    if ( !input_is_bam ) {
        call ReLU.GetLengthsFromFastq { input: fastq = input_file }
    }
    File rl_file = select_first([GetLengthsFromBam.read_lengths, GetLengthsFromFastq.read_lengths])

    # stats
    call ReLU.Dyst { input: read_lengths_txt = rl_file }
    call ReLU.Peaker { input: dyst_histogram = Dyst.histogram }
    call ReLU.ReverseYield { input: read_lengths_txt = rl_file }
    call ReLU.Skewness { input: read_lengths_txt = rl_file }
    call ReLU.GetNumReadsAndShorts { input: read_lengths_txt = rl_file, short_threshold = short_reads_threshold }
    String raw_pct = round(100 * GetNumReadsAndShorts.num_shorts/GetNumReadsAndShorts.num_seqs)

    call FF.FinalizeToFile as SaveRLArray { input: outdir = output_dir, file = GetNumReadsAndShorts.rl_bz2, name = id + ".readLen.txt.bz2" }
    call FF.FinalizeToFile as SaveHist    { input: outdir = output_dir, file = Dyst.histogram, name = id + ".readLen.hist.txt" }

    output {
        File read_len_hist = SaveHist.gcs_path
        Array[Int] read_len_peaks = Peaker.peaks
        Array[Int] read_len_deciles = ReverseYield.reverse_yield

        Map[String, String] read_len_summaries = {'shortie_pct': raw_pct + "%",
                                                  'shortie_threshold': short_reads_threshold,
                                                  'skew': Skewness.skew,
                                                  'raw_rl_file': SaveRLArray.gcs_path}
    }
}
