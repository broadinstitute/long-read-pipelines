version 1.0

import "../../structs/Structs.wdl"
# Run the Shasta assembler on a given set of data to produce a de-novo assembly.
#
# This WDL assumes that the `shasta` utility located in the following path of the given docker image:
#    /opt/shasta/build/shasta-install/bin/shasta
#
# Default values of assembly parameters are optimized for an assembly at coverage 60x.
# If your data have significantly different coverage, some changes in assembly parameters may be necessary to get good
# results.
#
# For an image that meets this criterion, see (on dockerhub):
#   jonnsmith/lrma_alignment_toolbox:latest
#
# For more information on the Shasta tool, see the official GitHub repository:
#   https://github.com/chanzuckerberg/shasta
#
# Description of inputs:
#
#   Required:
#     String docker                                               - Docker image in which to run this task.
#     File input_reads                                            - FASTA/FASTQ file containing reads from which to create an assembly.
#
#   Optional:
#     Int? reads_minReadLength                                    - Read length cutoff (=10000).
#     Int? reads_palindromicReads_maxSkip                         - Used for palindromic read detection. (=100)
#     Int? reads_palindromicReads_maxMarkerFrequency              - Used for palindromic read detection. (=10)
#     Float? reads_palindromicReads_alignedFractionThreshold      - Used for palindromic read detection. (=0.1)
#     Float? reads_palindromicReads_nearDiagonalFractionThreshold - Used for palindromic read detection. (=100)
#     Int? reads_palindromicReads_deltaThreshold                  - Used for palindromic read detection. (=100)
#     Int? kmers_k                                                - Length of marker k-mers (in run-length space). (=10)
#     Float? kmers_probability                                    - Probability that a k-mer is used as a marker. (0.1)
#     Boolean? kmers_suppressHighFrequencyMarkers                 - If set, high frequency k-mers are not used as markers. High frequency k-mers are those with enrichment greater than the value specified by Kmers.enrichmentThreshold.
#     Float? kmers_enrichmentThreshold                            - If Kmers.suppressHighFrequencyMarkers is set, this controls the enrichment threshold above which a k-mer is not considered as a possible marker. Enrichment is ratio of k-mer frequency in reads to random.  (=10.)
#     Int? minHash_m                                              - The number of consecutive markers that define a MinHash/LowHash feature. (=4)
#     Float? minHash_hashFraction                                 - Defines how low a hash has to be to be used with the LowHash algorithm. (=0.01)
#     Int? minHash_minHashIterationCount                          - The number of MinHash/LowHash iterations. (=10)
#     Int? minHash_maxBucketSize                                  - The maximum bucket size to be used by the MinHash/LowHash algorithm. (=10)
#     Int? minHash_minFrequency                                   - The minimum number of times a pair of reads must be found by the MinHash/LowHash algorithm in order to be considered a candidate alignment. (=2)
#     Int? align_maxSkip                                          - The maximum number of markers that an alignment is allowed to skip. (=30)
#     Int? align_maxTrim                                          - The maximum number of trim markers tolerated at the beginning and end of an alignment. (=30)
#     Int? align_maxMarkerFrequency                               - Marker frequency threshold. (=10)
#     Int? align_minAlignedMarkerCount                            - The minimum number of aligned markers for an alignment to be used. (=100)
#     Int? readGraph_maxAlignmentCount                            - The maximum alignments to be kept for each read. (=6)
#     Int? readGraph_minComponentSize                             - The minimum size (number of oriented reads) of a connected component of the read graph to be kept. (=100)
#     Int? readGraph_maxChimericReadDistance                      - Used for chimeric read detection. (=2)
#     Int? markerGraph_minCoverage                                - Minimum number of markers for a marker graph vertex. (=10)
#     Int? markerGraph_maxCoverage                                - Maximum number of markers for a marker graph vertex. (=100)
#     Int? markerGraph_lowCoverageThreshold                       - Used during approximate transitive reduction. (=0)
#     Int? markerGraph_highCoverageThreshold                      - Used during approximate transitive reduction. (=256)
#     Int? markerGraph_maxDistance                                - Used during approximate transitive reduction. (=30)
#     Int? markerGraph_edgeMarkerSkipThreshold                    - Used during approximate transitive reduction. (=100)
#     Int? markerGraph_pruneIterationCount                        - Number of prune iterations. (=6)
#     String? markerGraph_simplifyMaxLength                       - Maximum lengths (in markers) used at each iteration of simplifyMarkerGraph. (=10,100,1000)
#     Int? assembly_crossEdgeCoverageThreshold                    - Maximum average edge coverage for a cross edge of the assembly graph to be removed. (=0)
#     Int? assembly_markerGraphEdgeLengthThresholdForConsensus    - Controls assembly of long marker graph edges. (=1000)
#     String? assembly_consensusCaller                            - Selects the consensus caller for repeat counts. See the Shasta documentation for available choices. (=Bayesian:guppy-2.3.5-a)
#     Boolean? assembly_useMarginPhase                            - Used to turn on margin phase.
#     Boolean? assembly_storeCoverageData                         - Used to request storing coverage data.
#     Float? phasing_phasingSimilarityThreshold                   - The minimum phasing similarity for an edge to be added to the phasing graph. (=0.5)
#     Int? phasing_maxNeighborCount                               - The maximum number of phasing graph edges to be kept for each oriented read. (=6)
#
#   Runtime:
#     Int  mem                                                    - Amount of memory to give to the machine running each task in this workflow.
#     Int  preemptible_attempts                                   - Number of times to allow each task in this workflow to be preempted.
#     Int  disk_space_gb                                          - Amount of storage disk space (in Gb) to give to each machine running each task in this workflow.
#     Int  cpu                                                    - Number of CPU cores to give to each machine running each task in this workflow.
#     Int  boot_disk_size_gb                                      - Amount of boot disk space (in Gb) to give to each machine running each task in this workflow.
#
workflow Shasta {

    input {
        String docker = "us.gcr.io/broad-dsp-lrma/lr-shasta-marginpolish-helen:0.0.1"

        File input_reads

        Int? reads_minReadLength
        Int? reads_palindromicReads_maxSkip
        Int? reads_palindromicReads_maxMarkerFrequency
        Float? reads_palindromicReads_alignedFractionThreshold
        Float? reads_palindromicReads_nearDiagonalFractionThreshold
        Int? reads_palindromicReads_deltaThreshold
        Int? kmers_k
        Float? kmers_probability
        Boolean? kmers_suppressHighFrequencyMarkers
        Float? kmers_enrichmentThreshold
        Int? minHash_m
        Float? minHash_hashFraction
        Int? minHash_minHashIterationCount
        Int? minHash_maxBucketSize
        Int? minHash_minFrequency
        Int? align_maxSkip
        Int? align_maxTrim
        Int? align_maxMarkerFrequency
        Int? align_minAlignedMarkerCount
        Int? readGraph_maxAlignmentCount
        Int? readGraph_minComponentSize
        Int? readGraph_maxChimericReadDistance
        Int? markerGraph_minCoverage
        Int? markerGraph_maxCoverage
        Int? markerGraph_lowCoverageThreshold
        Int? markerGraph_highCoverageThreshold
        Int? markerGraph_maxDistance
        Int? markerGraph_edgeMarkerSkipThreshold
        Int? markerGraph_pruneIterationCount
        String? markerGraph_simplifyMaxLength
        Int? assembly_crossEdgeCoverageThreshold
        Int? assembly_markerGraphEdgeLengthThresholdForConsensus
        String? assembly_consensusCaller
        Boolean? assembly_useMarginPhase
        Boolean? assembly_storeCoverageData
        Float? phasing_phasingSimilarityThreshold
        Int? phasing_maxNeighborCount

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    call ShastaTask {
        input:
            docker_image                                         = docker,

            input_reads                                          = input_reads,

            reads_minReadLength                                  = reads_minReadLength,
            reads_palindromicReads_maxSkip                       = reads_palindromicReads_maxSkip,
            reads_palindromicReads_maxMarkerFrequency            = reads_palindromicReads_maxMarkerFrequency,
            reads_palindromicReads_alignedFractionThreshold      = reads_palindromicReads_alignedFractionThreshold,
            reads_palindromicReads_nearDiagonalFractionThreshold = reads_palindromicReads_nearDiagonalFractionThreshold,
            reads_palindromicReads_deltaThreshold                = reads_palindromicReads_deltaThreshold,
            kmers_k                                              = kmers_k,
            kmers_probability                                    = kmers_probability,
            kmers_suppressHighFrequencyMarkers                   = kmers_suppressHighFrequencyMarkers,
            kmers_enrichmentThreshold                            = kmers_enrichmentThreshold,
            minHash_m                                            = minHash_m,
            minHash_hashFraction                                 = minHash_hashFraction,
            minHash_minHashIterationCount                        = minHash_minHashIterationCount,
            minHash_maxBucketSize                                = minHash_maxBucketSize,
            minHash_minFrequency                                 = minHash_minFrequency,
            align_maxSkip                                        = align_maxSkip,
            align_maxTrim                                        = align_maxTrim,
            align_maxMarkerFrequency                             = align_maxMarkerFrequency,
            align_minAlignedMarkerCount                          = align_minAlignedMarkerCount,
            readGraph_maxAlignmentCount                          = readGraph_maxAlignmentCount,
            readGraph_minComponentSize                           = readGraph_minComponentSize,
            readGraph_maxChimericReadDistance                    = readGraph_maxChimericReadDistance,
            markerGraph_minCoverage                              = markerGraph_minCoverage,
            markerGraph_maxCoverage                              = markerGraph_maxCoverage,
            markerGraph_lowCoverageThreshold                     = markerGraph_lowCoverageThreshold,
            markerGraph_highCoverageThreshold                    = markerGraph_highCoverageThreshold,
            markerGraph_maxDistance                              = markerGraph_maxDistance,
            markerGraph_edgeMarkerSkipThreshold                  = markerGraph_edgeMarkerSkipThreshold,
            markerGraph_pruneIterationCount                      = markerGraph_pruneIterationCount,
            markerGraph_simplifyMaxLength                        = markerGraph_simplifyMaxLength,
            assembly_crossEdgeCoverageThreshold                  = assembly_crossEdgeCoverageThreshold,
            assembly_markerGraphEdgeLengthThresholdForConsensus  = assembly_markerGraphEdgeLengthThresholdForConsensus,
            assembly_consensusCaller                             = assembly_consensusCaller,
            assembly_useMarginPhase                              = assembly_useMarginPhase,
            assembly_storeCoverageData                           = assembly_storeCoverageData,
            phasing_phasingSimilarityThreshold                   = phasing_phasingSimilarityThreshold,
            phasing_maxNeighborCount                             = phasing_maxNeighborCount,

            mem_gb                                               = mem_gb,
            preemptible_attempts                                 = preemptible_attempts,
            disk_space_gb                                        = disk_space_gb,
            cpu                                                  = cpu,
            boot_disk_size_gb                                    = boot_disk_size_gb
    }

    output {
        File assembly_stranded_gfa              = ShastaTask.assembly_stranded_gfa
        File assembled_fasta                    = ShastaTask.assembled_fasta
        File assembly_gfa                       = ShastaTask.assembly_gfa
        File assembly_graph_dot                 = ShastaTask.assembly_graph_dot
        File assembly_chain_length_hist         = ShastaTask.assembly_chain_length_hist
        File assembly_summary_csv               = ShastaTask.assembly_summary_csv
        File assembly_summary_html              = ShastaTask.assembly_summary_html
        File assembly_summary_json              = ShastaTask.assembly_summary_json
        File binned_read_length_hist            = ShastaTask.binned_read_length_hist
        File marker_graph_edge_coverage_hist    = ShastaTask.marker_graph_edge_coverage_hist
        File marker_graph_vertex_coverage_hist  = ShastaTask.marker_graph_vertex_coverage_hist
        File palendromic_reads                  = ShastaTask.palendromic_reads
        File read_graph_components              = ShastaTask.read_graph_components
        File read_length_hist                   = ShastaTask.read_length_hist
        File read_summary                       = ShastaTask.read_summary
        File shasta_config                      = ShastaTask.shasta_config
        File timing_info                        = ShastaTask.timing_info
        File timing_info                        = ShastaTask.timing_info
    }
}

task ShastaTask {

    # ------------------------------------------------
    # Input args:
    input {
        # Docker:
        String docker_image

        # Required:
        File input_reads

        # Optional:
        Int? reads_minReadLength
        Int? reads_palindromicReads_maxSkip
        Int? reads_palindromicReads_maxMarkerFrequency
        Float? reads_palindromicReads_alignedFractionThreshold
        Float? reads_palindromicReads_nearDiagonalFractionThreshold
        Int? reads_palindromicReads_deltaThreshold
        Int? kmers_k
        Float? kmers_probability
        Boolean? kmers_suppressHighFrequencyMarkers
        Float? kmers_enrichmentThreshold
        Int? minHash_m
        Float? minHash_hashFraction
        Int? minHash_minHashIterationCount
        Int? minHash_maxBucketSize
        Int? minHash_minFrequency
        Int? align_maxSkip
        Int? align_maxTrim
        Int? align_maxMarkerFrequency
        Int? align_minAlignedMarkerCount
        Int? readGraph_maxAlignmentCount
        Int? readGraph_minComponentSize
        Int? readGraph_maxChimericReadDistance
        Int? markerGraph_minCoverage
        Int? markerGraph_maxCoverage
        Int? markerGraph_lowCoverageThreshold
        Int? markerGraph_highCoverageThreshold
        Int? markerGraph_maxDistance
        Int? markerGraph_edgeMarkerSkipThreshold
        Int? markerGraph_pruneIterationCount
        String? markerGraph_simplifyMaxLength
        Int? assembly_crossEdgeCoverageThreshold
        Int? assembly_markerGraphEdgeLengthThresholdForConsensus
        String? assembly_consensusCaller
        Boolean? assembly_useMarginPhase
        Boolean? assembly_storeCoverageData
        Float? phasing_phasingSimilarityThreshold
        Int? phasing_maxNeighborCount

        # Runtime Options:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    # ------------------------------------------------
    # Process input args:

    String bam_file_base_name = basename( input_reads )

    String kmers_suppressHighFrequencyMarkers_arg = if defined(kmers_suppressHighFrequencyMarkers) && kmers_suppressHighFrequencyMarkers then "--Kmers.suppressHighFrequencyMarkers" else ""
    String assembly_useMarginPhase_arg = if defined(assembly_useMarginPhase) && assembly_useMarginPhase then "--Assembly.useMarginPhase" else ""
    String assembly_storeCoverageData_arg = if defined(assembly_storeCoverageData) && assembly_storeCoverageData then "--Assembly.storeCoverageData" else ""

    String assembly_dir = "SHASTA_OUT"
    String out_assembly_stranded_gfa = assembly_dir + "/Assembly-BothStrands.gfa"
    String out_assembly_name = assembly_dir + "/Assembly.fasta"
    String out_assembly_gfa = assembly_dir + "/Assembly.gfa"
    String out_assembly_graph_dot = assembly_dir + "/AssemblyGraph-Final.dot"
    String out_assembly_chain_length_hist = assembly_dir + "/AssemblyGraphChainLengthHistogram.csv"
    String out_assembly_summary_csv = assembly_dir + "/AssemblySummary.csv"
    String out_assembly_summary_html = assembly_dir + "/AssemblySummary.html"
    String out_assembly_summary_json = assembly_dir + "/AssemblySummary.json"
    String out_binned_read_length_hist = assembly_dir + "/Binned-ReadLengthHistogram.csv"
    String out_marker_graph_edge_coverage_hist = assembly_dir + "/MarkerGraphEdgeCoverageHistogram.csv"
    String out_marker_graph_vertex_coverage_hist = assembly_dir + "/MarkerGraphVertexCoverageHistogram.csv"
    String out_palendromic_reads = assembly_dir + "/PalindromicReads.csv"
    String out_read_graph_components = assembly_dir + "/ReadGraphComponents.csv"
    String out_read_length_hist = assembly_dir + "/ReadLengthHistogram.csv"
    String out_read_summary = assembly_dir + "/ReadSummary.csv"
    String out_shasta_config = assembly_dir + "/shasta.conf"

    String timing_output_file = bam_file_base_name + ".timingInformation.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements:
    Int default_ram_mb = 3 * 1024

    Float reads_size_gb = size(input_reads, "GiB")
    Int default_disk_space_gb = ceil((reads_size_gb * 2) + 20)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    # ------------------------------------------------
    # Misc. helpers:
    String dollar = "$"

    # ------------------------------------------------
    # Run our command:
    command {

        # Find out how many threads to use:
        num_threads=~{cpu}
        found_count=false
        which lscpu &> /dev/null
        r=$?
        if [[ $r -eq 0 ]] ; then
            num_threads=$( lscpu | grep ^CPU\(s\) | sed 's#^CPU(s):[\t ]*##' )
            found_count=true
        fi

        if ! $found_count ; then
            if [[ -f /proc/cpuinfo ]] ; then
                num_threads=$( cat /proc/cpuinfo | grep ^processor | tail -n1 | sed 's#processor[ \t]*:[\t ]*##' )
            fi
        fi

        ############################

        # Now we can run marginPolish:
        set -e
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        sudo /opt/shasta/build/shasta-install/bin/shasta \
            --input ~{input_reads} \
            --assemblyDirectory ~{assembly_dir} \
            --threads $num_threads \
            ~{"--Reads.minReadLength" +  reads_minReadLength } \
            ~{"--Reads.palindromicReads.maxSkip" + reads_palindromicReads_maxSkip } \
            ~{"--Reads.palindromicReads.maxMarkerFrequency" + reads_palindromicReads_maxMarkerFrequency } \
            ~{"--Reads.palindromicReads.alignedFractionThreshold" + reads_palindromicReads_alignedFractionThreshold } \
            ~{"--Reads.palindromicReads.nearDiagonalFractionThreshold" + reads_palindromicReads_nearDiagonalFractionThreshold } \
            ~{"--Reads.palindromicReads.deltaThreshold" + reads_palindromicReads_deltaThreshold } \
            ~{"--Kmers.k" + kmers_k } \
            ~{"--Kmers.probability" + kmers_probability } \
            ~{kmers_suppressHighFrequencyMarkers_arg} \
            ~{"--Kmers.enrichmentThreshold" + kmers_enrichmentThreshold } \
            ~{"--MinHash.m" + minHash_m } \
            ~{"--MinHash.hashFraction" + minHash_hashFraction } \
            ~{"--MinHash.minHashIterationCount" + minHash_minHashIterationCount } \
            ~{"--MinHash.maxBucketSize" + minHash_maxBucketSize } \
            ~{"--MinHash.minFrequency" + minHash_minFrequency } \
            ~{"--Align.maxSkip" + align_maxSkip } \
            ~{"--Align.maxTrim" + align_maxTrim } \
            ~{"--Align.maxMarkerFrequency" + align_maxMarkerFrequency } \
            ~{"--Align.minAlignedMarkerCount" + align_minAlignedMarkerCount } \
            ~{"--ReadGraph.maxAlignmentCount" + readGraph_maxAlignmentCount } \
            ~{"--ReadGraph.minComponentSize" + readGraph_minComponentSize } \
            ~{"--ReadGraph.maxChimericReadDistance" + readGraph_maxChimericReadDistance } \
            ~{"--MarkerGraph.minCoverage" + markerGraph_minCoverage } \
            ~{"--MarkerGraph.maxCoverage" + markerGraph_maxCoverage } \
            ~{"--MarkerGraph.lowCoverageThreshold" + markerGraph_lowCoverageThreshold } \
            ~{"--MarkerGraph.highCoverageThreshold" + markerGraph_highCoverageThreshold } \
            ~{"--MarkerGraph.maxDistance" + markerGraph_maxDistance } \
            ~{"--MarkerGraph.edgeMarkerSkipThreshold" + markerGraph_edgeMarkerSkipThreshold } \
            ~{"--MarkerGraph.pruneIterationCount" + markerGraph_pruneIterationCount } \
            ~{"--MarkerGraph.simplifyMaxLength" + markerGraph_simplifyMaxLength } \
            ~{"--Assembly.crossEdgeCoverageThreshold" + assembly_crossEdgeCoverageThreshold } \
            ~{"--Assembly.markerGraphEdgeLengthThresholdForConsensus" + assembly_markerGraphEdgeLengthThresholdForConsensus } \
            ~{"--Assembly.consensusCaller" + assembly_consensusCaller } \
            ~{assembly_useMarginPhase_arg} \
            ~{assembly_storeCoverageData_arg} \
            ~{"--Phasing.phasingSimilarityThreshold" + phasing_phasingSimilarityThreshold } \
            ~{"--Phasing.maxNeighborCount" + phasing_maxNeighborCount }

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}
        elapsedTime=`python -c "print( $endTime - $startTime )"`
        echo "Elapsed Time: $elapsedTime" >> ~{timing_output_file}
    }

    # ------------------------------------------------
    # Runtime settings:
     runtime {
         docker: docker_image
         memory: machine_mem + " MB"
         disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
         bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
         preemptible: select_first([preemptible_attempts, 0])
         cpu: select_first([cpu, 1])
     }

    # ------------------------------------------------
    # Outputs:
    output {
        File assembly_stranded_gfa              = "${out_assembly_stranded_gfa}"
        File assembled_fasta                    = "${out_assembly_name}"
        File assembly_gfa                       = "${out_assembly_gfa}"
        File assembly_graph_dot                 = "${out_assembly_graph_dot}"
        File assembly_chain_length_hist         = "${out_assembly_chain_length_hist}"
        File assembly_summary_csv               = "${out_assembly_summary_csv}"
        File assembly_summary_html              = "${out_assembly_summary_html}"
        File assembly_summary_json              = "${out_assembly_summary_json}"
        File binned_read_length_hist            = "${out_binned_read_length_hist}"
        File marker_graph_edge_coverage_hist    = "${out_marker_graph_edge_coverage_hist}"
        File marker_graph_vertex_coverage_hist  = "${out_marker_graph_vertex_coverage_hist}"
        File palendromic_reads                  = "${out_palendromic_reads}"
        File read_graph_components              = "${out_read_graph_components}"
        File read_length_hist                   = "${out_read_length_hist}"
        File read_summary                       = "${out_read_summary}"
        File shasta_config                      = "${out_shasta_config}"
        File timing_info                        = "${timing_output_file}"
    }
 }
