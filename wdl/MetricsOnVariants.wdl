version 1.0

import "Structs.wdl"

workflow Metrics {
    input {
        File input_bam
        File input_bai
    }


}

# depth of coverage (autosome)
# depth of coverage (X)
# depth of coverage (Y)
# depth of coverage (MT)
# read lengths
# alignment stats
# yield
# number of passes
# error rate
# error rate / number of passes
# indel error lengths
# contamination
