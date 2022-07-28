#!/usr/bin/env awk -f

BEGIN {
    FS="\t"
    IFS="\t"
    OFS="\t"
    print "CCS Read Name\tBarcode\tConsensus Sequence\tICD1\tICD2\tICD3\tMapping Quality ICD1\tMapping Quality ICD2\tMapping Quality ICD3"
}
NR > 1 {
    print $1, $29, $48, $30, $36, $42, $35, $41, $47
}
END {

}