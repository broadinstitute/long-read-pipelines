version 1.0

task VerifyBamID {
    meta {
        desciption: "Uses VerifyBamID2 for human cross-individual contamination estimation. Assumes GRCh38."
    }

    input {
        File pileup
        File ref_fasta
        Boolean is_hgdp_sites
        Boolean is_100k_sites
    }

    String a = if is_hgdp_sites then 'hgdp' else '1000g.phase3'
    String b = if is_100k_sites then '100k' else  '10k'
    String resource_prefix = '~{a}.~{b}.b38.vcf.gz.dat'

    command <<<
        set -eux

        export VERIFY_BAM_ID_HOME='/VerifyBamID'

        time \
        ${VERIFY_BAM_ID_HOME}/bin/VerifyBamID \
            --SVDPrefix ${VERIFY_BAM_ID_HOME}/resource/~{resource_prefix} \
            --Reference ~{ref_fasta} \
            --PileupFile ~{pileup} \
            --NumThread 4 \
        > vbid2.out \
        2> vbid2.log

        cat vbid2.out
        tail -1 vbid2.out | awk -F ':' '{print $2}' | awk '{$1=$1};1' > "est_contam.txt"
    >>>

    output {
        File vbid2_log = "vbid2.log"
        File vbid2_out = "vbid2.out"
        Float contamination_est = read_float("est_contam.txt")
    }

    Int disk_size = 10 + ceil(size(pileup, "GiB"))
    runtime {
        cpu: 4
        memory: "8 GiB"
        disks: "local-disk ~{disk_size} SSD"
        docker: "us.gcr.io/broad-dsp-lrma/verifybamid2:v2.0.1"
    }
}