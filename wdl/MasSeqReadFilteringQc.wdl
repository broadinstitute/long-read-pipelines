version 1.0

import "tasks/Structs.wdl"
import "tasks/Utils.wdl" as Utils

workflow MasSeqReadFilteringQc {

    input {
        File annotated_aligned_array_element_bam
        File whitelist
    }

    call Utils.CountBamRecords as t_01_CountInputReads {
        input:
            bam = annotated_aligned_array_element_bam
    }

    # 1 - Primary:
    call Utils.FilterReadsBySamFlags as t_02_GetPrimaryOnly {
        input:
            bam = annotated_aligned_array_element_bam,
            sam_flags = "2308",
            prefix = "primary_only",
    }
    call Utils.CountBamRecords as t_03_CountPrimaryOnly {
        input:
            bam = t_02_GetPrimaryOnly.output_bam
    }

    # 2 - MQ > 0:
    call PrintReads as t_04_GetNonZeroMqReads {
        input:
            bam_file = t_02_GetPrimaryOnly.output_bam,
            bam_index = t_02_GetPrimaryOnly.output_bam_index,
            prefix = "primary_only.non_zero_mq",
            filters = ""
    }
    call Utils.CountBamRecords as t_05_CountNonZeroMq {
        input:
            bam = t_04_GetNonZeroMqReads.bam
    }

    # 2 - MQ > 20:
    call SamtoolsView as t_06_GetMqGt20Reads {
        input:
            bam_file = t_02_GetPrimaryOnly.output_bam,
            args = "-q 21",
            prefix = "primary_only.non_zero_mq.qual_gt_20",
    }
    call Utils.CountBamRecords as t_07_CountMqGt20 {
        input:
            bam = t_06_GetMqGt20Reads.bam
    }

    # 3 - length < 15kbp
    call Utils.Bamtools as t_08_GetLenLt15kbpReads {
        input:
            bamfile = t_06_GetMqGt20Reads.bam,
            cmd = "filter",
            args = '-length <15000',
            prefix = "primary_only.non_zero_mq.qual_gt_20.length_lt_15kbp",
    }
    call Utils.CountBamRecords as t_09_CountLenLt15kbp {
        input:
            bam = t_08_GetLenLt15kbpReads.bam_out
    }

    # 4 - < 1000 soft clipped bases:
    call PrintReads as t_10_GetSoftClipLt1000Reads {
        input:
            bam_file = t_08_GetLenLt15kbpReads.bam_out,
            prefix = "primary_only.non_zero_mq.qual_gt_20.length_lt_15kbp.soft_clip_lt_1000",
            filters = " --read-filter ExcessiveEndClippedReadFilter --max-clipped-bases 1000 "
    }
    call Utils.CountBamRecords as t_11_CountSoftClipLt1000 {
        input:
            bam = t_10_GetSoftClipLt1000Reads.bam
    }

    # 5 - Reads in CB whitelist:
    call FilterReadsByCellBarcodeWhitelist as t_12_GetCbWhitelistReads {
        input:
            bam_file = t_08_GetLenLt15kbpReads.bam_out,
            whitelist = whitelist,
            prefix = "primary_only.non_zero_mq.qual_gt_20.length_lt_15kbp.soft_clip_lt_1000.cb_whitelist",
    }

    # 6: Call output:
    call CreateCountsTable as t_13_CreateFinalCountsTable {
        input:
            numReadsIn = t_01_CountInputReads.num_records,
            numPrimary = t_03_CountPrimaryOnly.num_records,
            numMqGt0 = t_05_CountNonZeroMq.num_records,
            numMqGt20 = t_07_CountMqGt20.num_records,
            numLenLt15kbp = t_09_CountLenLt15kbp.num_records,
            numClipLenLt1000 = t_11_CountSoftClipLt1000.num_records,
            numMatchWhitelist = t_12_GetCbWhitelistReads.num_reads_out,
    }

    output {
        File final_bam = t_12_GetCbWhitelistReads.bam
        File counts_table = t_13_CreateFinalCountsTable.counts
    }
}

task CreateCountsTable {
    input {
        Int numReadsIn
        Int numPrimary
        Int numMqGt0
        Int numMqGt20
        Int numLenLt15kbp
        Int numClipLenLt1000
        Int numMatchWhitelist

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    String out_file = "~{prefix}.tsv"

    command <<<
        echo -e "#_Reads_In\t#_Primary\t#_Mq_>_0\t#_Mq_>_20\t#_len_<_15kbp\t#_end_clips_<_1000bp\t#_match_cb_whitelist" >> ~{out_file}
        echo -e "~{numReadsIn}\t~{numPrimary}\t~{numMqGt0}\t~{numMqGt20}\t~{numLenLt15kbp}\t~{numClipLenLt1000}\t~{numMatchWhitelist}" >> ~{out_file}
    >>>

    output {
        File counts = "~{out_file}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:  1,
        docker:             "ubuntu:18.04"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task SamtoolsView {
    input {
        File bam_file
        String args

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 1 + 10*ceil(2 * size(bam_file, "GiB"))

    command <<<

        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        samtools view -@{np} ~{args} ~{bam_file} > ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task FilterReadsByCellBarcodeWhitelist {
    input {
        File bam_file
        File whitelist

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 1 + 10*ceil(2 * size(bam_file, "GiB")) + ceil(size(whitelist, "GiB"))

    String out_file_name = "~{prefix}.bam"

    command <<<
python << END
        import pysam
        from tqdm import tqdm

        # Get our barcode list:
        with open(~{whitelist}, 'r') as f:
            whitelist_set = set()
            for line in f:
                barcode = line.strip()
                whitelist_set.add(barcode)

        # Filter the reads by those in the whitelist
        num_barcodes_in_whitelist = 0
        num_reads_in_file = 0
        with pysam.AlignmentFile("~{bam_file}", "rb", check_sq=False, require_index=False) as bam_file, tqdm(desc="Filtering barcode whitelisted reads...", unit="read") as pbar:
            with pysam.AlignmentFile("~{out_file_name}", 'wb', check_sq=False, header=bam_file.header) as output_file:
                for r in bam_file:
                    barcode = r.get_tag("CR")
                    if barcode in whitelist_set:
                        num_barcodes_in_whitelist += 1
                        output_file.write(r)

                    num_reads_in_file += 1
                    pbar.update(1)

        print(f"# Reads in file: {num_reads_in_file}")
        print(f"# Reads with whitelisted barcodes in file: {num_barcodes_in_whitelist}")
        print(f"# Reads without whitelisted barcodes in file: {num_reads_in_file - num_barcodes_in_whitelist}")

        with open('num_reads_in.txt', 'w') as f:
            f.write(f"{num_reads_in_file}\n")

        with open('num_reads_out.txt', 'w') as f:
            f.write(f"{num_barcodes_in_whitelist}\n")
END
    >>>

    output {
        File bam = "~{out_file_name}"
        Int num_reads_in = read_int("num_reads_in.txt")
        Int num_reads_out = read_int("num_reads_out.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.10"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task PrintReads {
    input {
        File bam_file
        File? bam_index

        String filters

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 1 + 10*ceil(2 * size(bam_file, "GiB")) + ceil(size(bam_index, "GiB"))

    command {
        /gatk/gatk PrintReads \
            -I ~{bam_file} \
            -O ~{prefix}.bam \
            --disable-read-filter WellformedReadFilter \
            ~{filters}
    }

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             2,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "broadinstitute/gatk-nightly:caa48f98c8207b688db6ee35fead3eafb7219e38"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}



task FilterMasSeqReadsWithGatk {
    input {
        File bam_file
        File bam_index

        Int maxReadLength = 15000
        Int maxEndClipping = 1000

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 1 + 10*ceil(2 * size(bam_file, "GiB")) + ceil(size(bam_index, "GiB"))

    command {
        /gatk/gatk PrintReads \
            -I ~{bam_file} \
            -O ~{prefix}.bam \
            --disable-read-filter WellformedReadFilter \
            --read-filter MappedReadFilter \
            --read-filter MappingQualityNotZeroReadFilter \
            --read-filter NotSecondaryAlignmentReadFilter \
            --read-filter NotSupplementaryAlignmentReadFilter \
            --read-filter ReadLengthReadFilter --max-read-length 15000 \
            --read-filter ExcessiveEndClippedReadFilter --max-clipped-bases 1000

        echo "PWD is:"
        pwd

        echo "PWD List:"
        ls -la

        echo "Outfile list:"
        ls -la ~{prefix}.bam*

        date
    }

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             2,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "broadinstitute/gatk-nightly:caa48f98c8207b688db6ee35fead3eafb7219e38"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}