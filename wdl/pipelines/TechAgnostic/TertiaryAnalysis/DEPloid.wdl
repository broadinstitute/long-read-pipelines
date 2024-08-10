version 1.0

import "../../../structs/Structs.wdl"

workflow DEPloid {

    meta {
        description: "Run DEPloid on a VCF file to determine complexity of infection and perform variant deconvolution."
    }
    parameter_meta {
        sample_name: "Sample name"
        vcf: "VCF file to process"
        allele_freq_tsv: "TSV file with allele frequencies for each alt allele at each locus."
        ref_map_file:  "Reference map file indicating reference sequence and auxillary file locations"

        ibd_mode: "Run DEploid in the newer IBD mode (slower, but potentially more accurate)."
    }

    input {
        String sample_name
        File vcf
        File allele_freq_tsv
        File ref_map_file

        Boolean ibd_mode = false
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    # 1 - normalize VCF file
    call NormalizeVcfForAlleleSubsetting as t01_NormalizeVcfForAlleleSubsetting {
        input:
            prefix = sample_name,
            vcf_in = vcf,
            ref_fasta = ref_map['fasta']
    }

    # 2 - create subsetted PLAF and VCF
    call CreateSiteAltVcfSubset as t02_CreateSiteAltVcfSubset {
        input:
            prefix = sample_name,
            vcf_in = t01_NormalizeVcfForAlleleSubsetting.normalized_vcf,
            allele_freq_tsv = allele_freq_tsv
    }

    # 3 - run DEPloid
    call RunDePloid as t03_RunDePloid {
        input:
            prefix = sample_name,
            vcf_in = t02_CreateSiteAltVcfSubset.filtered_vcf,
            plaf = t02_CreateSiteAltVcfSubset.filtered_plaf,
            ibd_mode = ibd_mode
    }

    output {
        Float effective_coi = t03_RunDePloid.effective_coi
        Int inferred_coi = t03_RunDePloid.inferred_coi
        Float adjusted_effective_coi = t03_RunDePloid.adjusted_effective_coi

        Array[String] proportions = t03_RunDePloid.proportions

        File log = t03_RunDePloid.log
        File trace_log = t03_RunDePloid.trace_log
        File prop = t03_RunDePloid.prop
        File llk = t03_RunDePloid.llk
        File hap = t03_RunDePloid.hap
        File vcf_haps = t03_RunDePloid.vcf
        File ibd_probs = t03_RunDePloid.ibd_probs
        File final_hap = t03_RunDePloid.final_hap
    }
}

task NormalizeVcfForAlleleSubsetting {
    input {
        String prefix
        File vcf_in
        File ref_fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 15*ceil(size(vcf_in, "GB")) + 2*ceil(size(ref_fasta, "GB"))

    command <<<
        input_file=~{vcf_in}

        echo "${input_file}" | grep -q '.gz$'
        r=$?
        if [ $r -eq 0 ] ; then
          gunzip ${input_file}
          input_file=$(echo "${1}" | sed 's@.gz$@@')
        fi

        set -euxo pipefail

        REF=~{ref_fasta}
        bn=~{prefix}

        # First fudge VQSLOD because deploid needs it.
        # We can swap in the `SCORE` field from VQSR.
        # This is not great, but it will work...
        chrom_line_num=$(grep -n '^#CHROM' ${input_file} | sed 's@:.*@@')
        sed -e "${chrom_line_num}i##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=\"FUDGED FIELD.  THIS IS THE SAME AS SCORE.  ORIGINAL DESCRIPTION: Log odds of being a true variant versus being false under the trained gaussian mixture model\">" \
          -e 's@;SCORE=\([-0-9\.]*\);@;SCORE=\1;VQSLOD=\1;@' ${input_file} > ${bn}.VQSLOD_FUDGE.vcf

        # Now remove potentially invalid fields and normalize the vcf.
        # Normalized to one alt allele per line
        bcftools annotate -x"INFO/HAPCOMP,INFO/HAPDOM,INFO/HEC" ${bn}.VQSLOD_FUDGE.vcf | bcftools view -f "PASS,." | bcftools norm -m -any --atom-overlaps . -f ${REF} -Ov -o ${bn}.VQSLOD_FUDGE.norm.vcf

        # Cleanup the intermediate file:
        rm ${bn}.VQSLOD_FUDGE.vcf
    >>>

    output {
        File normalized_vcf = "${prefix}.VQSLOD_FUDGE.norm.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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

task CreateSiteAltVcfSubset {
    input {
        String prefix
        File vcf_in
        File allele_freq_tsv

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 5*ceil(size(vcf_in, "GB")) + ceil(1.1 * size(allele_freq_tsv, "GB"))

    command <<<
        python3 <<CODE

        ################################################################################

        import sys
        import os
        import re

        from collections import defaultdict
        from tqdm import tqdm

        ################################################################################

        input_vcf = "~{vcf_in}"
        out_vcf   = "~{prefix}.sites_filtered.vcf"
        out_plaf  = "~{prefix}.sites_filtered.plaf.tsv"

        ################################################################################

        # First slurp the TSV into memory for later ease of use:
        def get_num_lines(file_path):
            def get_blocks(files, size=65536):
                while True:
                    b = files.read(size)
                    if not b: break
                    yield b
            with open(file_path, "r",encoding="utf-8",errors='ignore') as f:
                return sum(bl.count("\n") for bl in get_blocks(f))

        pf7_af_sites = defaultdict(dict)
        num_lines = get_num_lines(pf7_af_sites_tsv)

        with open(pf7_af_sites_tsv, 'r') as f:
            header = next(f)
            for line in tqdm(f, desc="Reading PF7 sites TSV", total=num_lines-1):
                chrom, pos, ref, alt, af = line.strip().split()
                pos = int(pos)
                af = float(af)
                pf7_af_sites[(chrom, pos, ref)][alt] = af

        # Now go through the VCF and generate the filtered VCF and filtered sites TSV:

        def get_pos_ref_alt_fast(vcf_line):
            #CHROM  POS ID  REF ALT
            i1 = vcf_line.find('\t')
            chrom = vcf_line[:i1]

            i2 = vcf_line.find('\t', i1+1)
            pos = int(vcf_line[i1+1:i2])

            # Skip ID:
            i3 = vcf_line.find('\t', i2+1)

            i4 = vcf_line.find('\t', i3+1)
            ref = vcf_line[i3+1:i4]

            i5 = vcf_line.find('\t', i4+1)
            alt = vcf_line[i4+1:i5]

            return chrom, pos, ref, alt

        # We need to accumulate lines for each position until we have gotten all of them
        # Then we can output the line that matches the TSV data:
        pos_lines = []

        last_chrom = None
        last_pos = None
        last_ref = None
        last_alt = None

        num_lines = get_num_lines(input_vcf)

        num_variants_emitted = 0
        tot_variants = 0
        tot_sites = 0
        sites_matched = 0

        with open(input_vcf, 'r') as f_vcf_in:
            with open(out_vcf, 'w') as f_vcf_out, open(out_plaf, 'w') as f_plaf_out:

                f_plaf_out.write(f"CHROM\tPOS\tPLAF\n")

                for line in tqdm(f_vcf_in, desc="Filtering VCF by site alleles", total=num_lines):
                    if line.startswith("#"):
                        f_vcf_out.write(line)
                    else:
                        # we just care about the pos / ref / alt here, so do a quick parse:
                        chrom, pos, ref, alt = get_pos_ref_alt_fast(line)

                        if last_chrom is None:
                            last_chrom = chrom
                            last_pos = pos
                            last_ref = ref
                            last_alt = alt
                            tot_sites += 1

                        elif chrom != last_chrom or pos != last_pos:
                            # Time to write out a VCF line and the TSV line:

                            # Only emit if we have AF data for this site:
                            if (last_chrom, last_pos, last_ref) in pf7_af_sites:
                                # Get the alts that we could emit:
                                variant_info_in_tsv = [d for d in pos_lines if d[3] in pf7_af_sites[(last_chrom, last_pos, last_ref)]]

                                highest_AF = sys.float_info.min
                                line_to_emit = None
                                alt_to_emit = None

                                for d in variant_info_in_tsv:
                                    v_alt = d[3]
                                    # get the info for this particular alt:
                                    if pf7_af_sites[(last_chrom, last_pos, last_ref)][v_alt] > highest_AF:
                                        line_to_emit = d[4]
                                        alt_to_emit = v_alt

                                if line_to_emit is not None:
                                    # now we have our line to emit, so write it out
                                    # and write out the corresponding TSV line:
                                    f_vcf_out.write(f"{line_to_emit}")
                                    f_plaf_out.write(f"{last_chrom}\t{last_pos}\t{pf7_af_sites[(last_chrom, last_pos, last_ref)][alt_to_emit]}\n")

                                    num_variants_emitted += 1

                                sites_matched += 1

                            tot_sites += 1

                            # Update to our next info:
                            last_chrom = chrom
                            last_pos = pos
                            last_ref = ref
                            last_alt = alt

                            pos_lines = []

                        pos_lines.append((chrom, pos, ref, alt, line))
                        tot_variants += 1

                # Now emit the last site:
                if (last_chrom, last_pos, last_ref) in pf7_af_sites:
                    # Get the alts that we could emit:
                    variant_info_in_tsv = [d for d in pos_lines if d[3] in pf7_af_sites[(last_chrom, last_pos, last_ref)]]

                    highest_AF = sys.float_info.min
                    line_to_emit = None
                    alt_to_emit = None

                    for d in variant_info_in_tsv:
                        v_alt = d[3]
                        # get the info for this particular alt:
                        if pf7_af_sites[(last_chrom, last_pos, last_ref)][v_alt] > highest_AF:
                            line_to_emit = d[4]
                            alt_to_emit = v_alt

                    if line_to_emit is not None:
                        # now we have our line to emit, so write it out
                        # and write out the corresponding TSV line:
                        f_vcf_out.write(f"{line_to_emit}")
                        f_plaf_out.write(f"{last_chrom}\t{last_pos}\t{pf7_af_sites[(last_chrom, last_pos, last_ref)][alt_to_emit]}\n")

                        num_variants_emitted += 1

                    sites_matched += 1

        print(f"Variants Emitted: {num_variants_emitted}/{tot_variants}\t({100*num_variants_emitted/tot_variants:4.2f}%)")
        print(f"Sites Matched: {sites_matched}/{tot_sites}\t({100*sites_matched/tot_sites:4.2f}%)")

        CODE

    >>>

    output {
        File filtered_vcf = "${prefix}.sites_filtered.vcf"
        File filtered_plaf = "${prefix}.sites_filtered.plaf.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-malaria:0.0.1"
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

task RunDePloid {
    input {
        String prefix

        File vcf_in
        File plaf

        Boolean ibd_mode = false

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 5*ceil(size(vcf_in, "GB")) + ceil(1.1 * size(plaf, "GB"))

    String ibd_mode_flag = if ibd_mode then "-ibd" else ""

    command <<<
        time /DEploid/dEploid \
          -plaf ~{plaf} \
          -noPanel \
          -o ~{prefix} \
          -seed 1 \
          -nSample 800 \
          -rate 5 \
          -burn 0.5 \
          ~{ibd_mode_flag} \
          -vcfOut \
          -vcf ~{vcf_in}

        grep Effective_K ~{prefix}.log | awk '{print $2}' > effective_k.txt
        grep Inferred_K ~{prefix}.log | awk '{print $2}' > inferred_k.txt
        grep Adjusted_effective_K ~{prefix}.log | awk '{print $2}' > adjusted_effective_k.txt
        grep '^Proportions:' -A 1 ~{prefix}.log | tail -n 1 | tr '\t' '\n' | tr -d ' ' > proportions.txt
    >>>

    output {
        Float effective_coi = read_float("effective_k.txt")
        Int inferred_coi = read_int("inferred_k.txt")
        Float adjusted_effective_coi = read_float("adjusted_effective_k.txt")
        Array[String] proportions = read_lines("proportions.txt")

        File log = "${prefix}.log"
        File trace_log = "${prefix}.trace.log"
        File prop = "${prefix}.classic.prop"
        File llk = "${prefix}.classic.llk"
        File hap = "${prefix}.classic.hap"
        File vcf = "${prefix}.classic.vcf"
        File ibd_probs = "${prefix}.ibd.probs"
        File final_hap = "${prefix}.final.hap"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/util-malaria-coi:0.0.1"
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