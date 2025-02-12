version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# that are specifically designed to work with long read data.
# Each individual task/algo. is directly callable, if so desired.
##########################################################################################

import "../../tasks/Utility/Utils.wdl"
import "../../tasks/Utility/VariantUtils.wdl"

import "../../tasks/VariantCalling/PBSV.wdl"
import "Sniffles.wdl"


workflow CallVariants {
    meta {
        descrition: "A workflow for calling small and/or structural variants from an aligned CLR BAM file."
    }

    parameter_meta {
        bam:               "input BAM from which to call SVs"
        bai:               "index accompanying the BAM"

        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        ref_dict:          "sequence dictionary accompanying the reference"

        call_small_variants: "if true, will attempt to call small variants"
        call_small_vars_on_mitochondria: "if false, will not attempt to call variants on mitochondria"
        fast_less_sensitive_sv: "if true, will run PBSV in a less sensitive mode, which is faster"
        sites_vcf:     "for use with Clair, the small variant caller"
        sites_vcf_tbi: "for use with Clair, the small variant caller"

        prefix:            "prefix for output files"

        tandem_repeat_bed: "BED file containing TRF finder (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
    }

    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix

        Boolean call_svs
        Boolean fast_less_sensitive_sv
        File? tandem_repeat_bed

        Boolean call_small_variants = false
        Boolean call_small_vars_on_mitochondria
    }

    call Utils.RandomZoneSpewer as arbitrary {input: num_of_zones = 3}

    ######################################################################
    # Block for small variants handling
    ######################################################################
    # todo: use NanoCaller, Clair isn't going to support CLR
    ######################################################################
    # Block for SV handling
    ######################################################################
    if (call_svs) {
        if (fast_less_sensitive_sv) {

            call Utils.MakeChrIntervalList {
            input:
                ref_dict = ref_dict,
                filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']
            }

            scatter (c in MakeChrIntervalList.chrs) {
                String contig_for_sv = c[0]

                call Utils.SubsetBam {
                    input:
                        bam = bam,
                        bai = bai,
                        locus = contig_for_sv
                }

                call PBSV.RunPBSV {
                    input:
                        bam = SubsetBam.subset_bam,
                        bai = SubsetBam.subset_bai,
                        ref_fasta = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,
                        prefix = prefix,
                        tandem_repeat_bed = tandem_repeat_bed,
                        is_ccs = false,
                        zones = arbitrary.zones
                }

                call Sniffles.Sniffles {
                    input:
                        bam    = SubsetBam.subset_bam,
                        bai    = SubsetBam.subset_bai,
                        chr    = contig_for_sv,
                        prefix = prefix
                }

                call Utils.InferSampleName {
                    input:
                        bam = SubsetBam.subset_bam,
                        bai = SubsetBam.subset_bai
                }
                call FixSnifflesVCF {
                    input:
                        vcf = Sniffles.vcf,
                        sample_name = InferSampleName.sample_name
                }
            }

            call VariantUtils.MergePerChrCalls as MergePBSVVCFs {
                input:
                    vcfs     = RunPBSV.vcf,
                    ref_dict = ref_dict,
                    prefix   = prefix + ".pbsv"
            }

            call CollectDefinitions as UnionHeadersSnifflesVCFs {
                input:
                    vcfs = FixSnifflesVCF.sortedVCF
            }
            call VariantUtils.MergeAndSortVCFs as MergeSnifflesVCFs {
                input:
                    vcfs   = FixSnifflesVCF.sortedVCF,
                    ref_fasta_fai = ref_fasta_fai,
                    header_definitions_file = UnionHeadersSnifflesVCFs.union_definitions,
                    prefix = prefix + ".sniffles"
            }
        }

        if (!fast_less_sensitive_sv) {

            call PBSV.RunPBSV as PBSVslow {
                input:
                    bam = bam,
                    bai = bai,
                    ref_fasta = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    prefix = prefix,
                    tandem_repeat_bed = tandem_repeat_bed,
                    is_ccs = false,
                    zones = arbitrary.zones
            }
            call VariantUtils.ZipAndIndexVCF as ZipAndIndexPBSV {input: vcf = PBSVslow.vcf }

            call Sniffles.Sniffles as SnifflesSlow {
                input:
                    bam    = bam,
                    bai    = bai,
                    prefix = prefix
            }
            call Utils.InferSampleName as infer {input: bam = bam, bai = bai}
            call FixSnifflesVCF as ZipAndIndexSniffles {input: vcf = SnifflesSlow.vcf, sample_name = infer.sample_name}
        }
    }

    output {
        File? sniffles_vcf = select_first([MergeSnifflesVCFs.vcf, ZipAndIndexSniffles.sortedVCF])
        File? sniffles_tbi = select_first([MergeSnifflesVCFs.tbi, ZipAndIndexSniffles.tbi])

        File? pbsv_vcf = select_first([MergePBSVVCFs.vcf, ZipAndIndexPBSV.vcfgz])
        File? pbsv_tbi = select_first([MergePBSVVCFs.tbi, ZipAndIndexPBSV.tbi])
    }
}

task CollectDefinitions {

    meta {
        description: "Collect (union) various definitions in vcf files, adddressing a bcftols bug: https://github.com/samtools/bcftools/issues/1629"
    }

    parameter_meta {
        vcfs: "List of VCF files to be concatenated"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        Array[File] vcfs

        RuntimeAttr? runtime_attr_override
    }

    Int sz = ceil(size(vcfs, 'GB'))
    Int disk_sz = if sz > 100 then 5 * sz else 375

    command <<<
        set -euxo pipefail

        zgrep "^##" ~{vcfs[0]} > header.txt
        grep -F '##fileformat' header.txt > tmp.0.txt
        grep -F '##fileDate=' header.txt > tmp.1.txt
        if grep -q -F '##source=' header.txt; then grep -F 'source=' header.txt > tmp.2.txt; fi
        touch tmp.2.txt
        grep -F '##contig=' header.txt > tmp.3.txt

        cat tmp*txt > easy.txt && rm tmp*txt

        touch tmp.alt.txt tmp.ft.txt tmp.info.txt tmp.format.txt
        for vcf in ~{sep=' ' vcfs}; do
            zgrep -F '##ALT=' "${vcf}" >> tmp.alt.txt
            zgrep -F '##FILTER=' "${vcf}" >> tmp.ft.txt
            zgrep -F '##INFO=' "${vcf}" >> tmp.info.txt
            zgrep -F '##FORMAT=' "${vcf}" >> tmp.format.txt
        done
        for txt in tmp*txt; do
            sort "${txt}" | uniq > "${txt}.union"
        done
        cat tmp.alt.txt.union tmp.ft.txt.union tmp.info.txt.union tmp.format.txt.union > hard.txt
        cat easy.txt hard.txt > definitions.union.txt
    >>>

    output {
        File union_definitions = "definitions.union.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_sz,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FixSnifflesVCF {

    meta {
        description: "Fixes the sample information in a VCF file and prepares to fix undefined VCF INFO/FT/FORMATs. It then proceeds to get the missing VCF headers for these undefined formats and filters. Specific to Sniffles-1"
    }

    parameter_meta {
        sample_name:    "Sniffles infers sample name from the BAM file name, so we fix it here"
        ref_fasta_fai:  "provide only when the contig section of the input vcf is suspected to be corrupted"
    }

    input {
        File vcf
        String sample_name
        File? ref_fasta_fai
        RuntimeAttr? runtime_attr_override
    }

    Boolean fix_contigs = defined(ref_fasta_fai)

    Boolean vcf_is_bgzipped = sub(vcf, ".gz", "") != sub(vcf, ".vcf.gz", "")
    String local_raw = if vcf_is_bgzipped then "to.be.fixed.vcf.gz" else "to.be.fixed.vcf"
    String local_sp_fixed = if vcf_is_bgzipped then "sample.fixed.vcf.gz" else "sample.fixed.vcf"

    String initial_grep_cmd = if vcf_is_bgzipped then "zgrep" else "grep"

    String prefix = if vcf_is_bgzipped then basename(vcf, ".vcf.gz") else basename(vcf, ".vcf")
    Int proposed_disk = 3*ceil(size(vcf, "GB")) + 1
    Int disk_size = if (proposed_disk > 100) then proposed_disk else 100

    command <<<
        set -euxo pipefail

        # 1. fix sample information (Sniffles derives VCF SM information from the path to the BAM ......)
        cp ~{vcf} ~{local_raw}
        echo ~{sample_name} > sample_names.txt
        bcftools reheader --samples sample_names.txt -o ~{local_sp_fixed} ~{local_raw}
        rm ~{vcf} && rm ~{local_raw}

        ####################################################################
        # 2. prep for fixing undefined VCF INFO/FT/FORMAT, also guard against when the VCF is empty
        ~{initial_grep_cmd} "^##" ~{local_sp_fixed} > header.txt
        ~{initial_grep_cmd} -v "^#" ~{local_sp_fixed} > body.txt || true
        if [[ ! -f body.txt ]] || [[ ! -s body.txt ]]; then
            echo "input VCF seem to contain only header, but I'll proceed anyway and give you only header"
            bcftools \
                sort \
                --temp-dir tm_sort \
                --output-type z \
                -o ~{prefix}.vcf.gz \
                ~{local_sp_fixed}
            bcftools index --tbi --force ~{prefix}.vcf.gz
            exit 0;
        fi

        ####################################################################
        # 2.1. more prep for fixing undefined VCF INFO/FT/FORMATs
        # get FORMATs in header
        if grep -q -F '##FORMAT=<' header.txt; then
            grep -F '##FORMAT=<' header.txt | awk -F ',' '{print $1}' | sed 's/##FORMAT=<ID=//' | sort > formats.in_header.txt
        else
            touch formats.in_header.txt
        fi
        # get FILTERs in header
        if grep -q -F '##FILTER=<' header.txt; then
            grep -F '##FILTER=<' header.txt | awk -F ',' '{print $1}' | sed 's/##FILTER=<ID=//' | sort > filters.in_header.txt
        else
            touch filters.in_header.txt
        fi
        # get non-flag INFO in header
        if grep -q -F '##INFO=<' header.txt; then
            grep -F '##INFO=<' header.txt | grep -vF 'Type=Flag' | awk -F ',' '{print $1}' | sed 's/##INFO=<ID=//' | sort > non_flag_info.in_header.txt
        else
            touch non_flag_info.in_header.txt
        fi
        # get     flag INFO in header
        if grep -q -F '##INFO=<' header.txt; then
            grep -F '##INFO=<' header.txt | grep  -F 'Type=Flag' | awk -F ',' '{print $1}' | sed 's/##INFO=<ID=//' | sort >     flag_info.in_header.txt
        else
            touch flag_info.in_header.txt
        fi

        # get FORMATs in practice
        awk '{print $9}' body.txt | sort | uniq | sed 's/:/\n/g' | sort | uniq > formats.in_vcf.txt
        # get FILTERs in practice, guard against no 'PASS'
        awk '{print $7}' body.txt | sort | uniq | grep -v "^PASS$" > filters.in_vcf.txt || touch filters.in_vcf.txt

        awk '{print $8}' body.txt | sed 's/;/\n/g' > tmp.info.entries.txt
        if grep -q -F '=' tmp.info.entries.txt; then
            # get non-flag INFO in practicez
            grep -F '=' tmp.info.entries.txt | awk -F '=' '{print $1}' | sort | uniq > non_flag_info.in_vcf.txt
        fi
        if grep -q -vF '=' tmp.info.entries.txt; then
            # get     flag INFO in practice
            awk '{print $8}' body.txt | sed 's/;/\n/g' | grep -vF '=' | sort | uniq > flag_info.in_vcf.txt
        fi
        touch non_flag_info.in_vcf.txt
        touch     flag_info.in_vcf.txt

        echo "I survived. More to go..."

        ####################################################################
        # 2.2. more prep for fixing undefined VCF INFO/FT/FORMATs
        comm -13 formats.in_header.txt formats.in_vcf.txt > missing.formats.txt
        while IFS= read -r line
        do
        echo "##FORMAT=<ID=${line},Number=.,Type=String,Description=\"CALLER DID NOT DEFINE THIS.\">" >> missing.formats.header
        done < missing.formats.txt

        comm -13 filters.in_header.txt filters.in_vcf.txt > missing.filters.txt
        while IFS= read -r line
        do
        echo "##FILTER=<ID=${line},Description=\"CALLER DID NOT DEFINE THIS.\">" >> missing.filters.header
        done < missing.filters.txt

        comm -13 non_flag_info.in_header.txt non_flag_info.in_vcf.txt > missing.non_flag_info.txt
        while IFS= read -r line
        do
        echo "##INFO=<ID=${line},Number=.,Type=String,Description=\"CALLER DID NOT DEFINE THIS.\">" >> missing.non_flag_info.header
        done < missing.non_flag_info.txt

        comm -13 flag_info.in_header.txt flag_info.in_vcf.txt > missing.flag_info.txt
        while IFS= read -r line
        do
        echo "##INFO=<ID=${line},Number=0,Type=Flag,Description=\"CALLER DID NOT DEFINE THIS.\">" >> missing.flag_info.header
        done < missing.flag_info.txt

        ####################################################################
        # 2. actually fix undefined VCF INFO/FT/FORMATs, if necessary
        if  find . -maxdepth 1 -type f -name "missing.*.header" 2>/dev/null | grep -q .; then
            grep "^##" ~{local_sp_fixed} | grep -v "^##[A-Z]" | grep -vF 'contig=' > first_lines.txt
            grep -F "##contig=<ID=" header.txt > contigs.txt
            grep "^#CHROM" ~{local_sp_fixed} > sample.line.txt
            grep "^##" ~{local_sp_fixed} | grep "^##[A-Z]" | sort > existing_definitions.txt
            cat existing_definitions.txt missing.*.header | sort > everything.defined.txt
            cat first_lines.txt contigs.txt everything.defined.txt sample.line.txt > fixed.header.txt
            # print to stdout for checking
            grep -vF "##contig=<ID=" fixed.header.txt

            cat fixed.header.txt body.txt > fixed.vcf
            rm ~{local_sp_fixed}
        else
            mv ~{local_sp_fixed} fixed.vcf
        fi

        ####################################################################
        # 3. fix contigs undefined (in later stages)
        if ~{fix_contigs}; then
            bcftools reheader \
                --fai ~{ref_fasta_fai} \
                -o fixed.and_contigs.vcf \
                fixed.vcf
            mv fixed.and_contigs.vcf fixed.vcf
        fi

        ####################################################################
        # 4. fix occationally unsorted VCF
        bcftools \
            sort \
            --temp-dir tm_sort \
            --output-type z \
            -o ~{prefix}.vcf.gz \
            fixed.vcf
        bcftools index --tbi --force ~{prefix}.vcf.gz
    >>>

    output {
        File sortedVCF = "~{prefix}.vcf.gz"
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             3,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
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
