version 1.0

import "../../structs/Structs.wdl"


task PrepStudyVcfForFlare {

    input {
        File joint_vcf
        File ref_fasta
        File ref_fasta_fai
        String prefix
        Float maf = 0.01
    }

    Int disk_size = 4 * ceil(size(joint_vcf, "GB")) + 50

    command <<<
        set -euxo pipefail

        bcftools view \
            -m2 -M2 -v snps \
            -Ob -o ~{prefix}.snps.pre.bcf \
            ~{joint_vcf}

        bcftools norm \
            -f ~{ref_fasta} \
            -m-any \
            -Ob \
            -o ~{prefix}.snps.bcf \
            ~{prefix}.snps.pre.bcf

        bcftools view \
            -q ~{maf}:minor \
            -Ob -o ~{prefix}.maf.bcf \
            ~{prefix}.snps.bcf

        bcftools index -c -f ~{prefix}.maf.bcf
    >>>

    output {
        File gt_bcf = "~{prefix}.maf.bcf"
        File gt_bcf_csi = "~{prefix}.maf.bcf.csi"
    }

    runtime {
        cpu: 8
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 10
        preemptible: 1
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}


task PrepRefVcfForFlare {

    input {
        File ref_vcf
        String prefix
    }

    Int disk_size = 4 * ceil(size(ref_vcf, "GB")) + 50

    command <<<
        set -euxo pipefail

        bcftools view \
            -m2 -M2 -v snps \
            -Ob -o ~{prefix}.bcf \
            ~{ref_vcf}

        bcftools index -c -f ~{prefix}.bcf
    >>>

    output {
        File ref_bcf = "~{prefix}.bcf"
        File ref_bcf_csi = "~{prefix}.bcf.csi"
    }

    runtime {
        cpu: 4
        memory: "32 GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 10
        preemptible: 1
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}


task IntersectVCFsForFlare {

    input {
        File gt_vcf
        File gt_vcf_index
        File ref_vcf
        File ref_vcf_index
        String prefix
    }

    Float gt_size_gb = size(gt_vcf, "GB")
    Float ref_size_gb = size(ref_vcf, "GB")
    Int disk_size = 4 * ceil(gt_size_gb + ref_size_gb) + 50

    command <<<
        set -euxo pipefail

        mkdir isec
        bcftools isec \
            -p isec \
            -c all \
            -n=2 \
            -w1,2 \
            ~{gt_vcf} \
            ~{ref_vcf}

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' isec/0000.vcf | sort -u > gt.sites
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' isec/0001.vcf | sort -u > ref.sites
        comm -12 gt.sites ref.sites > matched.sites

        n_sites=$(wc -l < matched.sites | tr -d ' ')
        if [ "$n_sites" -eq 0 ]; then
            echo "No overlapping sites with matching REF/ALT between study and reference VCFs" >&2
            echo "Study contigs:" >&2
            bcftools query -f '%CHROM\n' ~{gt_vcf} | sort -u | head -5 >&2
            echo "Reference contigs:" >&2
            bcftools query -f '%CHROM\n' ~{ref_vcf} | sort -u | head -5 >&2
            exit 1
        fi

        bcftools view -T matched.sites -Oz -o ~{prefix}.gt.vcf.gz isec/0000.vcf
        bcftools view -T matched.sites -Oz -o ~{prefix}.ref.vcf.gz isec/0001.vcf
        bcftools index -c -f ~{prefix}.gt.vcf.gz
        bcftools index -c -f ~{prefix}.ref.vcf.gz
    >>>

    output {
        File gt_vcf_out = "~{prefix}.gt.vcf.gz"
        File gt_vcf_csi = "~{prefix}.gt.vcf.gz.csi"
        File ref_vcf_out = "~{prefix}.ref.vcf.gz"
        File ref_vcf_csi = "~{prefix}.ref.vcf.gz.csi"
    }

    runtime {
        cpu: 4
        memory: "32 GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 10
        preemptible: 1
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}


task ThinVCFsForFlare {

    input {
        File gt_vcf
        File ref_vcf
        String prefix
        Int thin_bp = 20000
    }

    Float gt_size_gb = size(gt_vcf, "GB")
    Float ref_size_gb = size(ref_vcf, "GB")
    Int disk_size = 4 * ceil(gt_size_gb + ref_size_gb) + 50

    command <<<
        set -euxo pipefail

        bcftools query -f '%CHROM\t%POS\n' ~{gt_vcf} | \
            awk -v min=~{thin_bp} '
                {
                    if ($1 != c) { c = $1; last = -min }
                    if (last < 0 || $2 - last >= min) { print; last = $2 }
                }
            ' > keep.sites

        n_sites=$(wc -l < keep.sites | tr -d ' ')
        if [ "$n_sites" -eq 0 ]; then
            echo "No sites available for thinning; study VCF has no variant records" >&2
            exit 1
        fi

        bcftools view -T keep.sites -Oz -o ~{prefix}.gt.vcf.gz ~{gt_vcf}
        bcftools view -T keep.sites -Oz -o ~{prefix}.ref.vcf.gz ~{ref_vcf}
        bcftools index -c -f ~{prefix}.gt.vcf.gz
        bcftools index -c -f ~{prefix}.ref.vcf.gz
    >>>

    output {
        File gt_vcf_out = "~{prefix}.gt.vcf.gz"
        File gt_vcf_csi = "~{prefix}.gt.vcf.gz.csi"
        File ref_vcf_out = "~{prefix}.ref.vcf.gz"
        File ref_vcf_csi = "~{prefix}.ref.vcf.gz.csi"
    }

    runtime {
        cpu: 2
        memory: "16 GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 10
        preemptible: 1
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}


task FilterFlareReadySites {

    input {
        File gt_vcf
        File gt_vcf_index
        File ref_vcf
        File ref_vcf_index
        File ref_panel
        String prefix
    }

    Float gt_size_gb = size(gt_vcf, "GB")
    Float ref_size_gb = size(ref_vcf, "GB")
    Int disk_size = 4 * ceil(gt_size_gb + ref_size_gb) + 50

    command <<<
        set -euxo pipefail

        cut -f1 ~{ref_panel} | sort -u > ref.panel.samples
        bcftools query -l ~{gt_vcf} | sort -u > gt.samples
        comm -12 ref.panel.samples gt.samples > overlap.samples

        bcftools view -S ref.panel.samples -Ou ~{ref_vcf} | \
            bcftools query -f '%CHROM\t%POS[\t%GT]\n' | \
            awk '
                BEGIN { OFS = "\t" }
                {
                    bad = 0
                    for (i = 3; i <= NF; i++) {
                        if ($i == "") continue
                        if ($i ~ /\./ || $i ~ /\//) bad = 1
                    }
                    if (!bad) print $1, $2
                }
            ' | sort -k1,1 -k2,2n > ref.complete.sites

        if [ -s overlap.samples ]; then
            bcftools view -S ^overlap.samples -Ou ~{gt_vcf} | \
                bcftools query -f '%CHROM\t%POS[\t%GT]\n' | \
                awk '
                    BEGIN { OFS = "\t" }
                    {
                        bad = 0
                        for (i = 3; i <= NF; i++) {
                            if ($i == "") continue
                            if ($i ~ /\./ || $i ~ /\//) bad = 1
                        }
                        if (!bad) print $1, $2
                    }
                ' | sort -k1,1 -k2,2n > gt.complete.sites
        else
            bcftools view -Ou ~{gt_vcf} | \
                bcftools query -f '%CHROM\t%POS[\t%GT]\n' | \
                awk '
                    BEGIN { OFS = "\t" }
                    {
                        bad = 0
                        for (i = 3; i <= NF; i++) {
                            if ($i == "") continue
                            if ($i ~ /\./ || $i ~ /\//) bad = 1
                        }
                        if (!bad) print $1, $2
                    }
                ' | sort -k1,1 -k2,2n > gt.complete.sites
        fi

        comm -12 ref.complete.sites gt.complete.sites > complete.sites

        n_sites=$(wc -l < complete.sites | tr -d ' ')
        if [ "$n_sites" -eq 0 ]; then
            echo "No sites with complete phased genotypes in both reference panel and study VCF" >&2
            exit 1
        fi

        if [ -s overlap.samples ]; then
            n_study=$(bcftools view -S ^overlap.samples ~{gt_vcf} | bcftools query -l | wc -l | tr -d ' ')
        else
            n_study=$(bcftools query -l ~{gt_vcf} | wc -l | tr -d ' ')
        fi
        if [ "$n_study" -eq 0 ]; then
            echo "No study samples remain after excluding reference panel samples" >&2
            exit 1
        fi

        bcftools view -T complete.sites -Oz -o ~{prefix}.ref.vcf.gz ~{ref_vcf}
        if [ -s overlap.samples ]; then
            bcftools view -S ^overlap.samples -T complete.sites -Oz -o ~{prefix}.gt.vcf.gz ~{gt_vcf}
        else
            bcftools view -T complete.sites -Oz -o ~{prefix}.gt.vcf.gz ~{gt_vcf}
        fi
        bcftools index -c -f ~{prefix}.ref.vcf.gz
        bcftools index -c -f ~{prefix}.gt.vcf.gz
    >>>

    output {
        File gt_vcf_out = "~{prefix}.gt.vcf.gz"
        File gt_vcf_csi = "~{prefix}.gt.vcf.gz.csi"
        File ref_vcf_out = "~{prefix}.ref.vcf.gz"
        File ref_vcf_csi = "~{prefix}.ref.vcf.gz.csi"
    }

    runtime {
        cpu: 4
        memory: "16 GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 10
        preemptible: 1
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}


task MergeGlobalAncestry {

    input {
        Array[File] global_anc_files
        String prefix
    }

    Float input_size_gb = size(global_anc_files, "GB")
    Int disk_size = 2 * ceil(input_size_gb) + 10

    command <<<
        set -euxo pipefail

        zcat ~{sep=' ' global_anc_files} | awk '
            BEGIN { FS = OFS = "\t" }
            FNR == 1 {
                if (nfiles++ == 0) {
                    ncols = NF - 1
                    for (i = 2; i <= NF; i++) {
                        header[i - 1] = $i
                    }
                } else {
                    if (NF - 1 != ncols) {
                        print "Ancestry column count mismatch" > "/dev/stderr"
                        exit 1
                    }
                    for (i = 2; i <= NF; i++) {
                        if ($i != header[i - 1]) {
                            print "Ancestry header mismatch: " $i " vs " header[i - 1] > "/dev/stderr"
                            exit 1
                        }
                    }
                }
                next
            }
            {
                if (NF - 1 != ncols) {
                    print "Malformed row for sample " $1 > "/dev/stderr"
                    exit 1
                }
                sample = $1
                n[sample]++
                for (i = 2; i <= NF; i++) {
                    sum[sample, i - 1] += $i
                }
            }
            END {
                if (ncols == 0) {
                    print "No ancestry columns found in input" > "/dev/stderr"
                    exit 1
                }
                printf "SAMPLE"
                for (j = 1; j <= ncols; j++) {
                    printf "\t%s", header[j]
                }
                print ""
                for (sample in n) {
                    printf "%s", sample
                    for (j = 1; j <= ncols; j++) {
                        printf "\t%.6f", sum[sample, j] / n[sample]
                    }
                    print ""
                }
            }
        ' | sort -k1,1 > ~{prefix}.global_ancestry.tsv
    >>>

    output {
        File global_ancestry_tsv = "~{prefix}.global_ancestry.tsv"
    }

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 10
        preemptible: 1
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}


task Flare {

    meta {
        description: "Local Ancestry Inference"
    }


    input {
        File ref_vcf
        File ref_vcf_index
        File ref_panel
        File test_vcf
        File test_vcf_index
        File plink_map
        String output_prefix
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        Boolean em = true
        File? flare_model

        Int nthreads = 16
        Int mem_gb = 64

        RuntimeAttr? runtime_attr_override
    }

    String model_arg = if defined(flare_model) then "model=~{flare_model}" else ""

    command <<<
        set -euxo pipefail

        java -Xmx~{mem_gb}g -jar /LAI/flare.jar \
            ref=~{ref_vcf} \
            gt=~{test_vcf} \
            map=~{plink_map} \
            ref-panel=~{ref_panel} \
            out=~{output_prefix} \
            nthreads=~{nthreads} \
            em=~{em} \
            ~{model_arg}

    >>>

    output {
        File global_anc = "~{output_prefix}.global.anc.gz"
        File model = "~{output_prefix}.model"
        File anc_vcf = "~{output_prefix}.anc.vcf.gz"
        File log = "~{output_prefix}.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          nthreads,
        mem_gb:             mem_gb,
        disk_gb:            200,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "hangsuunc/flare:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
