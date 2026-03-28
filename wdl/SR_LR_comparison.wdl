version 1.0

struct RuntimeAttr {
    Int? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

workflow CompareLRSR {
    input {
        File sr_sv_vcf
        File sr_sv_vcf_tbi
        File lr_sv_vcf
        File lr_sv_vcf_tbi
        File reference_fa
        String prefix
        Boolean compare_samples = false
        String truvari_params_ins = "--pctseq 0 --pctsize 0.1 -r 100 --pctovl 0 --pick multi --dup-to-ins"
        String truvari_params_delinv_small = "--pctseq 0.5 --pctsize 0.5 -r 10000 --chunksize 10000 --pctovl 0.1 --pick multi --dup-to-ins"
        String truvari_params_delinv_large = "--pctseq 0.5 --pctsize 0.5 -r 20000 --chunksize 20000 --pctovl 0.1 --pick multi --dup-to-ins"
        String truvari_params_dup_small = "--pctseq 0.5 --pctsize 0.5 -r 10000 --chunksize 10000 --pctovl 0 --pick multi --dup-to-ins"
        String truvari_params_dup_large = "--pctseq 0.5 --pctsize 0.5 -r 20000 --chunksize 20000 --pctovl 0 --pick multi --dup-to-ins"
    }

    call FilterSRVCF {
        input:
            vcf = sr_sv_vcf,
            vcf_tbi = sr_sv_vcf_tbi,
            prefix = prefix,
            compare_samples = compare_samples,
            lr_vcf = lr_sv_vcf
    }

    call TruvariBench as TruvariBench_INS {
        input:
            truth_vcf = lr_sv_vcf,
            truth_vcf_index = lr_sv_vcf_tbi,
            comp_vcf = FilterSRVCF.ins_vcf,
            comp_vcf_index = FilterSRVCF.ins_vcf_tbi,
            ref_fa = reference_fa,
            params = truvari_params_ins
    }

    call TruvariBench as TruvariBench_largeDELINV {
        input:
            truth_vcf = lr_sv_vcf,
            truth_vcf_index = lr_sv_vcf_tbi,
            comp_vcf = FilterSRVCF.delinv_large_vcf,
            comp_vcf_index = FilterSRVCF.delinv_large_vcf_tbi,
            ref_fa = reference_fa,
            params = truvari_params_delinv_large
    }

    call TruvariBench as TruvariBench_smallDELINV {
        input:
            truth_vcf = lr_sv_vcf,
            truth_vcf_index = lr_sv_vcf_tbi,
            comp_vcf = FilterSRVCF.delinv_small_vcf,
            comp_vcf_index = FilterSRVCF.delinv_small_vcf_tbi,
            ref_fa = reference_fa,
            params = truvari_params_delinv_small
    }

    call TruvariBench as TruvariBench_largeDUP {
        input:
            truth_vcf = lr_sv_vcf,
            truth_vcf_index = lr_sv_vcf_tbi,
            comp_vcf = FilterSRVCF.dup_large_vcf,
            comp_vcf_index = FilterSRVCF.dup_large_vcf_tbi,
            ref_fa = reference_fa,
            params = truvari_params_dup_large
    }

    call TruvariBench as TruvariBench_smallDUP {
        input:
            truth_vcf = lr_sv_vcf,
            truth_vcf_index = lr_sv_vcf_tbi,
            comp_vcf = FilterSRVCF.dup_small_vcf,
            comp_vcf_index = FilterSRVCF.dup_small_vcf_tbi,
            ref_fa = reference_fa,
            params = truvari_params_dup_small
    }

    call FilterDups as FilterDupsLarge {
        input:
            vcf = TruvariBench_largeDUP.tp_base_vcf,
            vcf_index = TruvariBench_largeDUP.tp_base_vcf_index,
            vcf_type = "base",
            prefix = prefix + "_LARGEDUP_FILTERED"
    }

    call FilterDups as FilterDupsSmall {
        input:
            vcf = TruvariBench_smallDUP.tp_base_vcf,
            vcf_index = TruvariBench_smallDUP.tp_base_vcf_index,
            vcf_type = "base",
            prefix = prefix + "_SMALLDUP_FILTERED"
    }

    Array[String] annotation_labels = ["LARGEDUP", "SMALLDUP", "INS", "LARGEDELINV", "SMALLDELINV"]
    Array[File] base_annotation_vcfs = [
        FilterDupsLarge.filtered_vcf,
        FilterDupsSmall.filtered_vcf,
        TruvariBench_INS.tp_base_vcf,
        TruvariBench_largeDELINV.tp_base_vcf,
        TruvariBench_smallDELINV.tp_base_vcf
    ]
    Array[File] base_annotation_vcf_indices = [
        FilterDupsLarge.filtered_vcf_index,
        FilterDupsSmall.filtered_vcf_index,
        TruvariBench_INS.tp_base_vcf_index,
        TruvariBench_largeDELINV.tp_base_vcf_index,
        TruvariBench_smallDELINV.tp_base_vcf_index
    ]
    Array[File] comp_annotation_vcfs = [
        TruvariBench_largeDUP.tp_comp_vcf,
        TruvariBench_smallDUP.tp_comp_vcf,
        TruvariBench_INS.tp_comp_vcf,
        TruvariBench_largeDELINV.tp_comp_vcf,
        TruvariBench_smallDELINV.tp_comp_vcf
    ]
    Array[File] comp_annotation_vcf_indices = [
        TruvariBench_largeDUP.tp_comp_vcf_index,
        TruvariBench_smallDUP.tp_comp_vcf_index,
        TruvariBench_INS.tp_comp_vcf_index,
        TruvariBench_largeDELINV.tp_comp_vcf_index,
        TruvariBench_smallDELINV.tp_comp_vcf_index
    ]

    scatter (i in range(length(annotation_labels))) {
        call MakeAnnotations as BuildAnnotationTable {
            input:
                base_vcf = base_annotation_vcfs[i],
                base_vcf_index = base_annotation_vcf_indices[i],
                comp_vcf = comp_annotation_vcfs[i],
                comp_vcf_index = comp_annotation_vcf_indices[i],
                output_prefix = prefix + "_" + annotation_labels[i],
                require_sample_overlap = compare_samples
        }
    }

    call MergeAnnotations {
        input:
            annotation_tsvs = BuildAnnotationTable.annotation_tsv,
            prefix = prefix
    }

    call ApplyAnnotations {
        input:
            orig_vcf = lr_sv_vcf,
            orig_vcf_index = lr_sv_vcf_tbi,
            anno_tsv = MergeAnnotations.merged_annotation_tsv,
            prefix = prefix
    }

    output {
        File annotated_VCF = ApplyAnnotations.annotated_vcf
        File annotated_VCF_index = ApplyAnnotations.annotated_vcf_index
        File merged_annotation_tsv = MergeAnnotations.merged_annotation_tsv
        Array[File] per_bucket_annotation_tsvs = BuildAnnotationTable.annotation_tsv
    }
}

task FilterSRVCF {
    input {
        File vcf
        File vcf_tbi
        String prefix
        Boolean compare_samples
        File lr_vcf
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2 * ceil(size(vcf, "GB") + size(lr_vcf, "GB")) + 10

    command <<<
        set -euxo pipefail

        if [ ~{compare_samples} = true ]; then
            bcftools query -l ~{lr_vcf} | sort -u > LR_samples.txt
            bcftools query -l ~{vcf} | sort -u > SR_samples.txt
            comm -12 LR_samples.txt SR_samples.txt > shared_samples.txt

            if [ ! -s shared_samples.txt ]; then
                echo "ERROR: compare_samples=true but no overlapping sample IDs were found between LR and SR VCFs" >&2
                echo "LR samples:" >&2
                cat LR_samples.txt >&2
                echo "SR samples:" >&2
                cat SR_samples.txt >&2
                exit 1
            fi

            bcftools view -S shared_samples.txt ~{vcf} -Oz -o filtered.vcf.gz
            tabix -p vcf filtered.vcf.gz
        else
            ln -s ~{vcf} filtered.vcf.gz
            ln -s ~{vcf_tbi} filtered.vcf.gz.tbi
        fi

        bcftools view -Oz -f PASS -i 'SVTYPE=="DEL" | SVTYPE=="INV"' filtered.vcf.gz -o ~{prefix}.DELINV.vcf.gz
        bcftools view -Oz -f PASS -i 'SVTYPE=="DUP"' filtered.vcf.gz -o ~{prefix}.DUP.vcf.gz
        bcftools view -Oz -f PASS -i 'SVTYPE=="INS"' filtered.vcf.gz -o ~{prefix}.INS.vcf.gz

        bcftools view -Oz -i 'SVLEN>=5000 | SVLEN<=-5000' ~{prefix}.DELINV.vcf.gz -o ~{prefix}.DELINV_over5kb.vcf.gz
        bcftools view -Oz -i 'SVLEN>=5000 | SVLEN<=-5000' ~{prefix}.DUP.vcf.gz -o ~{prefix}.DUP_over5kb.vcf.gz
        bcftools view -Oz -i 'SVLEN<5000 & SVLEN>-5000' ~{prefix}.DELINV.vcf.gz -o ~{prefix}.DELINV_under5kb.vcf.gz
        bcftools view -Oz -i 'SVLEN<5000 & SVLEN>-5000' ~{prefix}.DUP.vcf.gz -o ~{prefix}.DUP_under5kb.vcf.gz

        tabix -p vcf ~{prefix}.DELINV_under5kb.vcf.gz
        tabix -p vcf ~{prefix}.DELINV_over5kb.vcf.gz
        tabix -p vcf ~{prefix}.DUP_under5kb.vcf.gz
        tabix -p vcf ~{prefix}.DUP_over5kb.vcf.gz
        tabix -p vcf ~{prefix}.INS.vcf.gz
    >>>

    output {
        File delinv_small_vcf = "~{prefix}.DELINV_under5kb.vcf.gz"
        File delinv_small_vcf_tbi = "~{prefix}.DELINV_under5kb.vcf.gz.tbi"
        File delinv_large_vcf = "~{prefix}.DELINV_over5kb.vcf.gz"
        File delinv_large_vcf_tbi = "~{prefix}.DELINV_over5kb.vcf.gz.tbi"
        File dup_small_vcf = "~{prefix}.DUP_under5kb.vcf.gz"
        File dup_small_vcf_tbi = "~{prefix}.DUP_under5kb.vcf.gz.tbi"
        File dup_large_vcf = "~{prefix}.DUP_over5kb.vcf.gz"
        File dup_large_vcf_tbi = "~{prefix}.DUP_over5kb.vcf.gz.tbi"
        File ins_vcf = "~{prefix}.INS.vcf.gz"
        File ins_vcf_tbi = "~{prefix}.INS.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 1,
        disk_gb: disk_size,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0,
        docker: "quay.io/ymostovoy/lr-utils-basic:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task TruvariBench {
    input {
        File truth_vcf
        File truth_vcf_index
        File comp_vcf
        File comp_vcf_index
        File ref_fa
        String params = ""
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(truth_vcf, "GB") + size(comp_vcf, "GB")) * 3 + 50

    command <<<
        set -euxo pipefail

        truvari bench -b ~{truth_vcf} -c ~{comp_vcf} -f ~{ref_fa} -o outdir ~{params}
        tabix -f -p vcf outdir/tp-base.vcf.gz
        tabix -f -p vcf outdir/tp-comp.vcf.gz
        tabix -f -p vcf outdir/fp.vcf.gz
        tabix -f -p vcf outdir/fn.vcf.gz
    >>>

    output {
        File tp_base_vcf = "outdir/tp-base.vcf.gz"
        File tp_base_vcf_index = "outdir/tp-base.vcf.gz.tbi"
        File tp_comp_vcf = "outdir/tp-comp.vcf.gz"
        File tp_comp_vcf_index = "outdir/tp-comp.vcf.gz.tbi"
        File fp_vcf = "outdir/fp.vcf.gz"
        File fp_vcf_index = "outdir/fp.vcf.gz.tbi"
        File fn_vcf = "outdir/fn.vcf.gz"
        File fn_vcf_index = "outdir/fn.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: disk_size,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0,
        docker: "quay.io/ymostovoy/lr-truvari:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task FilterDups {
    input {
        File vcf
        File vcf_index
        String vcf_type
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2 * ceil(size(vcf, "GB")) + 10

    command <<<
        set -euo pipefail

        python3 - ~{vcf} ~{prefix}.vcf ~{vcf_type} <<'EOF'
import sys
from pysam import VariantFile

vcf_in = VariantFile(sys.argv[1], 'r')
vcf_out = VariantFile(sys.argv[2], 'w', header=vcf_in.header)
truvari_type = sys.argv[3]

if truvari_type not in ['base', 'comp']:
    sys.stderr.write("Please enter either base or comp\n")
    sys.exit(1)

for rec in vcf_in.fetch():
    svtype = rec.info.get('SVTYPE')
    if isinstance(svtype, tuple):
        svtype = svtype[0]

    if svtype not in ['DUP', 'INS']:
        vcf_out.write(rec)
        continue

    svlen = rec.info.get('SVLEN')
    if isinstance(svlen, tuple):
        svlen = svlen[0]
    svlen = int(svlen)

    start_distance = rec.info.get('StartDistance', 0)
    size_diff = rec.info.get('SizeDiff', 0)
    if isinstance(start_distance, tuple):
        start_distance = start_distance[0]
    if isinstance(size_diff, tuple):
        size_diff = size_diff[0]

    if truvari_type == 'base':
        matched_start = rec.start - int(start_distance)
        matched_size = svlen - int(size_diff)
    else:
        matched_start = rec.start + int(start_distance)
        matched_size = svlen + int(size_diff)

    matched_end = matched_start + matched_size
    cur_end = rec.start + svlen

    if min(cur_end, matched_end) > max(rec.start, matched_start):
        vcf_out.write(rec)

vcf_in.close()
vcf_out.close()
EOF

        bcftools view -Oz -o ~{prefix}.filtered.vcf.gz ~{prefix}.vcf
        bcftools index -t ~{prefix}.filtered.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.filtered.vcf.gz"
        File filtered_vcf_index = "~{prefix}.filtered.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: disk_size,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0,
        docker: "quay.io/ymostovoy/lr-utils-basic:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task MakeAnnotations {
    input {
        File base_vcf
        File base_vcf_index
        File comp_vcf
        File comp_vcf_index
        String output_prefix
        Boolean require_sample_overlap = false
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(base_vcf, "GB") + size(comp_vcf, "GB")) * 2 + 10

    command <<<
        set -euo pipefail

        python3 - ~{base_vcf} ~{comp_vcf} ~{output_prefix}.annotations.tsv ~{require_sample_overlap} <<'EOF'
import sys
from collections import defaultdict
from pysam import VariantFile


def parse_match_id(value):
    if value is None:
        return None, None
    if isinstance(value, (tuple, list)):
        if len(value) >= 2:
            return str(value[0]), str(value[1])
        if len(value) == 1 and ',' in str(value[0]):
            parts = str(value[0]).split(',', 1)
            return parts[0], parts[1]
        return None, None
    text = str(value)
    if ',' not in text:
        return None, None
    parts = text.split(',', 1)
    return parts[0], parts[1]


def alt_samples(record):
    carriers = set()
    for sample, sample_data in record.samples.items():
        gt = sample_data.get('GT')
        if gt is None:
            continue
        for allele in gt:
            if allele not in (None, 0):
                carriers.add(sample)
                break
    return carriers


def format_af(record):
    if 'AF' in record.info:
        value = record.info['AF']
        if isinstance(value, (tuple, list)):
            return ','.join(str(v) for v in value)
        return str(value)
    if 'AC' in record.info and 'AN' in record.info:
        ac = record.info['AC']
        an = record.info['AN']
        if isinstance(ac, (tuple, list)):
            ac = ac[0]
        if isinstance(an, (tuple, list)):
            an = an[0]
        try:
            if float(an) != 0:
                return str(float(ac) / float(an))
        except Exception:
            pass
    return '.'


base_vcf_path = sys.argv[1]
comp_vcf_path = sys.argv[2]
out_tsv = sys.argv[3]
require_sample_overlap = sys.argv[4].lower() == 'true'

base_to_comp_ids = defaultdict(set)
comp_to_base_ids = defaultdict(set)
comp_records_by_base = defaultdict(list)
comp_records_by_comp = defaultdict(list)

comp_vcf = VariantFile(comp_vcf_path)
for record in comp_vcf.fetch():
    base_id, comp_id = parse_match_id(record.info.get('MatchId'))
    if base_id is None or comp_id is None:
        continue

    record_info = {
        'comp_id': comp_id,
        'sr_id': record.id if record.id is not None else '.',
        'sr_af': format_af(record),
        'samples': alt_samples(record)
    }
    base_key = (record.chrom, base_id)
    comp_key = (record.chrom, comp_id)

    base_to_comp_ids[base_key].add(comp_id)
    comp_to_base_ids[comp_key].add(base_id)
    comp_records_by_base[base_key].append(record_info)
    comp_records_by_comp[comp_key].append(record_info)
comp_vcf.close()

base_vcf = VariantFile(base_vcf_path)
base_records = []
for record in base_vcf.fetch():
    base_id, comp_id = parse_match_id(record.info.get('MatchId'))
    if base_id is None:
        continue

    base_key = (record.chrom, base_id)
    if comp_id is not None:
        comp_key = (record.chrom, comp_id)
        base_to_comp_ids[base_key].add(comp_id)
        comp_to_base_ids[comp_key].add(base_id)

    base_records.append(record)

with open(out_tsv, 'w') as out:
    for record in base_records:
        base_id, _ = parse_match_id(record.info.get('MatchId'))
        if base_id is None:
            continue

        base_key = (record.chrom, base_id)
        candidate_matches = []

        # Direct tp-comp records that explicitly point to this base MatchId.
        candidate_matches.extend(comp_records_by_base.get(base_key, []))

        # Additional tp-comp records recovered by walking through comp-side MatchIds
        # observed from either tp-base or tp-comp, which makes the logic robust to
        # unbalanced one-to-many reporting.
        for comp_id in base_to_comp_ids.get(base_key, set()):
            comp_key = (record.chrom, comp_id)
            candidate_matches.extend(comp_records_by_comp.get(comp_key, []))
            for linked_base_id in comp_to_base_ids.get(comp_key, set()):
                linked_base_key = (record.chrom, linked_base_id)
                candidate_matches.extend(comp_records_by_base.get(linked_base_key, []))

        if not candidate_matches:
            continue

        base_carriers = alt_samples(record)
        kept = []
        for match in candidate_matches:
            if require_sample_overlap and len(base_carriers.intersection(match['samples'])) == 0:
                continue
            kept.append((match['sr_id'], match['sr_af']))

        if not kept:
            continue

        kept = sorted(set(kept), key=lambda x: (x[0], x[1]))
        alt_value = ','.join(record.alts) if record.alts is not None else '.'
        sr_ids = ','.join(item[0] for item in kept)
        sr_afs = ','.join(item[1] for item in kept)
        lr_id = record.id if record.id is not None else '.'
        out.write(f"{record.chrom}\t{record.pos}\t{record.ref}\t{alt_value}\t{lr_id}\t{sr_ids}\t{sr_afs}\n")
base_vcf.close()
EOF
    >>>

    output {
        File annotation_tsv = "~{output_prefix}.annotations.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: disk_size,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0,
        docker: "quay.io/ymostovoy/lr-utils-basic:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task MergeAnnotations {
    input {
        Array[File] annotation_tsvs
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{sep=' ' annotation_tsvs} | awk 'NF > 0' | sort -u | sort -k1,1 -k2,2n > ~{prefix}_merged_annotation.tsv
    >>>

    output {
        File merged_annotation_tsv = "~{prefix}_merged_annotation.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0,
        docker: "quay.io/ymostovoy/lr-utils-basic:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task ApplyAnnotations {
    input {
        File orig_vcf
        File orig_vcf_index
        File anno_tsv
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2 * ceil(size(orig_vcf, "GB") + size(anno_tsv, "GB")) + 10

    command <<<
        set -euo pipefail

        python3 - ~{orig_vcf} ~{anno_tsv} ~{prefix}.annotated.vcf <<'EOF'
import sys
from pysam import VariantFile

vcf_in = VariantFile(sys.argv[1])
anno_tsv = sys.argv[2]
out_vcf = sys.argv[3]

annotations = {}
with open(anno_tsv, 'r') as handle:
    for line in handle:
        chrom, pos, ref, alt, lr_id, sr_ids, sr_afs = line.rstrip('\n').split('\t')
        annotations[(chrom, pos, ref, alt, lr_id)] = (tuple(sr_ids.split(',')), tuple(sr_afs.split(',')))

header = vcf_in.header.copy()
header.add_line('##INFO=<ID=SR_MATCH_IDS,Number=.,Type=String,Description="Matched short-read structural variant IDs from Truvari benchmarking">')
header.add_line('##INFO=<ID=SR_MATCH_AFS,Number=.,Type=String,Description="Allele frequencies of matched short-read structural variants, parallel to SR_MATCH_IDS">')
header.add_line('##INFO=<ID=SR_MATCH_COUNT,Number=1,Type=Integer,Description="Number of matched short-read structural variants">')

vcf_out = VariantFile(out_vcf, 'w', header=header)
for record in vcf_in:
    alt = ','.join(record.alts) if record.alts is not None else '.'
    lr_id = record.id if record.id is not None else '.'
    key = (record.chrom, str(record.pos), record.ref, alt, lr_id)
    if key in annotations:
        sr_ids, sr_afs = annotations[key]
        record.info['SR_MATCH_IDS'] = sr_ids
        record.info['SR_MATCH_AFS'] = sr_afs
        record.info['SR_MATCH_COUNT'] = len(sr_ids)
    vcf_out.write(record)

vcf_out.close()
vcf_in.close()
EOF

        bgzip -c ~{prefix}.annotated.vcf > ~{prefix}.annotated.vcf.gz
        tabix -f -p vcf ~{prefix}.annotated.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.annotated.vcf.gz"
        File annotated_vcf_index = "~{prefix}.annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: disk_size,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0,
        docker: "quay.io/ymostovoy/lr-utils-basic:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}
