version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/TertiaryAnalysis/FunctionalAnnotation.wdl" as FUNK
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ExpandedDrugResistanceMarkerExtraction {

    meta {
        author: "Jonn Smith"
        description: "Extract drug resistance marker information from a single-sample VCF."
    }

    parameter_meta {
        sample_name: "Name of the sample in the given single-sample VCF file."
        vcf: "Input VCF.gz file containing variant calls from which to create an expanded drug resistance marker report."
        snpeff_db: "SNPEff database file to use for functional annotations (if necessary)."

        protein_drug_resistance_list: "List of specific protein changes to include in the expanded drug resistance report."
        gene_drug_resistance_list: "List of genes over which to include any variants in the expanded drug resistance report."

        dir_prefix: "Directory prefix to use for finalized location."
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs."

        do_functional_annotation:    "If true, will perform functional annotation on the input VCF file prior to extracting drug resistance marker information."
    }

    input {
        String sample_name

        File vcf
        File snpeff_db
        File protein_drug_resistance_list
        File gene_drug_resistance_list

        String dir_prefix
        String gcs_out_root_dir

        Boolean do_functional_annotation = true
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ExpandedDrugResistanceMarkerExtraction/~{dir_prefix}"

    if (do_functional_annotation) {
        call FUNK.FunctionallyAnnotateVariants { input: vcf = vcf, snpeff_db = snpeff_db }
    }

    call CallDrugResistanceMutations {
        input:
            vcf = select_first([FunctionallyAnnotateVariants.annotated_vcf, vcf]),
            protein_drug_resistance_list = protein_drug_resistance_list,
            gene_drug_resistance_list = gene_drug_resistance_list,
            prefix = sample_name
    }

    # Finalize data
    String dir = outdir + "/reports"

    call FF.FinalizeToFile as FinalizeDRReportAllMarkers { input: outdir = dir, file = CallDrugResistanceMutations.all_markers }
    call FF.FinalizeToFile as FinalizeDRReportProteinMarkers { input: outdir = dir, file = CallDrugResistanceMutations.protein_coding_markers }

    if (do_functional_annotation) {
        call FF.FinalizeToFile as FinalizeAnnotatedVCF { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.annotated_vcf]) }
        call FF.FinalizeToFile as FinalizeAnnotatedVCFIndex { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.annotated_vcf_index]) }
        call FF.FinalizeToFile as FinalizeSnpEffSummary { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.snpEff_summary]) }
        call FF.FinalizeToFile as FinalizeSnpEffGenes { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.snpEff_genes]) }
    }

    output {
        File drug_res_report_all = FinalizeDRReportAllMarkers.gcs_path
        File drug_res_report_prot_only = FinalizeDRReportProteinMarkers.gcs_path

        File? annotated_vcf = FinalizeAnnotatedVCF.gcs_path
        File? annotated_vcf_index = FinalizeAnnotatedVCFIndex.gcs_path
        File? snpEff_summary = FinalizeSnpEffSummary.gcs_path
        File? snpEff_genes = FinalizeSnpEffGenes.gcs_path
    }
}

task CallDrugResistanceMutations {

    meta {
        author: "Jonn Smith"
        description: "Extract drug resistance marker information from a functionally annotated single-sample VCF."
    }

    parameter_meta {
        vcf: "Input VCF.gz file containing variant calls from which to create an expanded drug resistance marker report."
        protein_drug_resistance_list: "List of specific protein changes to include in the expanded drug resistance report."
        gene_drug_resistance_list: "List of genes over which to include any variants in the expanded drug resistance report."

        prefix: "Prefix to use for output files."

        runtime_attr_override: "Override for default runtime attributes."
    }

    input {
        File vcf
        File protein_drug_resistance_list
        File gene_drug_resistance_list

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 10*ceil(size([vcf, protein_drug_resistance_list, gene_drug_resistance_list], "GB"))

    command <<<
        set -euxo pipefail

        python3 <<CODE

        import gzip

        gene_info = []
        with open("~{gene_drug_resistance_list}", 'r') as f:
            for line in f:
                gene_name, gene_id = line.strip().split()
                gene_info.append((gene_name, gene_id))

        print(f"Gene info:")
        for gene_name, gene_id in gene_info:
            print(f"\t{gene_name}\t{gene_id}")
        print()

        p_change_marker_info = []
        with open("~{protein_drug_resistance_list}", 'r') as f:
            for line in f:
                gene_name, gene_id, prot_change = line.strip().split()
                p_change_marker_info.append((gene_name, gene_id, prot_change))

        print(f"Protein Change Drug Res Info:")
        for gene_name, gene_id, prot_change in p_change_marker_info:
            print(f"\t{gene_name}\t{gene_id}\t{prot_change}")
        print()

        gene_drug_report_all = "~{prefix}.expanded_drug_report.ALL.tsv"
        gene_drug_report_prot = "~{prefix}.expanded_drug_report.PROTEIN_CHANGES_ONLY.tsv"

        def make_ann_dict(ann_info, annotation_fields):
            ann_dicts = []

            vcf_annotations = ann_info.split(",")
            for v in vcf_annotations:
                ann_dict = dict()
                vcf_ann_fields = v.split("|")

                if len(vcf_ann_fields) != len(annotation_fields):
                    print(vcf_ann_fields)
                    raise RuntimeError(f"Error: Field length of annotation is not the length of the number of fields: {len(vcf_ann_fields)} != expected({len(annotation_fields)})")
                for i in range(len(vcf_ann_fields)):
                    ann_dict[annotation_fields[i]] = vcf_ann_fields[i]

                ann_dicts.append(ann_dict)

            return ann_dicts

        def ann_impact_sort_key(x):
            s = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
            return s.index(x[0])

        annotations = []
        with gzip.open("~{vcf}", 'rt') as f:
            for line in f:
                if line.startswith("##INFO=<ID=ANN"):
                    needle = 'Functional annotations: '
                    i1 = line.find(needle)
                    i2 = line.find(' ">')
                    annotation_fields = line[i1+len(needle):i2].replace("'", "").split(" | ")
                    continue
                elif line.startswith("#"):
                    continue
                chrom, pos, idd, ref, alt, qual, flt, info, gt_f, gt_d = line.strip().split("\t")
                infos = {}
                for i in info.split(";"):
                    if "=" in i:
                        k, v = i.split("=")
                        infos[k] = v
                ann_dicts = make_ann_dict(infos["ANN"], annotation_fields)

                # Get genotype info:
                gt_f = gt_f.split(":")
                gt_d = gt_d.split(":")

                i = gt_f.index("GT")
                gt = gt_d[i]

                base_out_data = [chrom, pos, ref, alt, gt]

                # Now check if our Gene level drug markers:
                gene_annotations = []
                for gene_name, gene_id in gene_info:
                    for ann_dict in ann_dicts:
                        if gene_id in ann_dict["Gene_Name"] or gene_id in ann_dict["Gene_ID"]:
                            gene_annotations.append(tuple(base_out_data + [ann_dict[a] for a in annotation_fields]))

                for g in gene_annotations:
                    annotations.append(g)

                # Now check for our protein change string drug markers:
                for gene_name, gene_id, prot_change in p_change_marker_info:
                    for ann_dict in ann_dicts:
                        if gene_id in ann_dict["Gene_Name"] or gene_id in ann_dict["Gene_ID"]:
                            if len(ann_dict["HGVS.p"]) > 0 and ann_dict["HGVS.p"] == prot_change:
                                annotations.append(tuple(base_out_data + [ann_dict[a] for a in annotation_fields]))

        header = "Chrom\tPos\tRef\tAlt\tGT\t" + "\t".join(annotation_fields)
        with open(gene_drug_report_all, 'w') as f:
            f.write(f"{header}\n")
            for a in annotations:
                f.write("\t".join(a))
                f.write("\n")

        with open(gene_drug_report_prot, 'w') as f:
            f.write(f"{header}\n")
            for a in annotations:
                f.write("\t".join(a))
                f.write("\n")

        print('Done')
        CODE

    >>>

    output {
        File all_markers = "~{prefix}.expanded_drug_report.ALL.tsv"
        File protein_coding_markers = "~{prefix}.expanded_drug_report.PROTEIN_CHANGES_ONLY.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/snpeff:5.1d--hdfd78af_0"
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
