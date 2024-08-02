version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/TertiaryAnalysis/FunctionalAnnotation.wdl" as FUNK
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PfalciparumTypeDrugResistanceMarkers {

    meta {
        description: "Workflow to generate a report of drug resistance markers"
    }
    parameter_meta {
        vcf: "VCF file to process"
        snpeff_db: "SnpEff database for functional annotation"
        drug_resistance_list: "List of drug resistance markers for which to search"

        dir_prefix: "Prefix for output directory"
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs.  If no bucket is given, outputs will not be finalized and instead will remain in their native execution location."

        do_functional_annotation: "Whether to perform functional annotation"
    }

    input {
        File vcf
        File snpeff_db
        File drug_resistance_list

        String dir_prefix
        String? gcs_out_root_dir

        Boolean do_functional_annotation = true
    }

    if (do_functional_annotation) {
        call FUNK.FunctionallyAnnotateVariants { input: vcf = vcf, snpeff_db = snpeff_db }
    }

    call CallDrugResistanceMutations {
        input:
            vcf = select_first([FunctionallyAnnotateVariants.annotated_vcf, vcf]),
            drug_resistance_list = drug_resistance_list
    }

    call CreateDrugResistanceSummary {
        input:
            raw_drug_resistance_report = CallDrugResistanceMutations.report,
            prefix = basename(basename(vcf, ".gz"), ".vcf")
    }

    # Finalize data
    if (defined(gcs_out_root_dir)) {

        String concrete_gcs_out_root_dir = select_first([gcs_out_root_dir])

        String outdir = sub(concrete_gcs_out_root_dir, "/$", "") + "/PfalciparumTypeDrugResistanceMarkers/~{dir_prefix}"
        String dir = outdir + "/reports"

        call FF.FinalizeToFile as FinalizeDRReport { input: outdir = dir, file = CallDrugResistanceMutations.report }
        call FF.FinalizeToFile as FinalizeDRSummary { input: outdir = dir, file = CreateDrugResistanceSummary.resistance_summary }

        if (do_functional_annotation) {
            call FF.FinalizeToFile as FinalizeAnnotatedVCF { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.annotated_vcf]) }
            call FF.FinalizeToFile as FinalizeAnnotatedVCFIndex { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.annotated_vcf_index]) }
            call FF.FinalizeToFile as FinalizeSnpEffSummary { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.snpEff_summary]) }
            call FF.FinalizeToFile as FinalizeSnpEffGenes { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.snpEff_genes]) }
        }

        File final_annotated_vcf = select_first([FinalizeAnnotatedVCF.gcs_path, FunctionallyAnnotateVariants.annotated_vcf])
        File final_annotated_vcf_index = select_first([FinalizeAnnotatedVCFIndex.gcs_path, FunctionallyAnnotateVariants.annotated_vcf_index])
        File final_snpeff_summary = select_first([FinalizeSnpEffSummary.gcs_path, FunctionallyAnnotateVariants.snpEff_summary])
        File final_snpeff_genes = select_first([FinalizeSnpEffGenes.gcs_path, FunctionallyAnnotateVariants.snpEff_genes])
    }

    output {
        File drug_resistance_summary = select_first([FinalizeDRSummary.gcs_path, CreateDrugResistanceSummary.resistance_summary])
        File raw_drug_res_report = select_first([FinalizeDRReport.gcs_path, CallDrugResistanceMutations.report])

        String chloroquine_status   = CreateDrugResistanceSummary.chloroquine_status
        String pyrimethamine_status = CreateDrugResistanceSummary.pyrimethamine_status
        String sulfadoxine_status   = CreateDrugResistanceSummary.sulfadoxine_status
        String mefloquine_status    = CreateDrugResistanceSummary.mefloquine_status
        String artemisinin_status   = CreateDrugResistanceSummary.artemisinin_status
        String piperaquine_status   = CreateDrugResistanceSummary.piperaquine_status

        File? annotated_vcf = final_annotated_vcf
        File? annotated_vcf_index = final_annotated_vcf_index
        File? snpEff_summary = final_snpeff_summary
        File? snpEff_genes = final_snpeff_genes
    }
}

task CallDrugResistanceMutations {
    input {
        File vcf
        File drug_resistance_list

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size([vcf, drug_resistance_list], "GB"))
    String prefix = basename(basename(vcf, ".gz"), ".vcf")

    command <<<
        set -x

        while read LINE; do
            GENE_NAME=$(echo $LINE | awk '{print $1}')
            GENE_ID=$(echo $LINE | awk '{print $2}')
            MUTATION=$(echo $LINE | awk '{print $3}')

            zcat ~{vcf} | grep $GENE_ID | grep $MUTATION | wc -l | \
                awk -v gene_name=$GENE_NAME -v gene_id=$GENE_ID -v mutation=$MUTATION \
                    '{ print gene_name, gene_id, mutation, ($1 > 0) ? "present" : "absent" }' | \
                tee -a ~{prefix}.drug_resistance_report.txt
        done < ~{drug_resistance_list}
    >>>

    output {
        File report = "~{prefix}.raw_drug_resistance_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
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

task CreateDrugResistanceSummary {
    meta {
        desciption: "Create a drug resistance report based on the given raw drug resistance loci report."
    }

    input {
        File raw_drug_resistance_report
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    String outfile_name = "~{prefix}.drug_resistance_summary.txt"

    Int disk_size = 1 + 4*ceil(size(raw_drug_resistance_report, "GB"))

    command <<<
        python3 <<CODE

        import re
        from enum import Enum

        # Set up some functions to determine drug sensitivity:
        pchange_regex = re.compile(r"""p\.([A-z]+)([0-9]+)([A-z]+)""")
        AA_3_2 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

        DrugSensitivity = Enum("DrugSensitivity", ["UNDETERMINED", "SENSITIVE", "RESISTANT"])

        def parse_pchange(pchange):
            old, pos, new = pchange_regex.match(pchange).groups()
            return AA_3_2[old.upper()], int(pos), AA_3_2[new.upper()]

        def get_chloroquine_sensitivity(dr_report):
        # Locus utilized: PF3D7_0709000 (crt)
        # Codon: 76
        # Workflow:
        # Step Genetic change Interpretation Classification
        # 1 76 K/T heterozygote Heterozygous mutant Undetermined
        # 2 76 missing Missing Undetermined
        # 3 K76 Wild type Sensitive
        # 4 76T Mutant Resistant
        # 5 otherwise Unknown mutant Undetermined
            with open(dr_report, 'r') as f:
                for line in f:
                    # We only care about this locus for this drug:
                    if line.startswith("pfcrt PF3D7_0709000"):
                        gene, locus, pchange, marker = line.strip().split(" ")
                        old, pos, new = parse_pchange(pchange)
                        if pos == 76:
                            if (old == "L") and (new == "T") and (marker == "absent"):
                                return DrugSensitivity.SENSITIVE
                            elif (new == "T") and (marker == "present"):
                                return DrugSensitivity.RESISTANT

            return DrugSensitivity.UNDETERMINED

        def get_pyrimethamine_sensitivity(dr_report):
        # Locus utilized: PF3D7_0417200 (dhfr)
        # Codon: 108
        # Workflow:
        # Step Genetic change Interpretation Classification
        # 1 108 S/N heterozygote Heterozygous mutant Undetermined
        # 2 108 missing Missing Undetermined
        # 3 S108 Wild type Sensitive
        # 4 108N Mutant Resistant
        # 5 otherwise Unknown mutant Undetermined

            with open(dr_report, 'r') as f:
                for line in f:
                    # We only care about this locus for this drug:
                    if line.startswith("pfdhfr PF3D7_0417200"):
                        gene, locus, pchange, marker = line.strip().split(" ")
                        old, pos, new = parse_pchange(pchange)
                        if pos == 108:
                            if (old == "S") and (new == "N") and (marker == "absent"):
                                return DrugSensitivity.SENSITIVE
                            elif (old == "S") and (new == "N") and (marker == "present"):
                                return DrugSensitivity.RESISTANT
                            elif (new == "N") and (marker == "present"):
                                return DrugSensitivity.RESISTANT
                            elif marker == "absent":
                                return DrugSensitivity.SENSITIVE

            return DrugSensitivity.UNDETERMINED


        def get_sulfadoxine_sensitivity(dr_report):
        # Locus utilized: PF3D7_0810800 (dhps)
        # Codon: 437
        # Workflow:
        # Step Genetic change Interpretation Classification
        # 1 437 A/G heterozygote Heterozygous mutant Undetermined
        # 2 437 missing Missing Undetermined
        # 3 A437 Wild type Sensitive
        # 4 437G Mutant Resistant
        # 5 otherwise Unknown mutant Undetermined

            with open(dr_report, 'r') as f:
                for line in f:
                    # We only care about this locus for this drug:
                    if line.startswith("pfdhps PF3D7_0810800"):
                        gene, locus, pchange, marker = line.strip().split(" ")
                        old, pos, new = parse_pchange(pchange)
                        if pos == 437:
                            if (old == "A") and (new == "G") and (marker == "absent"):
                                return DrugSensitivity.SENSITIVE
                            elif (old == "A") and (new == "G") and (marker == "present"):
                                return DrugSensitivity.RESISTANT
                            elif (new == "G") and (marker == "present"):
                                return DrugSensitivity.RESISTANT
                            elif marker == "absent":
                                return DrugSensitivity.SENSITIVE

            return DrugSensitivity.UNDETERMINED

        def get_mefloquine_sensitivity(dr_report):
        # Locus utilized: PF3D7_0523000 (mdr1)
        # Codons: Amplification status of whole gene
        # Workflow:
        # Step Genetic change Interpretation Classification
        # 1 Missing Missing Undetermined
        # 2 Heterozygous duplication Heterozygous mutant Undetermined
        # 3 Single copy Wild type Sensitive
        # 4 Multiple copies Mutant Resistant

            # Currently we can't determine this.
            # We need to get CNV calling working first.

            return DrugSensitivity.UNDETERMINED

        def get_artemisinin_sensitivity(dr_report):
        # Locus utilized: PF3D7_1343700 (kelch13)
        # Codons: 349-726 (BTB/POZ and propeller domains)
        # Workflow:
        # Step Genetic change Interpretation Classification
        # 1 Homozygous non-synonymous mutations in the kelch13 BTB/POZ and propeller
        # domain classified by the World Health Organisation as associated with delayed
        # parasite clearance
        # Mutant – associated with delayed clearance Resistant
        # 2 Heterozygous non-synonymous mutations in the kelch13 BTB/POZ and
        # propeller domain classified by the World Health Organisation as associated
        # with delayed parasite clearance
        # Mutant - heterozygous Undetermined
        # 3 578S as homozygous Mutant - not associated Sensitive
        # 4 Any missing call in amino acids 349-726 Missing Undetermined
        # 5 No non-reference calls in amino acids 349-726 Wild-type Sensitive
        # 6 otherwise Mutant - not in WHO list Undetermined

            with open(dr_report, 'r') as f:
                for line in f:
                    # We only care about this locus for this drug:
                    if line.startswith("pfkelch13 PF3D7_1343700"):
                        gene, locus, pchange, marker = line.strip().split(" ")
                        old, pos, new = parse_pchange(pchange)

                        has_non_ref = False
                        has_variants = False
                        if 349 <= pos <= 726:
                            if (old != new) and (marker == "present"):
                                return DrugSensitivity.RESISTANT
                            elif (new == "S") and (marker == "present"):
                                return DrugSensitivity.SENSITIVE
                            elif marker == "present":
                                has_non_ref = True
                            has_variants = True
            if (has_variants) and (not has_non_ref):
                return DrugSensitivity.SENSITIVE

            return DrugSensitivity.UNDETERMINED

        def get_piperaquine_sensitivity(dr_report):
        # Loci utilized: PF3D7_1408000 (plasmepsin 2) and PF3D7_1408100 (plasmepsin 3)
        # Codons: Amplification status of both genes
        # Workflow:
        # Step Genetic change Interpretation Classification
        # 1 Missing Missing Undetermined
        # 2 Heterozygous duplication Heterozygous mutant Undetermined
        # 3 Single copy Wild type Sensitive
        # 4 Multiple copies Mutant Resistant

            # Currently we can't determine this.
            # We need to get CNV calling working first.

            return DrugSensitivity.UNDETERMINED

        # Get the drug resistances:
        chloroquine = get_chloroquine_sensitivity("~{raw_drug_resistance_report}")
        pyrimethamine = get_pyrimethamine_sensitivity("~{raw_drug_resistance_report}")
        sulfadoxine = get_sulfadoxine_sensitivity("~{raw_drug_resistance_report}")
        mefloquine = get_mefloquine_sensitivity("~{raw_drug_resistance_report}")
        artemisinin = get_artemisinin_sensitivity("~{raw_drug_resistance_report}")
        piperaquine = get_piperaquine_sensitivity("~{raw_drug_resistance_report}")

        with open("~{outfile_name}", 'w') as f:
            f.write(f"#~{prefix} Drug Resistances:\n")
            f.write(f"Chloroquine: {chloroquine.name}\n")
            f.write(f"Pyrimethamine: {pyrimethamine.name}\n")
            f.write(f"Sulfadoxine: {sulfadoxine.name}\n")
            f.write(f"Mefloquine: {mefloquine.name}\n")
            f.write(f"Artemisinin: {artemisinin.name}\n")
            f.write(f"Piperaquine: {piperaquine.name}\n")
            f.write("\n")

        # Write our resistance status for each drug:
        with open("chloroquine_status.txt", 'w') as f:
            f.write(f"{chnoroquine.name}\n")
        with open("pyrimethamine_status.txt", 'w') as f:
            f.write(f"{pyrimethamine.name}\n")
        with open("sulfadoxine_status.txt", 'w') as f:
            f.write(f"{sulfadoxine.name}\n")
        with open("mefloquine_status.txt", 'w') as f:
            f.write(f"{mefloquine.name}\n")
        with open("artemisinin_status.txt", 'w') as f:
            f.write(f"{artemisinin.name}\n")
        with open("piperaquine_status.txt", 'w') as f:
            f.write(f"{piperaquine.name}\n")

        CODE
    >>>

    output {
        File resistance_summary = "~{outfile_name}"

        String chloroquine_status = read_string("chloroquine_status.txt")
        String pyrimethamine_status = read_string("pyrimethamine_status.txt")
        String sulfadoxine_status = read_string("sulfadoxine_status.txt")
        String mefloquine_status = read_string("mefloquine_status.txt")
        String artemisinin_status = read_string("artemisinin_status.txt")
        String piperaquine_status = read_string("piperaquine_status.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

