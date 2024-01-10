version 1.0

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PfalciparumDrugResistanceSummary {
    meta {
        desciption: "Create a drug resistance report based on the given raw drug resistance loci report."
    }

    input {
        File raw_drug_resistance_report

        String participant_name

        String gcs_out_root_dir
    }

    parameter_meta {
        raw_drug_resistance_report: "File containing a raw drug resistance report to use to determine drug resistance."
        participant_name:    "Participant (or sample) name for the given bam file."
        gcs_out_root_dir:    "Output folder into which to place the results of this workflow."
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PfalciparumDrugResistanceSummary/~{participant_name}"

    call CreateDrugResistanceSummary as CreateDrugResistanceSummary {
        input:
            raw_drug_resistance_report = raw_drug_resistance_report,
            prefix = participant_name
    }

    call FF.FinalizeToFile as FinalizeDrugResistanceSummary { input: outdir = outdir, file = CreateDrugResistanceSummary.resistance_summary }

    output {
        File drug_resistance_summary = FinalizeDrugResistanceSummary.gcs_path
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

        CODE
    >>>
    output {
        File resistance_summary = "~{outfile_name}"
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
