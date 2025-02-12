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
        vcf_index: "Index of the VCF file to process"

        snpeff_db: "SnpEff database for functional annotation"
        drug_resistance_list: "List of drug resistance markers for which to search"
        ref_map_file:    "Table indicating reference sequence, auxillary file locations, and metadata."

        dir_prefix: "Prefix for output directory"
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs.  If no bucket is given, outputs will not be finalized and instead will remain in their native execution location."

        do_functional_annotation: "Whether to perform functional annotation"
    }

    input {
        File vcf
        File vcf_index

        File snpeff_db
        File drug_resistance_list
        File ref_map_file

        String dir_prefix
        String? gcs_out_root_dir

        Boolean do_functional_annotation = true
    }

    # Read ref map into map data type so we can access its fields:
    Map[String, String] ref_map = read_map(ref_map_file)

    if (do_functional_annotation) {
        call FUNK.FunctionallyAnnotateVariants { input: vcf = vcf, snpeff_db = snpeff_db }
    }

    call CallDrugResistanceMutations {
        input:
            vcf = select_first([FunctionallyAnnotateVariants.annotated_vcf, vcf]),
            vcf_index = select_first([FunctionallyAnnotateVariants.annotated_vcf_index, vcf_index]),
            drug_resistance_list = drug_resistance_list,
            genes_gff = ref_map["gff"]
    }

    call CreateDrugResistanceSummary {
        input:
            raw_drug_resistance_info = CallDrugResistanceMutations.report,
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

        String predicted_drug_status_chloroquine   = CreateDrugResistanceSummary.predicted_chloroquine_status
        String predicted_drug_status_pyrimethamine = CreateDrugResistanceSummary.predicted_pyrimethamine_status
        String predicted_drug_status_sulfadoxine   = CreateDrugResistanceSummary.predicted_sulfadoxine_status
        String predicted_drug_status_mefloquine    = CreateDrugResistanceSummary.predicted_mefloquine_status
        String predicted_drug_status_artemisinin   = CreateDrugResistanceSummary.predicted_artemisinin_status
        String predicted_drug_status_piperaquine   = CreateDrugResistanceSummary.predicted_piperaquine_status

        File? annotated_vcf = final_annotated_vcf
        File? annotated_vcf_index = final_annotated_vcf_index
        File? snpEff_summary = final_snpeff_summary
        File? snpEff_genes = final_snpeff_genes
    }
}

task CallDrugResistanceMutations {
    input {
        File vcf
        File vcf_index
        File drug_resistance_list
        File genes_gff

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size([vcf, drug_resistance_list], "GB"))
    String prefix = basename(basename(vcf, ".gz"), ".vcf")
    String out_file_name = prefix + ".raw_drug_resistance_info.txt"

    command <<<
        set -x

        python3 <<CODE

        import gzip
        import os
        import multiprocessing
        import pysam

        from collections import defaultdict
        from tqdm import tqdm

        ################################################################################
        # Define some helpers here:
        ################################################################################

        class Open(object):
            """Wrapper for standard and gzipped files.

            Inspired by: https://stackoverflow.com/a/23796737"""

            def __init__(self, *args, **kwds):
                self.args = args
                self.kwds = kwds

                _, self.ext = os.path.splitext(self.args[0])

            def __enter__(self):

                if self.ext == '.gz':
                    self.file_obj = gzip.open(*self.args, **self.kwds)
                else:
                    self.file_obj = open(*self.args, **self.kwds)

                # return actual file object so we don't have to worry
                # about proxying
                return self.file_obj

            def __exit__(self, *args):
                self.file_obj.close()


        def get_num_lines(file_path):

            def get_blocks(files, size=65536):
                while True:
                    b = files.read(size)
                    if not b: break
                    yield b

            with Open(file_path, "rt", encoding="utf-8", errors='ignore') as f:
                return sum(bl.count("\n") for bl in get_blocks(f))


        def get_genes_from_gff_file(gff_file):
            # TSV Fields:
            # seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
            # source - name of the program that generated this feature, or the data source (database or project name)
            # feature - feature type name, e.g. Gene, Variation, Similarity
            # start - Start position* of the feature, with sequence numbering starting at 1.
            # end - End position* of the feature, with sequence numbering starting at 1.
            # score - A floating point value.
            # strand - defined as + (forward) or - (reverse).
            # frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
            # attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
            # *- Both, the start and end position are included. For example, setting start-end to 1-2 describes two bases, the first and second base in the sequence.

            gff_gene_dict = defaultdict(list)

            num_lines = get_num_lines(gff_file)
            with Open(gff_file, 'r') as gff:

                for line in tqdm(gff, desc="Loading Genes from GFF", total=num_lines):
                    # Skip headers:
                    if line.strip().startswith("#"):
                        continue

                    gff_fields = line.strip().split("\t")
                    seq_name = gff_fields[0]
                    ftype = gff_fields[2]
                    start = int(gff_fields[3])
                    end = int(gff_fields[4])

                    # Do our best to find genes.
                    # This is hard because teh spec allows anything...
                    if "gene" in ftype.lower():
                        # We have the contig, start, and end.
                        # If this gene has a name, we should get it.
                        # otherwise, we'll use the ID:
                        annotations = {raw.split("=")[0].lower():raw.split("=")[1] for raw in gff_fields[8].split(";")}
                        gene_id = annotations["id"]
                        if "name" in annotations:
                            gene_name = annotations["name"]
                        else:
                            gene_name = gene_id

                        gff_gene_dict[seq_name].append((gene_name, start, end, gene_id))

                # Make sure we sort by both contig and the positions of the genes:
                return {k:sorted(v, key=lambda x: (x[1], x[2])) for k, v in gff_gene_dict.items()}

        ################################################################################
        # Define some helpers here:
        ################################################################################

        single_sample_VCF = f"~{vcf}"
        drug_res_info_out_file = f"~{out_file_name}"

        gff_file = f"~{genes_gff}"
        drug_resistance_list = f"~{drug_resistance_list}"

        # Read in GFF file:
        pos_gene_dict = get_genes_from_gff_file(gff_file)
        gene_id_pos_dict = dict()
        for contig, gene_list in pos_gene_dict.items():
            for gene_name, start, end, gene_id in gene_list:
                if gene_id in gene_id_pos_dict:
                    raise RuntimeError(f"Duplicate Gene ID: {gene_id}")
                else:
                    gene_id_pos_dict[gene_id] = (contig, start, end)

        # Read in drug resistance list:
        drug_res_info = []
        drug_res_info_by_gene = defaultdict(set)
        drug_res_info_by_gene_id = defaultdict(set)

        num_lines = get_num_lines(drug_resistance_list)
        with Open(drug_resistance_list, 'r') as f:
            for line in tqdm(f, desc="Loading drug resistance markers", total=num_lines):
                gene, gene_id, aa_change = line.strip().split("\t")
                drug_res_info.append((gene, gene_id, aa_change))
                drug_res_info_by_gene[gene].add(aa_change)
                drug_res_info_by_gene_id[gene_id].add(aa_change)

        # Extract variants of interest from the VCF:
        regions = [gene_id_pos_dict[gene_id] for gene_id in drug_res_info_by_gene_id.keys()]

        raw_variants_of_interest = []
        with pysam.VariantFile(single_sample_VCF, 'r') as vcf:
            for contig, start, end in tqdm(regions,desc="Extracting variants from drug res loci"):
                for v in vcf.fetch(contig, start, end):
                    raw_variants_of_interest.append(v)

        # Now filter the variants of interest to the specific genes and amino acid changes in our list:
        variants_of_interest = []
        for v in raw_variants_of_interest:
            # Get the VEP annotations:
            for ann in v.info['ANN']:
                # Get the gene id and protein change string from this annotation:
                # We know they are always fields 3 and 10, respectively:
                g_id = ann.split("|")[3]
                p_change = ann.split("|")[10]
                if p_change in drug_res_info_by_gene_id[g_id]:
                    variants_of_interest.append(v)
                    break

        print(f"Raw variants of interest extracted: {len(raw_variants_of_interest)}")
        print(f"Variants of interest extracted: {len(variants_of_interest)}")

        # Now for each drug res marker, go through and check if that marker is present.
        # If the marker is present, extract the genotype and mark it as het / hom var
        # If the marker is absent, since this is a single-sample VCF, mark it as absent
        #   Strictly speaking this could be because we had no coverage at this site, but
        #   for single-sample VCFs it doesn't matter.
        drug_res_summary_info = []
        num_found = 0
        for gene, gene_id, aa_change in tqdm(drug_res_info, desc="Checking variants for drug res markers"):
            status = "absent"
            for v in variants_of_interest:
                for ann in v.info['ANN']:
                    g_id = ann.split("|")[3]
                    p_change = ann.split("|")[10]
                    if gene_id == g_id and aa_change == p_change:
                        # We have the annotation we're looking for.
                        # Now get the genotype!
                        gt_a, gt_b = v.samples.values()[0]["GT"]

                        # Check for hom:
                        if gt_a == gt_b:
                            if gt_a != 0:
                                # hom var - update status:
                                status = "hom_var"
                        else:
                            # het of some kind:
                            status = "het"

                        num_found += 1
                        break

            drug_res_summary_info.append((gene, gene_id, aa_change, status))
        print(f"Num drug resistance markers found: {num_found}")
        print("")
        print("Drug resistance marker table:")
        for summary_info in drug_res_summary_info:
            print("\t".join(summary_info))

        with open(drug_res_info_out_file, 'w') as f:
            for summary_info in tqdm(drug_res_summary_info, desc="Writing drug resistance summary file"):
                f.write("\t".join(summary_info))
                f.write("\n")
        CODE
    >>>

    output {
        File report = "~{out_file_name}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
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
        File raw_drug_resistance_info
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    String outfile_name = "~{prefix}.drug_resistance_summary.txt"

    Int disk_size = 1 + 4*ceil(size(raw_drug_resistance_info, "GB"))

    command <<<
        python3 <<CODE

        import re
        from enum import Enum

        # Debugging code:
        OUT_PATH_PREFIX = ""

        # Set up some functions to determine drug sensitivity:
        pchange_regex = re.compile(r"""p\.([A-z]+)([0-9]+)([A-z]+)""")
        AA_3_2 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

        DrugSensitivity = Enum("DrugSensitivity", ["UNDETERMINED", "SENSITIVE", "RESISTANT"])

        def parse_pchange(pchange):
            old, pos, new = pchange_regex.match(pchange).groups()
            return AA_3_2[old.upper()], int(pos), AA_3_2[new.upper()]

        def line_is_drug_locus(line, gene_name, gene_id):
            line_re = re.compile(f"^{gene_name}\s+{gene_id}")
            return line_re.match(line) is not None

        def get_chloroquine_sensitivity(dr_info_file):
        # Locus utilized: PF3D7_0709000 (crt)
        # Codon: 76
        # Workflow:
        # Step Genetic change Interpretation Classification
        # 1 76 K/T heterozygote Heterozygous mutant Undetermined
        # 2 76 missing Missing Undetermined
        # 3 K76 Wild type Sensitive
        # 4 76T Mutant Resistant
        # 5 otherwise Unknown mutant Undetermined
            with open(dr_info_file, 'r') as f:
                for line in f:
                    # We only care about this locus for this drug:
                    if line_is_drug_locus(line, "pfcrt", "PF3D7_0709000"):
                        gene, locus, pchange, marker = re.split("\s", line.strip())
                        old, pos, new = parse_pchange(pchange)
                        if pos == 76:
                            if (old == "L") and (new == "T") and (marker == "absent"):
                                return DrugSensitivity.SENSITIVE
                            elif (new == "T") and ((marker == "present") or (marker == "hom_var")):
                                return DrugSensitivity.RESISTANT

            return DrugSensitivity.UNDETERMINED

        def get_pyrimethamine_sensitivity(dr_info_file):
        # Locus utilized: PF3D7_0417200 (dhfr)
        # Codon: 108
        # Workflow:
        # Step Genetic change Interpretation Classification
        # 1 108 S/N heterozygote Heterozygous mutant Undetermined
        # 2 108 missing Missing Undetermined
        # 3 S108 Wild type Sensitive
        # 4 108N Mutant Resistant
        # 5 otherwise Unknown mutant Undetermined

            with open(dr_info_file, 'r') as f:
                for line in f:
                    # We only care about this locus for this drug:
                    if line_is_drug_locus(line, "pfdhfr", "PF3D7_0417200"):
                        gene, locus, pchange, marker = re.split("\s", line.strip())
                        old, pos, new = parse_pchange(pchange)
                        if pos == 108:
                            if (old == "S") and (new == "N") and (marker == "absent"):
                                return DrugSensitivity.SENSITIVE
                            elif (old == "S") and (new == "N") and ((marker == "present") or (marker == "hom_var")):
                                return DrugSensitivity.RESISTANT
                            elif (new == "N") and ((marker == "present") or (marker == "hom_var")):
                                return DrugSensitivity.RESISTANT
                            elif marker == "absent":
                                return DrugSensitivity.SENSITIVE

            return DrugSensitivity.UNDETERMINED


        def get_sulfadoxine_sensitivity(dr_info_file):
        # Locus utilized: PF3D7_0810800 (dhps)
        # Codon: 437
        # Workflow:
        # Step Genetic change Interpretation Classification
        # 1 437 A/G heterozygote Heterozygous mutant Undetermined
        # 2 437 missing Missing Undetermined
        # 3 A437 Wild type Sensitive
        # 4 437G Mutant Resistant
        # 5 otherwise Unknown mutant Undetermined

            with open(dr_info_file, 'r') as f:
                for line in f:
                    # We only care about this locus for this drug:
                    if line_is_drug_locus(line, "pfdhps", "PF3D7_0810800"):
                        gene, locus, pchange, marker = re.split("\s", line.strip())
                        old, pos, new = parse_pchange(pchange)
                        if pos == 437:
                            if (old == "A") and (new == "G") and (marker == "absent"):
                                return DrugSensitivity.SENSITIVE
                            elif (old == "A") and (new == "G") and ((marker == "present") or (marker == "hom_var")):
                                return DrugSensitivity.RESISTANT
                            elif (new == "G") and ((marker == "present") or (marker == "hom_var")):
                                return DrugSensitivity.RESISTANT
                            elif marker == "absent":
                                return DrugSensitivity.SENSITIVE

            return DrugSensitivity.UNDETERMINED

        def get_mefloquine_sensitivity(dr_info_file):
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

        def get_artemisinin_sensitivity(dr_info_file):
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

            has_variants = False
            has_non_ref = False
            with open(dr_info_file, 'r') as f:
                for line in f:
                    # We only care about this locus for this drug:
                    if line_is_drug_locus(line, "pfkelch13", "PF3D7_1343700"):
                        gene, locus, pchange, marker = re.split("\s", line.strip())
                        old, pos, new = parse_pchange(pchange)
                        if 349 <= pos <= 726:
                            if (old != new) and ((marker == "present") or (marker == "hom_var")):
                                return DrugSensitivity.RESISTANT
                            elif (new == "S") and ((marker == "present") or (marker == "hom_var")):
                                return DrugSensitivity.SENSITIVE
                            elif marker == "present":
                                has_non_ref = True
                            has_variants = True

            if (has_variants) and (not has_non_ref):
                return DrugSensitivity.SENSITIVE

            return DrugSensitivity.UNDETERMINED

        def get_piperaquine_sensitivity(dr_info_file):
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
        print("Determining drug resistances...")
        chloroquine = get_chloroquine_sensitivity("~{raw_drug_resistance_info}")
        pyrimethamine = get_pyrimethamine_sensitivity("~{raw_drug_resistance_info}")
        sulfadoxine = get_sulfadoxine_sensitivity("~{raw_drug_resistance_info}")
        mefloquine = get_mefloquine_sensitivity("~{raw_drug_resistance_info}")
        artemisinin = get_artemisinin_sensitivity("~{raw_drug_resistance_info}")
        piperaquine = get_piperaquine_sensitivity("~{raw_drug_resistance_info}")

        print("Drug resistances:")
        print(f"Chloroquine: {chloroquine.name}")
        print(f"Pyrimethamine: {pyrimethamine.name}")
        print(f"Sulfadoxine: {sulfadoxine.name}")
        print(f"Mefloquine: {mefloquine.name}")
        print(f"Artemisinin: {artemisinin.name}")
        print(f"Piperaquine: {piperaquine.name}")
        print()

        print("Writing out summary file...")
        with open(f"~{outfile_name}", 'w') as f:
            f.write(f"#~{prefix} Drug Resistances:\n")
            f.write(f"Chloroquine: {chloroquine.name}\n")
            f.write(f"Pyrimethamine: {pyrimethamine.name}\n")
            f.write(f"Sulfadoxine: {sulfadoxine.name}\n")
            f.write(f"Mefloquine: {mefloquine.name}\n")
            f.write(f"Artemisinin: {artemisinin.name}\n")
            f.write(f"Piperaquine: {piperaquine.name}\n")
            f.write("\n")
        print()

        # Write our resistance status for each drug:
        print("Writing individual status files...")
        with open(f"{OUT_PATH_PREFIX}chloroquine_status.txt", 'w') as f:
            f.write(f"{chloroquine.name}\n")
        with open(f"{OUT_PATH_PREFIX}pyrimethamine_status.txt", 'w') as f:
            f.write(f"{pyrimethamine.name}\n")
        with open(f"{OUT_PATH_PREFIX}sulfadoxine_status.txt", 'w') as f:
            f.write(f"{sulfadoxine.name}\n")
        with open(f"{OUT_PATH_PREFIX}mefloquine_status.txt", 'w') as f:
            f.write(f"{mefloquine.name}\n")
        with open(f"{OUT_PATH_PREFIX}artemisinin_status.txt", 'w') as f:
            f.write(f"{artemisinin.name}\n")
        with open(f"{OUT_PATH_PREFIX}piperaquine_status.txt", 'w') as f:
            f.write(f"{piperaquine.name}\n")


        CODE
    >>>

    output {
        File resistance_summary = "~{outfile_name}"

        String predicted_chloroquine_status = read_string("chloroquine_status.txt")
        String predicted_pyrimethamine_status = read_string("pyrimethamine_status.txt")
        String predicted_sulfadoxine_status = read_string("sulfadoxine_status.txt")
        String predicted_mefloquine_status = read_string("mefloquine_status.txt")
        String predicted_artemisinin_status = read_string("artemisinin_status.txt")
        String predicted_piperaquine_status = read_string("piperaquine_status.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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

