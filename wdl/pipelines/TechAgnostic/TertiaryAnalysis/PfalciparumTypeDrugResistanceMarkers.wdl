version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/TertiaryAnalysis/FunctionalAnnotation.wdl" as FUNK

workflow PfalciparumTypeDrugResistanceMarkers {

    meta {
        description: "Workflow to generate a report of drug resistance markers"
    }
    parameter_meta {
        vcf: "VCF file to process"
        vcf_index: "Index of the VCF file to process"

        gvcf: "GVCF file to process to use for coverage analysis of drug resistance markers"
        gvcf_index: "Index of the GVCF file to process"

        snpeff_db: "SnpEff database for functional annotation"
        snpeff_db_identifier: "Identifier for the SnpEff database to use"
        drug_resistance_list: "List of drug resistance markers for which to search"
        ref_map_file:    "Table indicating reference sequence, auxillary file locations, and metadata."

        do_functional_annotation: "Whether to perform functional annotation"
    }

    input {
        File vcf
        File vcf_index

        File gvcf
        File gvcf_index

        File snpeff_db
        String snpeff_db_identifier
        File drug_resistance_list
        File ref_map_file

        Boolean do_functional_annotation = true
    }

    # Read ref map into map data type so we can access its fields:
    Map[String, String] ref_map = read_map(ref_map_file)

    if (do_functional_annotation) {
        call FUNK.FunctionallyAnnotateVariants { input: vcf = vcf, snpeff_db = snpeff_db, snpeff_db_identifier = snpeff_db_identifier }
    }

    call CallDrugResistanceMutations {
        input:
            vcf = select_first([FunctionallyAnnotateVariants.annotated_vcf, vcf]),
            vcf_index = select_first([FunctionallyAnnotateVariants.annotated_vcf_index, vcf_index]),
            gvcf = gvcf,
            gvcf_index = gvcf_index,
            drug_resistance_list = drug_resistance_list,
            genes_gff = ref_map["gff"]
    }

    call CreateDrugResistanceSummary {
        input:
            raw_drug_resistance_info = CallDrugResistanceMutations.report,
            prefix = basename(basename(vcf, ".gz"), ".vcf")
    }

    output {
        File drug_resistance_summary = CreateDrugResistanceSummary.resistance_summary
        File drug_resistance_markers = CallDrugResistanceMutations.report

        String predicted_drug_status_chloroquine   = CreateDrugResistanceSummary.predicted_chloroquine_status
        String predicted_drug_status_pyrimethamine = CreateDrugResistanceSummary.predicted_pyrimethamine_status
        String predicted_drug_status_sulfadoxine   = CreateDrugResistanceSummary.predicted_sulfadoxine_status
        String predicted_drug_status_mefloquine    = CreateDrugResistanceSummary.predicted_mefloquine_status
        String predicted_drug_status_artemisinin   = CreateDrugResistanceSummary.predicted_artemisinin_status
        String predicted_drug_status_piperaquine   = CreateDrugResistanceSummary.predicted_piperaquine_status

        File annotated_vcf = select_first([FunctionallyAnnotateVariants.annotated_vcf, vcf])
        File annotated_vcf_index = select_first([FunctionallyAnnotateVariants.annotated_vcf_index, vcf_index])

        File? snpEff_summary = FunctionallyAnnotateVariants.snpEff_summary
        File? snpEff_genes = FunctionallyAnnotateVariants.snpEff_genes

        # Pull out the drug resistance markers from the raw drug resistance report:
        
    }
}

task CallDrugResistanceMutations {
    input {
        File vcf
        File vcf_index

        File gvcf
        File gvcf_index

        File drug_resistance_list
        File genes_gff

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size([vcf, gvcf, drug_resistance_list], "GB"))
    String prefix = basename(basename(vcf, ".gz"), ".vcf")
    String out_file_name = prefix + ".raw_drug_resistance_info.txt"

    command <<<
        set -x

        python3 <<CODE

        import gzip
        import os
        import multiprocessing
        import re
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

        def get_full_gene_info_from_gff_file(gff_file):
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

            gene_cds_transcript_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

            num_lines = get_num_lines(gff_file)
            with Open(gff_file, 'r') as gff:

                for line in tqdm(gff, desc="Loading Gene CDS from GFF", total=num_lines):
                    # Skip headers:
                    if line.strip().startswith("#"):
                        continue

                    gff_fields = line.strip().split("\t")
                    contig = gff_fields[0]
                    ftype = gff_fields[2]
                    start = int(gff_fields[3])
                    end = int(gff_fields[4])

                    # Do our best to find genes.
                    # This is hard because the spec allows anything...
                    if "CDS" == ftype:
                        # We have the contig, start, and end.
                        # If this gene has a name, we should get it.
                        # otherwise, we'll use the ID:
                        annotations = {raw.split("=")[0].lower():raw.split("=")[1] for raw in gff_fields[8].split(";")}
                        gene_id = annotations["gene_id"]
                        transcript_id = annotations['parent']
                        
                        gene_cds_transcript_dict[contig][gene_id][transcript_id].append((start, end))

            # Make sure we sort by both contig and the positions of the genes:
            sorted_gene_cds_transcript_dict = defaultdict(dict)
            for contig in sorted(gene_cds_transcript_dict.keys()):
                gene_dict = gene_cds_transcript_dict[contig]
                for gene in gene_dict.keys():
                    transcript_dict = gene_dict[gene]
                    sorted_gene_cds_transcript_dict[contig][gene] = {k:sorted(v, key=lambda x: x[0]) for k, v in transcript_dict.items()}
                
            return sorted_gene_cds_transcript_dict


        def convert_protein_change_strings_to_genomic_positions(drug_res_info_by_gene_id, gene_cds_transcript_dict, pos_gene_dict, gene_id_pos_dict):

            p_change_regex = re.compile('p\.([A-z]+)([0-9]+)([A-z]+)')

            # drug_res_info_by_gene_id: gene_id -> set(protein_change_strings)
            # gene_cds_transcript_dict: contig -> {gene_id -> list(tuple(cds_start, cds_end))}
            # pos_gene_dict: contig -> list(tuple(gene_name, gene_start, gene_end, gene_id, direction (+/-) ))
            # gene_id_pos_dict: gene_id -> tuple(contig, gene_start, gene_end, direction (+/-) )

            # Now get genomic positions for each protein change we're looking for.
            # We need these positions so we can check for coverage.
            drug_res_info_by_gene_id_with_genome_pos = defaultdict(lambda: defaultdict(tuple))
            for gene_id, protein_changes in drug_res_info_by_gene_id.items():

                # Get the contig and direction for this gene:
                contig, _, _, direction = gene_id_pos_dict[gene_id]

                # Get the transcripts for this gene:
                transcript_dict = gene_cds_transcript_dict[contig][gene_id]

                for p in protein_changes:
                    ref_aa, aa_pos, alt_aa = p_change_regex.match(p).groups()
                    aa_pos = int(aa_pos)

                    # Check all transcripts in this gene for this protein change:
                    for transcript, cds_regions in transcript_dict.items():

                        # Get the total number of base pairs in the coding sequence:
                        total_cds_length_bp = sum([cds[1] - cds[0] + 1 for cds in cds_regions])

                        # For each CDS region, compute the protein coordinates:
                        # Note: We know this is divisible by 3, so we can convert to int without worrying.
                        total_num_aas_in_transcript = int(total_cds_length_bp/3)

                        if aa_pos > total_num_aas_in_transcript:
                            # This transcript doesn't work for this protein position,
                            # so we need to go to the next one:
                            continue

                        # We could be clever and do some math, or we could special case the reverse strand.
                        # We do the latter because it's easier to implement:
                        if direction == "+":

                            # Get the start pos of the transcript on the genome:
                            transcript_start = cds_regions[0][0]

                            # Get the end position of the transcript on the genome:
                            transcript_end = cds_regions[-1][1]

                            # Get the length of hte transcript in base pairs:
                            transcript_length_bp = sum([e - s for (s,e) in cds_regions])

                            last_cds_end_pos_bp = 0
                            for i, (cds_start, cds_end) in enumerate(cds_regions):
                                cds_bp_start = last_cds_end_pos_bp + 1
                                cds_bp_end = cds_end - cds_start + last_cds_end_pos_bp

                                cds_aa_start = cds_bp_start / 3
                                cds_aa_end = cds_bp_end / 3

                                if cds_aa_start <= aa_pos <= cds_aa_end:

                                    # now convert the AA Pos into a genomic position
                                    # based on the start / end position of this CDS:
                                    aa_genomic_pos = int(cds_start + (aa_pos - cds_aa_start) * 3)

                                    # Add the genomic interval of the amino acid to the new dictionary:
                                    # Add only 2 because of inclusive positions:
                                    drug_res_info_by_gene_id_with_genome_pos[gene_id][p] = (contig, aa_genomic_pos, aa_genomic_pos+2)

                                    # We're done, let's go to the next protein change:
                                    break

                                last_cds_end_pos_bp = cds_bp_end
                        else:                    
                            # Get the start pos of the transcript on the genome:
                            transcript_start = cds_regions[-1][1]

                            # Get the end position of the transcript on the genome:
                            transcript_end = cds_regions[0][0]

                            # Get the length of hte transcript in base pairs:
                            transcript_length_bp = sum([e - s for (s,e) in cds_regions])
                            
                            # Add one here to account for the -1 below:
                            last_cds_start_pos = 0
                            for i, (cds_start, cds_end) in enumerate(cds_regions[::-1]):
                                cds_bp_end = last_cds_start_pos + 1
                                cds_bp_start = cds_bp_end + (cds_end-cds_start) - 1

                                cds_aa_start = (cds_bp_start) / 3
                                cds_aa_end = (cds_bp_end) / 3
                                
                                # Because end is before start, we have to check the opposite
                                # of the relationship than we normally would:
                                if cds_aa_end <= aa_pos <= cds_aa_start:
                                    
                                    # now convert the AA Pos into a genomic position
                                    # based on the start / end position of this CDS:
                                    aa_genomic_pos = int(cds_end - (cds_aa_end - aa_pos) * 3)

                                    # Add the genomic interval of the amino acid to the new dictionary:
                                    # Add only 2 because of inclusive positions:
                                    drug_res_info_by_gene_id_with_genome_pos[gene_id][p] = (contig, aa_genomic_pos-2, aa_genomic_pos)

                                    break
                                
                                last_cds_start_pos = cds_bp_start


            return drug_res_info_by_gene_id_with_genome_pos


        def check_coverage_of_ref_protein_changes(protein_change_info_dict, gvcf_file, min_GQ=20):
            
            gene_protein_change_presence_dict = defaultdict(lambda: defaultdict(str))
            
            for gene_id, prot_changes in protein_change_info_dict.items():
                for prot_change, (contig, start, end) in prot_changes.items():
                    with pysam.VariantFile(gvcf_file, 'r') as vcf:
                        sample = next(iter(vcf.header.samples))
                        gqs = []
                        for v in vcf.fetch(contig, start, end):
                            # Get the GQ and add it to our list of GQs for the regions:
                            # Note: We do not have to check for "<NON_REF>" alleles here.
                            #       We are simply interested in whether we saw any data over
                            #       this locus.
                            gqs.append(v.samples[sample]['GQ'])
                            
                        # Aggregate genotype qualities.
                        # A simple average is probably sufficient.
                        if len(gqs) > 1:
                            final_GQ = sum(gqs)/len(gqs)
                        elif len(gqs) == 1:
                            final_GQ = gqs[0]
                            
                        # Does our final GQ meet or exceed the threshold for "good" data?
                        if final_GQ >= min_GQ:
                            pchange_status = "LOCUS_COVERED"
                        else:
                            pchange_status = "MISSING"
                            
                    gene_protein_change_presence_dict[gene_id][prot_change] = pchange_status
                
            return gene_protein_change_presence_dict

        ################################################################################
        # Define some helpers here:
        ################################################################################

        single_sample_VCF = f"~{vcf}"
        single_sample_GVCF = f"~{gvcf}"
        drug_res_info_out_file = f"~{out_file_name}"

        gff_file = f"~{genes_gff}"
        drug_resistance_list = f"~{drug_resistance_list}"

        # Read in GFF file:
        pos_gene_dict = get_genes_from_gff_file(gff_file)
        gene_cds_transcript_dict = get_full_gene_info_from_gff_file(gff_file)
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
        
        AMBIGUOUS_DATA_LABEL = "absent"

        drug_res_summary_info = []
        num_found = 0
        for gene, gene_id, aa_change in tqdm(drug_res_info, desc="Checking variants for drug res markers"):
            status = AMBIGUOUS_DATA_LABEL
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

        print("Checking for missing or hom_ref protein changes...")

        # Get the genomic coordinates of every protein change:
        drug_res_info_by_gene_id_with_genome_pos = convert_protein_change_strings_to_genomic_positions(drug_res_info_by_gene_id, gene_cds_transcript_dict, pos_gene_dict, gene_id_pos_dict)

        # Get the presence of each protein change using the GVCF for genotype quality measurements:
        gene_id_prot_change_presence_dict = check_coverage_of_ref_protein_changes(drug_res_info_by_gene_id_with_genome_pos, single_sample_GVCF, min_GQ=20)

        # Now we can review our drug_res_summary_info table and resolve anything we see as absent:
        final_drug_res_summary_info = []
        for gene, gene_id, aa_change, status in drug_res_summary_info:
            if status == AMBIGUOUS_DATA_LABEL:
                if gene_id_prot_change_presence_dict[gene_id][aa_change] == "LOCUS_COVERED":
                    # We have coverage.
                    # This is a hom_ref variant:
                    status = "hom_ref"
                elif gene_id_prot_change_presence_dict[gene_id][aa_change] == "MISSING":
                    # We have no coverage!
                    # This variant is missing!
                    status = "missing"
                else:
                    raise RuntimeError(f"Unknown coverage status ({gene_id_prot_change_presence_dict[gene_id][aa_change]})!  This should never happen!")
                    
            final_drug_res_summary_info.append((gene, gene_id, aa_change, status))

        print("Drug resistance marker table:")
        for summary_info in final_drug_res_summary_info:
            print("\t".join(summary_info))

        with open(drug_res_info_out_file, 'w') as f:
            for summary_info in tqdm(final_drug_res_summary_info, desc="Writing drug resistance summary file"):
                f.write("\t".join(summary_info))
                f.write("\n")

        # Now write out the drug resistance markers to their own files:
        for gene, gene_id, aa_change, status in tqdm(final_drug_res_summary_info, desc="Writing drug resistance markers"):
            marker_file_name = f"{gene.lower()}.{gene_id}.{aa_change}.txt"
            print(f"Writing marker file: {marker_file_name} -> {status}")
            with open(f"{marker_file_name}", 'w') as f:
                f.write(f"{status}")
        CODE

        # Now we need to make sure all our output files exist so that they can be read from,
        # irrespective of whether the associated markers are present in the input file.
        touch pfFd.PF3D7_1318100.p.Asp193Tyr.txt
        touch pfaat1.PF3D7_0629500.p.Gln454Glu.txt
        touch pfaat1.PF3D7_0629500.p.Lys541Asn.txt
        touch pfaat1.PF3D7_0629500.p.Phe313Ser.txt
        touch pfaat1.PF3D7_0629500.p.Ser258Leu.txt
        touch pfap2-mu.PF3D7_1218300.p.Ile592Thr.txt
        touch pfarps10.PF3D7_1460900.p.Val127Met.txt
        touch pfatg18.PF3D7_1012900.p.Thr38Ile.txt
        touch pfcarl.PF3D7_0321900.p.Ile1139Lys.txt
        touch pfcarl.PF3D7_0321900.p.Leu830Val.txt
        touch pfcarl.PF3D7_0321900.p.Ser1076Asn.txt
        touch pfcarl.PF3D7_0321900.p.Ser1076Ile.txt
        touch pfcarl.PF3D7_0321900.p.Val1103Leu.txt
        touch pfcoronin.PF3D7_1251200.p.Arg100Lys.txt
        touch pfcoronin.PF3D7_1251200.p.Glu107Val.txt
        touch pfcoronin.PF3D7_1251200.p.Gly50Glu.txt
        touch pfcoronin.PF3D7_1251200.p.Pro76Ser.txt
        touch pfcrt.PF3D7_0709000.p.Asn75Glu.txt
        touch pfcrt.PF3D7_0709000.p.Cys101Phe.txt
        touch pfcrt.PF3D7_0709000.p.Cys72Ser.txt
        touch pfcrt.PF3D7_0709000.p.Gly353Val.txt
        touch pfcrt.PF3D7_0709000.p.His97Tyr.txt
        touch pfcrt.PF3D7_0709000.p.Lys76Thr.txt
        touch pfcrt.PF3D7_0709000.p.Met343Leu.txt
        touch pfcrt.PF3D7_0709000.p.Met74Ile.txt
        touch pfcrt.PF3D7_0709000.p.Phe145Ile.txt
        touch pfcrt.PF3D7_0709000.p.Ser350Arg.txt
        touch pfdhfr.PF3D7_0417200.p.Asn51Ile.txt
        touch pfdhfr.PF3D7_0417200.p.Cys50Arg.txt
        touch pfdhfr.PF3D7_0417200.p.Cys59Arg.txt
        touch pfdhfr.PF3D7_0417200.p.Ile164Lys.txt
        touch pfdhfr.PF3D7_0417200.p.Ser108Asn.txt
        touch pfdhps.PF3D7_0810800.p.Ala581Gly.txt
        touch pfdhps.PF3D7_0810800.p.Ala613Ser.txt
        touch pfdhps.PF3D7_0810800.p.Ala613Thr.txt
        touch pfdhps.PF3D7_0810800.p.Lys437Gly.txt
        touch pfdhps.PF3D7_0810800.p.Lys540Glu.txt
        touch pfdhps.PF3D7_0810800.p.Ser436Ala.txt
        touch pfexo.PF3D7_1362500.p.Glu415Gly.txt
        touch pfkelch13.PF3D7_1343700.p.Ala675Val.txt
        touch pfkelch13.PF3D7_1343700.p.Arg539Thr.txt
        touch pfkelch13.PF3D7_1343700.p.Arg561His.txt
        touch pfkelch13.PF3D7_1343700.p.Arg633Ile.txt
        touch pfkelch13.PF3D7_1343700.p.Asn458Tyr.txt
        touch pfkelch13.PF3D7_1343700.p.Cys580Tyr.txt
        touch pfkelch13.PF3D7_1343700.p.Ile543Thr.txt
        touch pfkelch13.PF3D7_1343700.p.Met476Ile.txt
        touch pfkelch13.PF3D7_1343700.p.Met579Ile.txt
        touch pfkelch13.PF3D7_1343700.p.Phe446Ile.txt
        touch pfkelch13.PF3D7_1343700.p.Phe553Leu.txt
        touch pfkelch13.PF3D7_1343700.p.Phe574Leu.txt
        touch pfkelch13.PF3D7_1343700.p.Phe673Ile.txt
        touch pfkelch13.PF3D7_1343700.p.Pro441Ile.txt
        touch pfkelch13.PF3D7_1343700.p.Pro553Leu.txt
        touch pfkelch13.PF3D7_1343700.p.Pro574Leu.txt
        touch pfkelch13.PF3D7_1343700.p.Tyr493His.txt
        touch pfkelch13.PF3D7_1343700.p.Val568Gly.txt
        touch pfmdr1.PF3D7_0523000.p.Asn1024Asp.txt
        touch pfmdr1.PF3D7_0523000.p.Asn86Tyr.txt
        touch pfmdr1.PF3D7_0523000.p.Asp1246Tyr.txt
        touch pfmdr1.PF3D7_0523000.p.Ser1034Cys.txt
        touch pfmdr1.PF3D7_0523000.p.Tyr184Phe.txt
        touch pfmdr2.PF3D7_1447900.p.Thr484Ile.txt
        touch pfubp1.PF3D7_0104300.p.Val3275Phe.txt
    >>>

    output {
        File report = "~{out_file_name}"

        # Pull out the drug resistance markers from the raw drug resistance report:
        # Note: We have to list them individually here because Terra can't handle programmatic output naming.
        String pfFd_Asp_193_Tyr =      read_string("pfFd.PF3D7_1318100.p.Asp193Tyr.txt")
        String pfaat1_Gln_454_Glu =    read_string("pfaat1.PF3D7_0629500.p.Gln454Glu.txt")
        String pfaat1_Lys_541_Asn =    read_string("pfaat1.PF3D7_0629500.p.Lys541Asn.txt")
        String pfaat1_Phe_313_Ser =    read_string("pfaat1.PF3D7_0629500.p.Phe313Ser.txt")
        String pfaat1_Ser_258_Leu =    read_string("pfaat1.PF3D7_0629500.p.Ser258Leu.txt")
        String pfap2_mu_Ile_592_Thr =  read_string("pfap2-mu.PF3D7_1218300.p.Ile592Thr.txt")
        String pfarps10_Val_127_Met =  read_string("pfarps10.PF3D7_1460900.p.Val127Met.txt")
        String pfatg18_Thr_38_Ile =    read_string("pfatg18.PF3D7_1012900.p.Thr38Ile.txt")
        String pfcarl_Ile_1139_Lys =   read_string("pfcarl.PF3D7_0321900.p.Ile1139Lys.txt")
        String pfcarl_Leu_830_Val =    read_string("pfcarl.PF3D7_0321900.p.Leu830Val.txt")
        String pfcarl_Ser_1076_Asn =   read_string("pfcarl.PF3D7_0321900.p.Ser1076Asn.txt")
        String pfcarl_Ser_1076_Ile =   read_string("pfcarl.PF3D7_0321900.p.Ser1076Ile.txt")
        String pfcarl_Val_1103_Leu =   read_string("pfcarl.PF3D7_0321900.p.Val1103Leu.txt")
        String pfcoronin_Arg_100_Lys = read_string("pfcoronin.PF3D7_1251200.p.Arg100Lys.txt")
        String pfcoronin_Glu_107_Val = read_string("pfcoronin.PF3D7_1251200.p.Glu107Val.txt")
        String pfcoronin_Gly_50_Glu =  read_string("pfcoronin.PF3D7_1251200.p.Gly50Glu.txt")
        String pfcoronin_Pro_76_Ser =  read_string("pfcoronin.PF3D7_1251200.p.Pro76Ser.txt")
        String pfcrt_Asn_75_Glu =      read_string("pfcrt.PF3D7_0709000.p.Asn75Glu.txt")
        String pfcrt_Cys_101_Phe =     read_string("pfcrt.PF3D7_0709000.p.Cys101Phe.txt")
        String pfcrt_Cys_72_Ser =      read_string("pfcrt.PF3D7_0709000.p.Cys72Ser.txt")
        String pfcrt_Gly_353_Val =     read_string("pfcrt.PF3D7_0709000.p.Gly353Val.txt")
        String pfcrt_His_97_Tyr =      read_string("pfcrt.PF3D7_0709000.p.His97Tyr.txt")
        String pfcrt_Lys_76_Thr =      read_string("pfcrt.PF3D7_0709000.p.Lys76Thr.txt")
        String pfcrt_Met_343_Leu =     read_string("pfcrt.PF3D7_0709000.p.Met343Leu.txt")
        String pfcrt_Met_74_Ile =      read_string("pfcrt.PF3D7_0709000.p.Met74Ile.txt")
        String pfcrt_Phe_145_Ile =     read_string("pfcrt.PF3D7_0709000.p.Phe145Ile.txt")
        String pfcrt_Ser_350_Arg =     read_string("pfcrt.PF3D7_0709000.p.Ser350Arg.txt")
        String pfdhfr_Asn_51_Ile =     read_string("pfdhfr.PF3D7_0417200.p.Asn51Ile.txt")
        String pfdhfr_Cys_50_Arg =     read_string("pfdhfr.PF3D7_0417200.p.Cys50Arg.txt")
        String pfdhfr_Cys_59_Arg =     read_string("pfdhfr.PF3D7_0417200.p.Cys59Arg.txt")
        String pfdhfr_Ile_164_Lys =    read_string("pfdhfr.PF3D7_0417200.p.Ile164Lys.txt")
        String pfdhfr_Ser_108_Asn =    read_string("pfdhfr.PF3D7_0417200.p.Ser108Asn.txt")
        String pfdhps_Ala_581_Gly =    read_string("pfdhps.PF3D7_0810800.p.Ala581Gly.txt")
        String pfdhps_Ala_613_Ser =    read_string("pfdhps.PF3D7_0810800.p.Ala613Ser.txt")
        String pfdhps_Ala_613_Thr =    read_string("pfdhps.PF3D7_0810800.p.Ala613Thr.txt")
        String pfdhps_Lys_437_Gly =    read_string("pfdhps.PF3D7_0810800.p.Lys437Gly.txt")
        String pfdhps_Lys_540_Glu =    read_string("pfdhps.PF3D7_0810800.p.Lys540Glu.txt")
        String pfdhps_Ser_436_Ala =    read_string("pfdhps.PF3D7_0810800.p.Ser436Ala.txt")
        String pfexo_Glu_415_Gly =     read_string("pfexo.PF3D7_1362500.p.Glu415Gly.txt")
        String pfkelch13_Ala_675_Val = read_string("pfkelch13.PF3D7_1343700.p.Ala675Val.txt")
        String pfkelch13_Arg_539_Thr = read_string("pfkelch13.PF3D7_1343700.p.Arg539Thr.txt")
        String pfkelch13_Arg_561_His = read_string("pfkelch13.PF3D7_1343700.p.Arg561His.txt")
        String pfkelch13_Arg_633_Ile = read_string("pfkelch13.PF3D7_1343700.p.Arg633Ile.txt")
        String pfkelch13_Asn_458_Tyr = read_string("pfkelch13.PF3D7_1343700.p.Asn458Tyr.txt")
        String pfkelch13_Cys_580_Tyr = read_string("pfkelch13.PF3D7_1343700.p.Cys580Tyr.txt")
        String pfkelch13_Ile_543_Thr = read_string("pfkelch13.PF3D7_1343700.p.Ile543Thr.txt")
        String pfkelch13_Met_476_Ile = read_string("pfkelch13.PF3D7_1343700.p.Met476Ile.txt")
        String pfkelch13_Met_579_Ile = read_string("pfkelch13.PF3D7_1343700.p.Met579Ile.txt")
        String pfkelch13_Phe_446_Ile = read_string("pfkelch13.PF3D7_1343700.p.Phe446Ile.txt")
        String pfkelch13_Phe_553_Leu = read_string("pfkelch13.PF3D7_1343700.p.Phe553Leu.txt")
        String pfkelch13_Phe_574_Leu = read_string("pfkelch13.PF3D7_1343700.p.Phe574Leu.txt")
        String pfkelch13_Phe_673_Ile = read_string("pfkelch13.PF3D7_1343700.p.Phe673Ile.txt")
        String pfkelch13_Pro_441_Ile = read_string("pfkelch13.PF3D7_1343700.p.Pro441Ile.txt")
        String pfkelch13_Pro_553_Leu = read_string("pfkelch13.PF3D7_1343700.p.Pro553Leu.txt")
        String pfkelch13_Pro_574_Leu = read_string("pfkelch13.PF3D7_1343700.p.Pro574Leu.txt")
        String pfkelch13_Tyr_493_His = read_string("pfkelch13.PF3D7_1343700.p.Tyr493His.txt")
        String pfkelch13_Val_568_Gly = read_string("pfkelch13.PF3D7_1343700.p.Val568Gly.txt")
        String pfmdr1_Asn_1024_Asp =   read_string("pfmdr1.PF3D7_0523000.p.Asn1024Asp.txt")
        String pfmdr1_Asn_86_Tyr =     read_string("pfmdr1.PF3D7_0523000.p.Asn86Tyr.txt")
        String pfmdr1_Asp_1246_Tyr =   read_string("pfmdr1.PF3D7_0523000.p.Asp1246Tyr.txt")
        String pfmdr1_Ser_1034_Cys =   read_string("pfmdr1.PF3D7_0523000.p.Ser1034Cys.txt")
        String pfmdr1_Tyr_184_Phe =    read_string("pfmdr1.PF3D7_0523000.p.Tyr184Phe.txt")
        String pfmdr2_Thr_484_Ile =    read_string("pfmdr2.PF3D7_1447900.p.Thr484Ile.txt")
        String pfubp1_Val_3275_Phe =   read_string("pfubp1.PF3D7_0104300.p.Val3275Phe.txt")
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
            line_re = re.compile(f"^{gene_name}\s+{gene_id}", re.IGNORECASE)
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

