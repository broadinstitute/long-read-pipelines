version 1.0

import "tasks/ttmars.wdl"
import "tasks/VariantUtils.wdl"

workflow evaluateSVsWithTTMars {
    meta {
        description: "Adapted from WDL by Jean Monlong. Evaluate structural variants against truth assemblies using TT-Mars (https://github.com/ChaissonLab/TT-Mars). If necessary, lift-over files can be made from the reference genome and the two assembled haplotypes. "
    }
    parameter_meta {
        referenceFile: "reference fasta (.fa or .fa.gz)"
        svsFile: "VCF with the SVs to evaluate (.vcf.gz). Must be indexed"
        hap1File: "truth assembly for haplotype 1 (.fa or .fa.gz)"
        hap2File: "truth assembly for haplotype 1 (.fa or .fa.gz)"
        trfFile: "BED file listing all the regions in the reference with simple repeats"
        centromereFile: "centromere location of the reference (see https://github.com/ChaissonLab/TT-Mars/blob/main/centromere_hg38.txt)"
        threads: "how many threads to use for the LRA alignment (when lifting over is needed)"
        nbXchr: "specify if sample is XY (nbXchr=1) or XX (nbXchr=2)"
        nonCovReg1File: "lift-over coverage of haplotype 1 on the reference. Provided at https://github.com/ChaissonLab/TT-Mars or made by liftover.sh"
        nonCovReg2File: "lift-over coverage of haplotype 2 on the reference. Provided at https://github.com/ChaissonLab/TT-Mars or made by liftover.sh"
        loPosAssem1File: "lift-over information of reference on haplotype 1. Provided at https://github.com/ChaissonLab/TT-Mars or made by liftover.sh"
        loPosAssem2File: "lift-over information of reference on haplotype 2. Provided at https://github.com/ChaissonLab/TT-Mars or made by liftover.sh"
        loPosAssem10File: "lift-over information of haplotype 1 on the reference. Provided at https://github.com/ChaissonLab/TT-Mars or made by liftover.sh"
        loPosAssem20File: "lift-over information of haplotype 2 on the reference. Provided at https://github.com/ChaissonLab/TT-Mars or made by liftover.sh"
        arguments: "optional additional arguments accepted by tt-mars, e.g. -v -g -p -s"
        contigs: "which contigs/chromosomes to evaluate"
    }

    input {
        File referenceFile
        File svsFile
        File hap1File
        File hap2File
        File trfFile
        File centromereFile
        Int threads=16
        Int nbXchr=2
        File? nonCovReg1File
        File? nonCovReg2File
        File? loPosAssem1File
        File? loPosAssem2File
        File? loPosAssem10File
        File? loPosAssem20File
        String? arguments
        Array[String] contigs
    }

    ## if any of those input files are missing, the assemblies have to be lifted to the reference sequence
    if(!defined(nonCovReg1File) || !defined(nonCovReg2File) || !defined(loPosAssem1File) || !defined(loPosAssem2File) || !defined(loPosAssem10File) || !defined(loPosAssem20File)){
        ## align each haplotype to the reference with LRA
        call ttmars.lra_t as lra1 {
            input:
            reference_file=referenceFile,
            hap_file=hap1File,
            threads=threads
        }
        call ttmars.lra_t as lra2 {
            input:
            reference_file=referenceFile,
            hap_file=hap2File,
            threads=threads
        }
        ## make liftover files for TTmars
        call ttmars.liftover_t {
            input:
            hap1_lra_file=lra1.bamOutput,
            hap2_lra_file=lra2.bamOutput
        }
    }

    File non_cov_reg_1_file = select_first([nonCovReg1File, liftover_t.non_cov_reg_1_file])
    File non_cov_reg_2_file = select_first([nonCovReg2File, liftover_t.non_cov_reg_2_file])
    File lo_pos_assem1_file = select_first([loPosAssem1File, liftover_t.lo_pos_assem1_file])
    File lo_pos_assem2_file = select_first([loPosAssem2File, liftover_t.lo_pos_assem2_file])
    File lo_pos_assem1_0_file = select_first([loPosAssem10File, liftover_t.lo_pos_assem1_0_file])
    File lo_pos_assem2_0_file = select_first([loPosAssem20File, liftover_t.lo_pos_assem2_0_file])

    ## call TT-Mars on each chromosome
    scatter(contig in contigs) {
        call VariantUtils.SubsetVCF { input: vcf_gz = svsFile, vcf_tbi = svsFile + ".tbi", locus = contig }

        call ttmars.ttmars_t {
            input:
            reference_file=referenceFile,
            svs_file=SubsetVCF.subset_vcf,
            hap1_file=hap1File,
            hap2_file=hap2File,
            centromere_file=centromereFile,
            trf_file=trfFile,
            non_cov_reg_1_file=non_cov_reg_1_file,
            non_cov_reg_2_file=non_cov_reg_2_file,
            lo_pos_assem1_file=lo_pos_assem1_file,
            lo_pos_assem2_file=lo_pos_assem2_file,
            lo_pos_assem1_0_file=lo_pos_assem1_0_file,
            lo_pos_assem2_0_file=lo_pos_assem2_0_file,
            nb_x_chr=nbXchr,
            arguments=arguments
        }
    }

    call ttmars.combine_ttmars_t {input: tsvs=ttmars_t.tsvOutput}

    output {
        File tsvOutput = combine_ttmars_t.tsvOutput
        File? liftoverNonCovReg1File = liftover_t.non_cov_reg_1_file
        File? liftoverNonCovReg2File = liftover_t.non_cov_reg_2_file
        File? liftoverLoPosAssem1File = liftover_t.lo_pos_assem1_file
        File? liftoverLoPosAssem2File = liftover_t.lo_pos_assem2_file
        File? liftoverLoPosAssem10File = liftover_t.lo_pos_assem1_0_file
        File? liftoverLoPosAssem20File = liftover_t.lo_pos_assem2_0_file
    }
}