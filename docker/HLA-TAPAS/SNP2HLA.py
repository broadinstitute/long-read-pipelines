# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
from platform import platform



########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))



def SNP2HLA(_input, _reference_panel, _out,
            _mem = "2000m", _marker_window_size=1000, _tolerated_diff=.15,
            _beagle_NTHREADS=1, _beagle_ITER=5, _beagle_MAP=None,
            _dependency="./"):



    ### Optional Arguments check.

    if not os.path.exists(_dependency):
        print(std_ERROR_MAIN_PROCESS_NAME + "The path(folder) of depedency('{}') doesn't exist. Please check it again.".format(_dependency))
        sys.exit()

    p_Mb = re.compile(r'\d+m')
    p_Gb = re.compile(r'\d+[gG]')

    if p_Mb.match(_mem):
        pass # No problem.
    elif p_Gb.match(_mem):
        _mem = re.sub(r'[gG]', '000m', _mem)
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Given Java memory value('{}') has bizzare representation. Please check it again.".format(_mem))
        sys.exit()



    ### Check dependencies.

    _plink = os.path.join(_dependency, "plink_mac" if not bool(re.search(pattern="Linux", string=platform())) else "plink") #plink v1.9
    _beagle = os.path.join(_dependency, "beagle.jar")   # Beagle(v4.1)
    #_linkage2beagle = os.path.join(_dependency, "linkage2beagle.jar")
    #_beagle2linkage = os.path.join(_dependency, "beagle2linkage.jar")
    _vcf2gprobs = os.path.join(_dependency, "vcf2gprobs.jar")
    _merge_table = os.path.join("src/merge_tables.pl")
    _parse_dosage = os.path.join("src/ParseDosage.csh")


    if not os.path.exists(_plink):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare PLINK(v1.90) in 'dependency/' folder.")
        sys.exit()
    if not os.path.exists(_beagle):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare Beagle(v4.1) in 'dependency/' folder.")
        sys.exit()
    #if not os.path.exists(_linkage2beagle):
    #    print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'linkage2beagle.jar' in 'dependency/' folder.")
    #    sys.exit()
    #if not os.path.exists(_beagle2linkage):
    #    print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'beagle2linkage.jar' in 'dependency/' folder.")
    #    sys.exit()
    if not os.path.exists(_merge_table):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'merge_tables.pl' in 'src/' folder.")
        sys.exit()
    if not os.path.exists(_parse_dosage):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'ParseDosage.csh' in 'src/' folder.")
        sys.exit()



    ### Intermediate path.

    OUTPUT = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(OUTPUT)):
        INTERMEDIATE_PATH = os.path.dirname(OUTPUT)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)
    else:
        # If `os.path.dirname(OUTPUT)` doesn't exist, then it means the output of MakeReference should be genrated in current directory.
        INTERMEDIATE_PATH = "./"


    JAVATMP = _out+".javatmpdir"
    os.system("mkdir -p " + JAVATMP)



    ### Setting commands

    PLINK = ' '.join([_plink, "--silent", "--allow-no-sex"]) # "--noweb" won't be included because it is Plink1.9
    BEAGLE = ' '.join(["java", "-Djava.io.tmpdir="+JAVATMP, "-Xmx"+_mem, "-jar", _beagle])
    #LINKAGE2BEAGLE = ' '.join(["java", "-Djava.io.tmpdir="+JAVATMP, "-Xss5M -Xmx"+_mem, "-jar", _linkage2beagle])
    #BEAGLE2LINKAGE = ' '.join(["java", "-Djava.io.tmpdir="+JAVATMP, "-Xmx"+_mem, "-jar", _beagle2linkage])

    # MERGE = ' '.join(["perl", _merge_table])
    MERGE = _merge_table
    PARSEDOSAGE = _parse_dosage



    ### Control Flags
    EXTRACT_MHC = 1
    FLIP = 1
    IMPUTE = 1


    print("SNP2HLA: Performing HLA imputation for dataset {}".format(_input))
    print("- Java memory = {}(Mb)".format(_mem))
    print("- Beagle(v4.1) window size = \"{}\" markers".format(_marker_window_size))


    index= 1
    __MHC__ = _out+".MHC"


    if EXTRACT_MHC:

        print("[{}] Extracting SNPs from the MHC.".format(index)); index += 1
        #MAF >1% as imputation threshold
        command = ' '.join([PLINK, "--bfile", _input, "--chr 6", "--from-mb 28 --to-mb 34", "--maf 0.01", "--make-bed", "--out", __MHC__])
        print(command)
        os.system(command)


    if FLIP:

        print("[{}] Performing SNP quality control.".format(index)); index += 1

        ### Identifying non-A/T non-C/G SNPs to flip
        command = ' '.join(["echo", "SNP 	POS	A1	A2", ">", OUTPUT+".tmp1"])
        print(command)
        os.system(command)
        command = ' '.join(["cut", "-f2,4-", __MHC__+".bim", ">>", OUTPUT+".tmp1"])
        print(command)
        os.system(command)

        command = ' '.join(["echo", "SNP    POSR    A1R A2R", ">", OUTPUT+".tmp2"])
        print(command)
        os.system(command)
        command = ' '.join(["cut", "-f2,4-", _reference_panel+".bim", ">>", OUTPUT+".tmp2"])
        print(command)
        os.system(command)

        command = ' '.join([MERGE, OUTPUT+".tmp2", OUTPUT+".tmp1", "SNP", "|", "grep -v -w NA", ">", OUTPUT+".SNPS.alleles"])
        print(command)
        os.system(command)



        ### < Major flip 1 > ###

        command = ' '.join(["awk", "'{if ($3 != $6 && $3 != $7){print $1}}'", OUTPUT+".SNPS.alleles", ">", OUTPUT+".SNPS.toflip1"])
        print(command)
        os.system(command)

        command = ' '.join([PLINK, "--bfile", __MHC__, "--flip", OUTPUT+".SNPS.toflip1", "--make-bed", "--out", __MHC__+".FLP"])
        print(command)
        os.system(command)

        ## Calculating allele freqeuncy
        command = ' '.join([PLINK, "--bfile", __MHC__+".FLP", "--freq", "--out", __MHC__+".FLP.FRQ"])
        print(command)
        os.system(command)


        command = ' '.join(["sed 's/A1/A1I/g'", __MHC__+".FLP.FRQ.frq", "|", "sed 's/A2/A2I/g'", "|", "sed 's/MAF/MAF_I/g'", ">", OUTPUT+".tmp"])
        print(command)
        os.system(command)



        command = ' '.join(["mv", OUTPUT+".tmp", __MHC__+".FLP.FRQ"])
        print(command)
        os.system(command)

        command = ' '.join([MERGE, _reference_panel+".FRQ.frq", __MHC__+".FLP.FRQ.frq", "SNP", "|", "grep -v -w NA", ">", OUTPUT+".SNPS.frq"])
        print(command)
        os.system(command)




        ### < Major flip 2 > ### (*.parsed file)
        command = ' '.join(["sed 's/ /\t/g'", OUTPUT+".SNPS.frq", "|",
                            'awk \'{if ($3 != $8){print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $9 "\t" $8 "\t" 1-$10 "\t*"}else{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 "\t."}}\'',
                            ">", OUTPUT+".SNPS.frq.parsed"])
        print(command)
        os.system(command)



        ### < Major flip 3 > ###
        # Finding A/T and C/G SNPs
        command = ' '.join(['awk \'{if (($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C")){if ($4 > $7){diff=$4 - $7; if ($4 > 1-$7){corrected=$4-(1-$7)}else{corrected=(1-$7)-$4}}else{diff=$7-$4;if($7 > (1-$4)){corrected=$7-(1-$4)}else{corrected=(1-$4)-$7}};print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" diff "\t" corrected}}\'',
                            OUTPUT+".SNPS.frq.parsed", ">", OUTPUT+".SNPS.ATCG.frq"])
        print(command)
        os.system(command)



        ### < Major flip 4 > ###

        # Identifying A/T and C/G SNPs to flip or remove
        command = ' '.join(["awk '{if ($10 < $9 && $10 < ", str(_tolerated_diff) ,"){print $1}}'", OUTPUT+".SNPS.ATCG.frq", ">", OUTPUT+".SNPS.toflip2"]) # modified by Saori
        print(command)
        os.system(command)

        command = ' '.join(["awk '{if ($4 > 0.5){print $1}}'", OUTPUT+".SNPS.ATCG.frq", ">", OUTPUT+".SNPS.toremove"]) # increased from 0.4 by Saori
        print(command)
        os.system(command)


        ## Identifying non A/T and non C/G SNPs to remove
        command = ' '.join(['awk \'{if (!(($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C"))){if ($4 > $7){diff=$4 - $7;}else{diff=$7-$4}; if (diff > \'%f\'){print $1}}}\''%(_tolerated_diff),
                            OUTPUT+".SNPS.frq.parsed", ">>", OUTPUT+".SNPS.toremove"])
        print(command)
        os.system(command)


        command = ' '.join(['awk \'{if (($2 != "A" && $2 != "C" && $2 != "G" && $2 != "T") || ($3 != "A" && $3 != "C" && $3 != "G" && $3 != "T")){print $1}}\'',
                            OUTPUT+".SNPS.frq.parsed", ">>", OUTPUT+".SNPS.toremove"])
        print(command)
        os.system(command)

        command = ' '.join(['awk \'{if (($2 == $5 && $3 != $6) || ($3 == $6 && $2 != $5)){print $1}}\'',
                            OUTPUT+".SNPS.frq.parsed", ">>", OUTPUT+".SNPS.toremove"])
        print(command)
        os.system(command)
        
        command = ' '.join(["sort", OUTPUT+".SNPS.toremove", "|","uniq > temp" ])
        print(command)
        os.system(command)
 
        command = ' '.join(['mv temp',OUTPUT+".SNPS.toremove" ] )
        print(command)
        os.system(command)

        ## Making QCd SNP file
        command = ' '.join([PLINK, "--bfile", __MHC__+".FLP", "--geno 0.2", "--exclude", OUTPUT+".SNPS.toremove", "--flip", OUTPUT+".SNPS.toflip2", "--make-bed", "--out", __MHC__+".QC"])
        print(command)
        os.system(command)

        command = ' '.join([PLINK, "--bfile", __MHC__+".QC", "--freq", "--out", __MHC__+".QC.FRQ"])
        print(command)
        os.system(command)

        command = ' '.join(["sed 's/A1/A1I/g'", __MHC__+".QC.FRQ.frq", "|", "sed 's/A2/A2I/g'", "|", "sed 's/MAF/MAF_I/g'", ">", OUTPUT+".tmp"])
        print(command)
        os.system(command)

        command = ' '.join(["mv", OUTPUT+".tmp", __MHC__+".QC.FRQ.frq"])
        print(command)
        os.system(command)

        command = ' '.join([MERGE, _reference_panel+".FRQ.frq", __MHC__+".QC.FRQ.frq", "SNP", "|", "grep -v -w NA", ">", OUTPUT+".SNPS.QC.frq"])
        print(command)
        os.system(command)


        command = ' '.join(["cut -f2", OUTPUT+".SNPS.QC.frq", "|", "awk '{if (NR > 1){print $1}}'", ">", OUTPUT+".SNPS.toinclude"])
        print(command)
        os.system(command)

        command = ' '.join(['echo "SNP     POS    A1    A2"', ">", OUTPUT+".tmp1"])
        print(command)
        os.system(command)

        command = ' '.join(["cut -f2,4-", __MHC__+".QC.bim", ">>", OUTPUT+".tmp1"])
        print(command)
        os.system(command)

        command = ' '.join([MERGE, OUTPUT+".tmp2", OUTPUT+".tmp1", "SNP", "|", 'awk \'{if (NR > 1){if ($5 != "NA"){pos=$5}else{pos=$2}; print "6\t" $1 "\t0\t" pos "\t" $3 "\t" $4}}\'',
                            ">", __MHC__+".QC.bim"])
        print(command)
        os.system(command)



        ### < Making *.vcf file for imputation (Beagle v4.x.x.) > ###

        """
        # Extracting SNPs and recoding QC'd file as vcf
        plink --bfile $MHC.QC --extract $OUTPUT.SNPS.toinclude --make-bed --out $MHC.QC.reorder
        plink --bfile $MHC.QC.reorder --recode vcf-iid --a1-allele $REFERENCE.markers 4 1 --out $MHC.QC
        
        """

        # Extracting SNPs and recoding QC'd file as vcf
        command = ' '.join([PLINK, "--bfile", __MHC__+".QC", "--extract", OUTPUT+".SNPS.toinclude", "--make-bed --out", __MHC__+".QC.reorder" ])
        print(command)
        os.system(command)

        command = ' '.join([PLINK, "--bfile", __MHC__+".QC.reorder", "--recode vcf-iid --a1-allele", _reference_panel+".markers 4 1", "--out", __MHC__+".QC" ])
        print(command)
        os.system(command)


        # Just in case of storage problem.
        command = ' '.join(["gzip -f", __MHC__+".QC.vcf"])
        print(command)
        os.system(command)



        ## Remove temporary files.
        os.system(' '.join(["rm ", OUTPUT+".tmp{1,2}"]))
        os.system(' '.join(["rm ", __MHC__+".FLP.*"]))
        os.system(' '.join(["rm ", __MHC__+".QC.{ped,map}"]))
        os.system(' '.join(["rm ", OUTPUT+".SNPS.*"]))
	
        # Beagle v4.x.x.
        os.system(' '.join(["rm ", __MHC__+".{bed,bim,fam,log}"]))
        os.system(' '.join(["rm ", __MHC__+".QC.reorder.*"]))
        os.system(' '.join(["rm ", __MHC__+".QC.{bed,bim,fam,log}"]))
        os.system(' '.join(["rm ", __MHC__+".QC.FRQ.{frq,log}"]))




    if IMPUTE:

        """
        if ($#argv >= 8) then
            beagle ref=$REFERENCE.bgl.vcf.gz gt=$MHC.QC.vcf impute=true gprobs=true nthreads=$THREAD chrom=6 niterations=$ITER lowmem=true out=$OUTPUT.bgl map=$MAP
	    else
            beagle ref=$REFERENCE.bgl.vcf.gz gt=$MHC.QC.vcf impute=true gprobs=true nthreads=$THREAD chrom=6 niterations=$ITER lowmem=true out=$OUTPUT.bgl
        """

        print("[{}] Performing HLA imputation.".format(index)); index += 1

        ## new beagle (>v4), assuming 4 threads and 10 interations
        command=' '.join([BEAGLE,
                          "gt=" + __MHC__+".QC.vcf.gz",
                          'ref='+_reference_panel+".bgl.phased.vcf.gz",
                          'impute=true',
                          'gprobs=true',
                          'nthreads={}'.format(_beagle_NTHREADS),
                          'chrom=6',
                          'niterations={}'.format(_beagle_ITER),
                          'lowmem=true',
                          ('map={}'.format(_beagle_MAP) if bool(_beagle_MAP) else ''),
                          'out='+ OUTPUT+".bgl.phased"])

        print(command)
        os.system(command)

        __IMPUTED__ = OUTPUT+".bgl.phased.vcf.gz"



        """
        (1) Imputation result in *.vcf.gz file
        (2) Imputation result in *.{bed,bim,fam} files (*.vcf.gz => *.{bed,bim,fam})
        (2) Dosage file (*.gprobs => *.dosage)
        """


        # (2) Imputation result in *.{bed,bim,fam} files (*.vcf.gz => *.{bed,bim,fam})
        command = ' '.join([PLINK, "--make-bed", "--vcf", __IMPUTED__, "--a1-allele {} 4 1".format(_reference_panel+".markers"), "--out", OUTPUT])
        #print(command)
        #os.system(command)


        # (3) Dosage file
        command = ' '.join(["gunzip -c", __IMPUTED__, "|", "cat", "|", "java -jar {} > {}".format(_vcf2gprobs, OUTPUT+".bgl.gprobs")])
        #print(command)
        #os.system(command)

        __gprobs__ = OUTPUT+".bgl.gprobs"


        command = ' '.join(["tail -n +2 {}".format(__gprobs__), "|",
                            PARSEDOSAGE, "- > {}".format(OUTPUT+".dosage")])
        #print(command)
        #os.system(command)



        os.system(' '.join(["rm ", __MHC__+".QC.vcf.gz"]))
        os.system(' '.join(["rm -rf", JAVATMP]))

        


    print("Done\n")

    return 0



if __name__ == "__main__" :

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        < SNP2HLA.py >
    
        SNP2HLA: Imputation of HLA amino acids and classical alleles from SNP genotypes
    
        Author: Sherman Jia (xiaomingjia@gmail.com)
                + Small modifications by Buhm Han (buhmhan@broadinstitute.org): 8/7/12
                + Extensive modifications by Phil Stuart (pstuart@umich.edu) on 4/19/16 to allow use of Beagle 4.1:
                    verified to work with 22Feb16, 15Apr16, and 03May16 versions of Beagle.
                + Small modifications by Yang Luo (yangluo@broadinstitute.org): 09/30/16: verfiied working with Bealge 4.1 27Jun16.
                + Recoded to Python and updated by Wanson Choi(wschoi.bhlab@gmail.com) : 2019/02/06 
                
                
        DESCRIPTION: This script runs imputation of HLA amino acids and classical alleles using SNP data.        
        
        INPUTS:
        1. Plink dataset (*.bed/bim/fam)
        2. Reference dataset (*.bgl.phased.vcf.gz(Beagle 4.1), *.markers(Beagle 3.0.4); *.fam/.bim/.FRQ.frq in PLINK format)
    
        DEPENDENCIES: (download and place in the same folder as this script)
        1. PLINK (1.9)  (Will not work with older Plink 1.07)
        2. Beagle (4.1) (Need to rename java executable as beagle.jar)
        3. vcf2gprobs.jar (Beagle utility for generating a Beagle v3 genotypes probability file from a Beagle 4.1 vcf file with GT field data)
        4. [Optional] If genetic_map_file argument is specified, PLINK format genetic map on cM scale 
            (plink.chr6.GRCh36.map, downloaded from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)

        USAGE:
        python3 SNP2HLA.py 
                --input `DATA (.bed/.bim/.fam)` 
                --reference `REFERENCE (.bgl.phased.vcf.gz/.markers/.fam/.bim/.bed/.FRQ.frq)` 
                --out `OUTPUT`
    
        (ex1)
        python3 SNP2HLA.py 
                --input data/1958BC 
                --reference data/Reference_Panel_bglv4/HM_CEU_REF.hg18.imgt3320.bglv4 
                --out TEST_ex1_SNP2HLA
    
        (ex2)
        python3 SNP2HLA.py 
                --input data/1958BC 
                --reference data/Reference_Panel_bglv4/HM_CEU_REF.hg18.imgt3320.bglv4 
                --out TEST_ex2_SNP2HLA 
                --dependency dependency/ 
                --java-mem 10G 
                --nthreads 4 
                --iter 10

    #################################################################################################
                                     '''),
                                     add_help=False)


    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--input", "-i", help="\nInput Plink data file prefix(.bed/.bim/.fam)\n\n", required=True)
    parser.add_argument("--out", "-o", help="\nOutput file prefix\n\n", required=True)
    parser.add_argument("--reference", "-rf", help="\nThe file prefix of reference panel for imputation.\n\n", required=True)

    parser.add_argument("--tolerated-diff", help="\nTolerated diff (default : 0.15).\n\n", default=0.15)
    parser.add_argument("--dependency", help="\nPath(folder) to dependecy software.\n\n", default="./") # Default : the folder where SNP2HLA.py is implemented.

    # Beagle(v4).
    parser.add_argument("--java-mem", "-mem", help="\nJava memory allocation(ex. '2000m', '2g', or '2G').\n\n", default='2000m')
    parser.add_argument("--marker-window", help="\n(Beagle4.1) Marker window size for imputation (default: 1000).\n\n", default=1000)
    parser.add_argument("--nthreads", help="\n(Beagle4.1) The number of threads to be used in imputation. (default: 1)\n\n", default=1)
    parser.add_argument("--iter", help="\n(Beagle4.1) The number of iteration in imputation (default: 5).\n\n", default=5)
    parser.add_argument("--plink-genetic-map", help="\n(Beagle4.1) Plink genetic map file to be utilized in imputation (default: None).\n\n", default=None)


    ##### <for Test> #####


    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)



    SNP2HLA(args.input, args.reference, args.out,
            _mem=args.java_mem, _marker_window_size=args.marker_window, _tolerated_diff=float(args.tolerated_diff),
            _beagle_NTHREADS=args.nthreads, _beagle_ITER=args.iter, _beagle_MAP=args.plink_genetic_map,
            _dependency=args.dependency)
