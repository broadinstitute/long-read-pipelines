version 1.0

workflow Benchmark {

    meta {
        description: "A workflow to calculate sensitivity and precision of a germline variant calling pipeline by comparing a 'call' vcf produced by the pipeline to a gold standard 'truth' vcf.  Allows for stratification based on interval lists, bed files, or variant types defined according to GATK SelectVariants.  Borrowed and adapted from the Broad Institute's Hydrogen/Palantir repo, courtesy of Michael Gatzen (https://github.com/broadinstitute/palantir-workflows/tree/mg_benchmark_compare/BenchmarkVCFs ; permalink: https://github.com/broadinstitute/palantir-workflows/blob/0bf48efc6de818364993e46d89591a035cfd80c7/BenchmarkVCFs/BenchmarkVCFs.wdl)."
    }

    parameter_meta {
        evalVcf: {description: "vcfs to be evaluated"}
        evalLabel: {description: "label to identify vcf to be evaluated"}
        evalVcfIndex: {description: "vcf index for evalVcf"}
        evalBam: {description: "bam file contaning the reads that generated the evalVcf"}
        evalBamLabel: {description: "label to use for the evalBam in IGV"}
        truthVcf: {description: "truth vcf against which to evaluate"}
        truthLabel: {description: "label by which to indentify truth set"}
        truthBam: {description: "bam file contaning the reads that generated the truthVcf"}
        truthBamLabel: {description: "label to use for the truthBam in IGV"}
        confidenceInterval: {description: "confidence interval for truth set (can be bed or picard interval_list)"}
        ref_map_file: {description: "table indicating reference sequence and auxillary file locations" }
        hapMap: {description: "reference haplotype map for CrosscheckFingerprints"}
        stratIntervals: {description: "intervals for stratifiction (can be picard interval_list or bed format)"}
        stratLabels: {description: "labels by which to identify stratification intervals (must be same length as stratIntervals)"}
        jexlVariantSelectors: {description: "variant types to select over (defined by jexl fed to GATK SelectVariants)"}
        variantSelectorLabels: {description: "labels by which to identify variant selectors (must be same length as jexlVariantSelectors)"}
        doIndelLengthStratification: {description: "whether or not to perform stratification by indel length"}
        requireMatchingGenotypes: {description: "whether to require genotypes to match in order to be a true positive"}
        truthIsSitesOnlyVcf: {description: "whether the truth VCF is a sites-only VCF file without any sample information"}
        gatkTag: {description: "version of gatk docker to use.  Defaults to 4.0.11.0"}
        analysisRegion: {description: "if provided (gatk format, single interval e.g., 'chr20', or 'chr20:1-10') all the analysis will be performed only within the region."}
        passingOnly: {description:"Have vcfEval only consider the passing variants"}
        vcfScoreField: {description:"Have vcfEval use this field for making the roc-plot. If this is an info field (like VSQLOD) it should be provided as INFO.VQSLOD, otherewise it is assumed to be a format field."}
        gatkJarForAnnotation: {description:"GATK jar that can calculate necessary annotations for jexl Selections when using VCFEval."}
        annotationNames: {description:"Annotation arguments to GATK (-A argument, multiple OK)"}
        dummyInputForTerraCallCaching: {description:"When running on Terra, use workspace.name as this input to ensure that all tasks will only cache hit to runs in your own workspace. This will prevent call caching from failing with 'Cache Miss (10 failed copy attempts)'. Outside of Terra this can be left empty. This dummy input is only needed for tasks that have no inputs specific to the sample being run (such as CreateIntervalList which does not take in any sample data)."}
    }

    input{

        File evalVcf
        String evalLabel
        File evalVcfIndex
        File? evalBam
        String? evalBamLabel

        File truthVcf
        String truthLabel
        File truthVcfIndex
        File? truthBam
        String? truthBamLabel

        File confidenceInterval

        File ref_map_file

        String? analysisRegion
        File? hapMap

        Array[File] stratIntervals = []
        Array[String] stratLabels = []
        Array[String]? jexlVariantSelectors
        Array[String]? variantSelectorLabels

        Int? threadsVcfEval = 2
        Boolean doIndelLengthStratification = true
        Int? preemptible
        String gatkTag="4.0.11.0"
        Boolean requireMatchingGenotypes = true
        Boolean truthIsSitesOnlyVcf = false
        File? gatkJarForAnnotation
        Array[String]? annotationNames
        Boolean enableRefOverlap = false
        Boolean passingOnly = true
        String? vcfScoreField
        String? dummyInputForTerraCallCaching
    }

    # Get ref info:
    Map[String, String] ref_map = read_map(ref_map_file)

    if (defined(analysisRegion)) {
        call CreateIntervalList {
            input:
                reference = ref_map["fasta"],
                reference_index = ref_map["fai"],
                reference_dict = ref_map["dict"],
                interval_string = select_first([analysisRegion]),
                gatkTag = gatkTag,
                dummyInputForTerraCallCaching = dummyInputForTerraCallCaching
        }
    }

    Array[File] actualStratIntervals = flatten([[""], stratIntervals])
    Array[String] actualStratLabels = flatten([[""], stratLabels])
    Array[String] actualSelectorLabels = select_first([variantSelectorLabels,[""]])
    Array[String] actualSelectorJEXL = select_first([jexlVariantSelectors,[""]])

    #check that lengths of different arrays are compatible
    if (length(actualStratLabels)!= length(actualStratIntervals)) {
         call ErrorWithMessage as Error6 {
             input:
                 message="Stratification vcf list is length "+length(actualStratIntervals)+" while stratification labels list is length "+length(actualStratLabels)
         }
    }

    if (length(actualSelectorLabels) != length(actualSelectorJEXL)) {
        call ErrorWithMessage as Error7 {
            input:
                message="Variant selector list is length "+length(actualSelectorJEXL)+" while labels list is "+length(actualSelectorLabels)
        }
    }

    if (defined(hapMap)) {
        call MatchEvalTruth as Match {
            input:
                evalVcf = evalVcf,
                truthVcf = truthVcf,
                evalVcfIndex = evalVcfIndex,
                truthVcfIndex = truthVcfIndex,
                hapMap = select_first([hapMap]),
                gatkTag = gatkTag,
                preemptible = preemptible
        }
    }
    Array[String] indelLabels=["deletion","insertion","indel_fine_m20","indel_fine_m19","indel_fine_m18","indel_fine_m17","indel_fine_m16","indel_fine_m15",
                               "indel_fine_m14","indel_fine_m13","indel_fine_m12","indel_fine_m11","indel_fine_m10","indel_fine_m9","indel_fine_m8","indel_fine_m7",
                               "indel_fine_m6","indel_fine_m5","indel_fine_m4","indel_fine_m3","indel_fine_m2","indel_fine_m1","indel_fine_1","indel_fine_2","indel_fine_3",
                               "indel_fine_4","indel_fine_5","indel_fine_6","indel_fine_7","indel_fine_8","indel_fine_9","indel_fine_10","indel_fine_11","indel_fine_12",
                               "indel_fine_13","indel_fine_14","indel_fine_15","indel_fine_16","indel_fine_17","indel_fine_18","indel_fine_19","indel_fine_20","indel_coarse_m30.0",
                               "indel_coarse_m25.0","indel_coarse_m20.0","indel_coarse_m15.0","indel_coarse_m10.0","indel_coarse_m5.0","indel_coarse_0.0","indel_coarse_5.0",
                               "indel_coarse_10.0","indel_coarse_15.0","indel_coarse_20.0","indel_coarse_25.0","indel_coarse_30.0"]

    Array[String] indelJexl=["vc.isSimpleIndel()  && vc.getIndelLengths().0<0","vc.isSimpleIndel()  && vc.getIndelLengths().0>0","vc.isSimpleIndel()  && vc.getIndelLengths().0==-20",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==-19","vc.isSimpleIndel()  && vc.getIndelLengths().0==-18","vc.isSimpleIndel()  && vc.getIndelLengths().0==-17",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==-16","vc.isSimpleIndel()  && vc.getIndelLengths().0==-15","vc.isSimpleIndel()  && vc.getIndelLengths().0==-14",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==-13","vc.isSimpleIndel()  && vc.getIndelLengths().0==-12","vc.isSimpleIndel()  && vc.getIndelLengths().0==-11",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==-10","vc.isSimpleIndel()  && vc.getIndelLengths().0==-9","vc.isSimpleIndel()  && vc.getIndelLengths().0==-8",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==-7","vc.isSimpleIndel()  && vc.getIndelLengths().0==-6","vc.isSimpleIndel()  && vc.getIndelLengths().0==-5",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==-4","vc.isSimpleIndel()  && vc.getIndelLengths().0==-3","vc.isSimpleIndel()  && vc.getIndelLengths().0==-2",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==-1","vc.isSimpleIndel()  && vc.getIndelLengths().0==1","vc.isSimpleIndel()  && vc.getIndelLengths().0==2",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==3","vc.isSimpleIndel()  && vc.getIndelLengths().0==4","vc.isSimpleIndel()  && vc.getIndelLengths().0==5",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==6","vc.isSimpleIndel()  && vc.getIndelLengths().0==7","vc.isSimpleIndel()  && vc.getIndelLengths().0==8",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==9","vc.isSimpleIndel()  && vc.getIndelLengths().0==10","vc.isSimpleIndel()  && vc.getIndelLengths().0==11",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==12","vc.isSimpleIndel()  && vc.getIndelLengths().0==13","vc.isSimpleIndel()  && vc.getIndelLengths().0==14",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==15","vc.isSimpleIndel()  && vc.getIndelLengths().0==16","vc.isSimpleIndel()  && vc.getIndelLengths().0==17",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0==18","vc.isSimpleIndel()  && vc.getIndelLengths().0==19","vc.isSimpleIndel()  && vc.getIndelLengths().0==20",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0<-27.5 && vc.getIndelLengths().0>-32.5","vc.isSimpleIndel()  && vc.getIndelLengths().0<-22.5 && vc.getIndelLengths().0>-27.5",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0<-17.5 && vc.getIndelLengths().0>-22.5","vc.isSimpleIndel()  && vc.getIndelLengths().0<-12.5 && vc.getIndelLengths().0>-17.5",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0<-7.5 && vc.getIndelLengths().0>-12.5","vc.isSimpleIndel()  && vc.getIndelLengths().0<-2.5 && vc.getIndelLengths().0>-7.5",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0<2.5 && vc.getIndelLengths().0>-2.5","vc.isSimpleIndel()  && vc.getIndelLengths().0<7.5 && vc.getIndelLengths().0>2.5",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0<12.5 && vc.getIndelLengths().0>7.5","vc.isSimpleIndel()  && vc.getIndelLengths().0<17.5 && vc.getIndelLengths().0>12.5",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0<22.5 && vc.getIndelLengths().0>17.5","vc.isSimpleIndel()  && vc.getIndelLengths().0<27.5 && vc.getIndelLengths().0>22.5",
                             "vc.isSimpleIndel()  && vc.getIndelLengths().0<32.5 && vc.getIndelLengths().0>27.5"]

    scatter (indel in zip(indelLabels,indelJexl)) {
        VariantSelector indelSelectors = object{ jexl : indel.right,
                                         label : indel.left
                                         }
    }

    if (defined(jexlVariantSelectors)) {
        scatter (select in zip(actualSelectorLabels,actualSelectorJEXL)) {
            VariantSelector variantSelectors = object{ jexl: select.right,
                                                label : select.left
                                                }
        }
    }

    Array[VariantSelector] defaultVS = [object{ jexl: "vc.isIndel() && vc.getHetCount() == 1", label: "HetIndel" },
                                        object{jexl: "vc.isIndel() && vc.getHomVarCount() == 1", label: "HomVarIndel"},
                                        object{jexl: "vc.isSNP() && vc.getHetCount() == 1", label: "HetSNP"},
                                        object{jexl: "vc.isSNP() && vc.getHomVarCount() == 1", label: "HomVarSNP"}]
    Array[VariantSelector] actualVariantSelectors = flatten(select_all([defaultVS,variantSelectors]))

    if (defined(stratIntervals)) {
        scatter (actStratIL in actualStratIntervals) {

            # Required because of `miniwdl check`:
            String string_conversion_of_actStratIL = actStratIL

            if(string_conversion_of_actStratIL != "") {
                call ConvertIntervals as StratConvertIntervals {
                    input:
                        inputIntervals = actStratIL,
                        refDict = ref_map["dict"],
                        gatkTag = gatkTag,
                        subset_interval = CreateIntervalList.interval_list,
                        preemptible = preemptible,
                        dummyInputForTerraCallCaching = dummyInputForTerraCallCaching

                }
            }
        }
    }
    Array[File] stratBeds = select_all(flatten(select_all([[""],StratConvertIntervals.bed])))
    Array[File] stratILs = select_all(flatten(select_all([[""],StratConvertIntervals.intervalList])))

    scatter (strat in zip(zip(stratILs,stratBeds),actualStratLabels)) {
        Stratifier stratifiers = object {intervalList : strat.left.left,
                                    bed : strat.left.right,
                                    label : strat.right
                                    }
    }

    call ConvertIntervals as ConfidenceConvertIntervals {
        input:
            inputIntervals = confidenceInterval,
            refDict = ref_map["dict"],
            gatkTag = gatkTag,
            preemptible = preemptible,
            subset_interval = CreateIntervalList.interval_list,
            dummyInputForTerraCallCaching = dummyInputForTerraCallCaching
    }

    scatter (stratifier in stratifiers) {

        # Required because of `miniwdl check`
        String tmp_strat_interval_list = select_first([stratifier.intervalList])

        if (defined(stratifier.label) && defined(tmp_strat_interval_list)) {
            String stratLabel = select_first([stratifier.label,""])
            File stratIL = select_first([stratifier.intervalList,""])
            File stratBed = select_first([stratifier.bed,""])
            String outputPreStrat = evalLabel+"_"+truthLabel+"_"+stratLabel
        }
        String outputPrefix = select_first([outputPreStrat,evalLabel+"_"+truthLabel])

        call CheckForVariants as CheckForVariantsEval {
            input:
                vcf = evalVcf,
                vcfIndex = evalVcfIndex,
                confidenceIL = ConfidenceConvertIntervals.intervalList,
                stratIL = stratIL,
                gatkTag = gatkTag,
                preemptible = preemptible
        }

        call CheckForVariants as CheckForVariantsTruth {
            input:
                vcf = truthVcf,
                vcfIndex = truthVcfIndex,
                confidenceIL = ConfidenceConvertIntervals.intervalList,
                stratIL = stratIL,
                gatkTag = gatkTag,
                preemptible = preemptible
        }

        if (CheckForVariantsTruth.variantsFound && CheckForVariantsEval.variantsFound) {
            call VcfEval as StandardVcfEval {
                input:
                    truthVCF = truthVcf,
                    truthVCFIndex = truthVcfIndex,
                    evalVCF = evalVcf,
                    evalVCFIndex = evalVcfIndex,
                    confidenceBed = ConfidenceConvertIntervals.bed,
                    stratBed = stratBed,
                    ref = ref_map["fasta"],
                    refDict = ref_map["dict"],
                    refIndex = ref_map["fai"],
                    outputPre = outputPrefix+"_vcfeval",
                    threads = threadsVcfEval,
                    preemptible = preemptible,
                    requireMatchingGenotypes = requireMatchingGenotypes,
                    truthIsSitesOnlyVcf = truthIsSitesOnlyVcf,
                    passingOnly = passingOnly,
                    vcfScoreField = vcfScoreField,
                    enableRefOverlap = enableRefOverlap
            }

            call WriteXMLfile as VcfEvalWriteXMLfile {
                input:
                    input_files = flatten([select_all([StandardVcfEval.outVcf,ConfidenceConvertIntervals.bed,stratifier.bed]), select_all([evalBam, truthBam])]),
                    input_names = flatten([select_all([outputPrefix+"_vcfeval","confidence_intervals",stratifier.label]), select_all([evalBamLabel, truthBamLabel])]),
                    reference_version = ref_map["fasta"],
                    file_name = outputPrefix+"_vcfeval"
            }

            call CountUNKVcfEval {
                input:
                    vcf = StandardVcfEval.outVcf,
                    vcfIndex = StandardVcfEval.outVcfIndex,
                    gatkTag = gatkTag,
                    preemptible = preemptible
            }
        }

        String areVariants = if(CheckForVariantsTruth.variantsFound && CheckForVariantsEval.variantsFound) then "yes" else "no"
        call SummariseVcfEval {
            input:
                evalLabel = evalLabel,
                truthLabel = truthLabel,
                stratLabel = stratLabel,
                summaryFile = StandardVcfEval.outSummary,
                igvSession = VcfEvalWriteXMLfile.igv_session,
                areVariants = areVariants,
                unkSNP = CountUNKVcfEval.UNK_SNP,
                unkINDEL = CountUNKVcfEval.UNK_INDEL,
                preemptible = preemptible
        }
    }

    scatter ( i in range(length(stratifiers)) ) {
        AnnotatedVcfs annotatedVcfsList = object{vcfVcfEval : StandardVcfEval.outVcf[i],
                                            vcfVcfEvalIndex : StandardVcfEval.outVcfIndex[i],
                                            stratLabel : stratifiers[i].label,
                                            evalLabel : evalLabel,
                                            truthLabel : truthLabel,
                                            stratBed : stratBed[i],
                                            confidenceBed : ConfidenceConvertIntervals.bed,
                                            namePrefix : outputPrefix[i]
                                            }
    }


    scatter (indelCombo in cross(annotatedVcfsList,indelSelectors)) {
        EvalStratSelectorCombo evalStratIndelCombos = object{annotatedVcfs : indelCombo.left,
                                                    variantSelector : indelCombo.right
                                                    }
        }

    scatter (evalStratIndelCombo in evalStratIndelCombos) {
        String jexl = evalStratIndelCombo.variantSelector.jexl
        File? vcfVcfEval = evalStratIndelCombo.annotatedVcfs.vcfVcfEval
        File? vcfVcfEvalIndex = evalStratIndelCombo.annotatedVcfs.vcfVcfEvalIndex
        String evalIndelLabel = evalStratIndelCombo.annotatedVcfs.evalLabel
        String truthIndelLabel = evalStratIndelCombo.annotatedVcfs.truthLabel
        String? stratIndelLabel = evalStratIndelCombo.annotatedVcfs.stratLabel
        String indelLabel = evalStratIndelCombo.variantSelector.label
        File? stratIndelBed = evalStratIndelCombo.annotatedVcfs.stratBed
        File? confidenceBed = evalStratIndelCombo.annotatedVcfs.confidenceBed
        String namePrefix = evalStratIndelCombo.annotatedVcfs.namePrefix+"_"+indelLabel

        if (defined(vcfVcfEval) && defined(vcfVcfEvalIndex) && doIndelLengthStratification) {
            call EvalForVariantSelection as EvalIndelLengthVcfEval {
                input:
                    vcf = vcfVcfEval,
                    vcfIndex = vcfVcfEvalIndex,
                    jexl = jexl,
                    engine="VcfEval",
                    selectTPCall="CALL == 'TP'",
                    selectTPBase="BASE == 'TP'",
                    selectFN="(BASE == 'FN' || BASE == 'FN_CA')",
                    selectFP="(CALL == 'FP' || CALL == 'FP_CA')",
                    sampleCall="CALLS",
                    sampleBase="BASELINE",
                    gatkTag = gatkTag,
                    preemptible = preemptible,
                    gatkJarForAnnotation = gatkJarForAnnotation,
                    annotationNames = annotationNames,
                    reference = ref_map["fasta"],
                    refDict = ref_map["dict"],
                    refIndex = ref_map["fai"]
            }

            call WriteXMLfile as VcfEvalIndelWriteXMLfile {
                        input:
                            input_files = flatten([select_all([EvalIndelLengthVcfEval.selectedTPCall,EvalIndelLengthVcfEval.selectedTPBase,EvalIndelLengthVcfEval.selectedFP,EvalIndelLengthVcfEval.selectedFN,vcfVcfEval,confidenceBed,stratIndelBed]), select_all([evalBam, truthBam])]),
                            input_names = flatten([select_all(["TP_Eval","TP_Base","FP","FN","All_Variants","confidence_intervals",stratIndelLabel]), select_all([evalBamLabel, truthBamLabel])]),
                            reference_version = ref_map["fasta"],
                            file_name = namePrefix+"_vcfeval"
            }

            call SummariseForIndelSelection as VcfEvalSummariseForIndelSelection {
                        input:
                            evalLabel = evalIndelLabel,
                            truthLabel = truthIndelLabel,
                            stratLabel = stratIndelLabel,
                            indelLabel = indelLabel,
                            engine="VcfEval",
                            igvSession = VcfEvalIndelWriteXMLfile.igv_session,
                            TP_CALL = EvalIndelLengthVcfEval.TP_CALL,
                            TP_BASE = EvalIndelLengthVcfEval.TP_BASE,
                            FP = EvalIndelLengthVcfEval.FP,
                            FN = EvalIndelLengthVcfEval.FN,
                            preemptible = preemptible
            }
        }
    }




    scatter (selectorCombo in cross(annotatedVcfsList,actualVariantSelectors)) {
       EvalStratSelectorCombo evalStratSelectorCombos = object{annotatedVcfs : selectorCombo.left,
                                                    variantSelector : selectorCombo.right
                                                    }
    }
    scatter (evalStratSelectorCombo in evalStratSelectorCombos) {
                if (defined(evalStratSelectorCombo.annotatedVcfs.vcfVcfEval) && defined(evalStratSelectorCombo.annotatedVcfs.vcfVcfEvalIndex)) {
                    call EvalForVariantSelection as EvalSelectorVcfEval {
                        input:
                            vcf = evalStratSelectorCombo.annotatedVcfs.vcfVcfEval,
                            vcfIndex = evalStratSelectorCombo.annotatedVcfs.vcfVcfEvalIndex,
                            jexl = evalStratSelectorCombo.variantSelector.jexl,
                            engine="VcfEval",
                            selectTPCall="CALL == 'TP'",
                            selectTPBase="BASE == 'TP'",
                            selectFN="(BASE == 'FN' || BASE == 'FN_CA')",
                            selectFP="(CALL == 'FP' || CALL == 'FP_CA')",
                            sampleCall="CALLS",
                            sampleBase="BASELINE",
                            gatkTag = gatkTag,
                            preemptible = preemptible,
                            gatkJarForAnnotation = gatkJarForAnnotation,
                            annotationNames = annotationNames,
                            reference = ref_map["fasta"],
                            refDict = ref_map["dict"],
                            refIndex = ref_map["fai"]
                    }
                    call WriteXMLfile as VcfEvalSelectorWriteXMLfile {
                                input:
                                    input_files = flatten([select_all([EvalSelectorVcfEval.selectedTPCall,EvalSelectorVcfEval.selectedTPBase,EvalSelectorVcfEval.selectedFP,EvalSelectorVcfEval.selectedFN,
                                        evalStratSelectorCombo.annotatedVcfs.vcfVcfEval,evalStratSelectorCombo.annotatedVcfs.confidenceBed,evalStratSelectorCombo.annotatedVcfs.stratBed]), select_all([evalBam, truthBam])]),
                                    input_names = flatten([select_all(["TP_Eval","TP_Base","FP","FN","All_Variants","confidence_intervals",evalStratSelectorCombo.annotatedVcfs.stratLabel]), select_all([evalBamLabel, truthBamLabel])]),
                                    reference_version = ref_map["fasta"],
                                    file_name = evalStratSelectorCombo.annotatedVcfs.namePrefix+"_"+evalStratSelectorCombo.variantSelector.label+"_vcfeval"
                    }

                    call SummariseForVariantSelection as VcfEvalSummariseForVariantSelection {
                                input:
                                    evalLabel = evalStratSelectorCombo.annotatedVcfs.evalLabel,
                                    truthLabel = evalStratSelectorCombo.annotatedVcfs.truthLabel,
                                    stratLabel = evalStratSelectorCombo.annotatedVcfs.stratLabel,
                                    variantLabel = evalStratSelectorCombo.variantSelector.label,
                                    engine="VcfEval",
                                    igvSession = VcfEvalSelectorWriteXMLfile.igv_session,
                                    TP_CALL = EvalSelectorVcfEval.TP_CALL,
                                    TP_BASE = EvalSelectorVcfEval.TP_BASE,
                                    FP = EvalSelectorVcfEval.FP,
                                    FN = EvalSelectorVcfEval.FN,
                                    preemptible = preemptible
                    }
                }
    }

    Array[File] summaries = flatten([SummariseVcfEval.summaryOut,select_all(VcfEvalSummariseForVariantSelection.summaryOut),
                                    select_all(VcfEvalSummariseForIndelSelection.summaryOut)])

    call CombineSummaries {
        input:
            summaries = summaries,
            preemptible = preemptible
    }

    ################################################################################

    output {
        File summary = CombineSummaries.summaryOut
        Float snpPrecision = SummariseVcfEval.snpPrecision[0]
        Float indelPrecision = SummariseVcfEval.indelPrecision[0]
        Float snpRecall = SummariseVcfEval.snpRecall[0]
        Float indelRecall = SummariseVcfEval.indelRecall[0]
        Float snpF1Score = SummariseVcfEval.snpF1Score[0]
        Float indelF1Score = SummariseVcfEval.indelF1Score[0]
        Array[File?] snpRocs = StandardVcfEval.outSnpRoc
        Array[File?] nonSnpRocs = StandardVcfEval.outNonSnpRoc
    }
}

################################################################################
################################################################################
################################################################################

struct EvalTruthMatch {
        File truthVcf
        File truthVcfIndex
        File confidenceIntervals
        String truthLabel
        File evalVcf
        File evalVcfIndex
        String evalLabel
}

struct VariantSelector {
    String jexl
    String label
}

struct Stratifier {
    File? intervalList
    File? bed
    String? label
}

struct EvalStratCombo {
    EvalTruthMatch evalTruthMatch
    Stratifier stratifier
}

struct AnnotatedVcfs {
    File? vcfVcfEval
    File? vcfVcfEvalIndex
    String? stratLabel
    String evalLabel
    String truthLabel
    File? stratBed
    File? confidenceBed
    String namePrefix
}
struct EvalStratSelectorCombo {
    AnnotatedVcfs annotatedVcfs
    VariantSelector variantSelector
}

#Check to see if there are variants in the given vcf which overlap the confidence and stratification intervals
task CheckForVariants {
    input{
        File vcf
        File vcfIndex
        File confidenceIL
        File? stratIL
        Int? preemptible
        Int? memoryMaybe
        String gatkTag
        }
        Int memoryDefault = 16
        Int memoryJava = select_first([memoryMaybe,memoryDefault])
        Int memoryRam = memoryJava+2

        Int disk_size = 10 + ceil(size(vcf, "GB") + size(vcfIndex, "GB") + size(confidenceIL, "GB") + size(stratIL, "GB"))

    command <<<
        set -xeuo pipefail

    nVariants="$(gatk --java-options "-Xmx~{memoryJava}G" CountVariants -V ~{vcf} -L ~{confidenceIL} ~{"-L " + stratIL} -isr INTERSECTION | tail -1)"
    if [ "$nVariants" -gt "0" ]; then echo "true" > outBool.txt; else echo "false" > outBool.txt; fi
    >>>

    runtime {
            docker: "us.gcr.io/broad-gatk/gatk:"+gatkTag
            preemptible: select_first([preemptible,0])
            disks: "local-disk " + disk_size + " HDD"
            bootDiskSizeGb: "16"
            memory: memoryRam + " GB"
    }

    output {
        Boolean variantsFound = read_boolean("outBool.txt")
    }
}

#Evaluate evalVCF against truthVCF using vcfeval
task VcfEval {
    input{
        File truthVCF
        File truthVCFIndex

        File evalVCF
        File evalVCFIndex

        File confidenceBed

        File? stratBed

        File ref
        File refDict
        File refIndex

        String outputPre
        Boolean passingOnly

        String? vcfScoreField
        Int? preemptible
        String? memUser
        Int? threads

        Boolean requireMatchingGenotypes
        Boolean truthIsSitesOnlyVcf
        Boolean enableRefOverlap = false
    }
    String memDefault="16 GB"
    String mem = select_first([memUser,memDefault])

    Int cpu = select_first([threads,1])
    Int disk_size = 50 + ceil(size(truthVCF, "GB") + size(truthVCFIndex, "GB") + 2.2 * size(evalVCF, "GB") + size(evalVCFIndex, "GB") + size(confidenceBed, "GB") + size(stratBed, "GB") + size(ref, "GB") + size(refDict, "GB") + size(refIndex, "GB"))

    command <<<
    set -xeuo pipefail

    /bin/rtg-tools/rtg format -o rtg_ref ~{ref}
    /bin/rtg-tools/rtg vcfeval \
        ~{false="--all-records" true="" passingOnly} \
        ~{"--vcf-score-field=" + vcfScoreField} \
        ~{false="--squash-ploidy" true="" requireMatchingGenotypes} \
        ~{false="" true="--sample ALT" truthIsSitesOnlyVcf} \
        ~{true="--ref-overlap" false="" enableRefOverlap} \
        -b ~{truthVCF} -c ~{evalVCF} \
        -e ~{confidenceBed} ~{"--bed-regions " + stratBed} \
        ~{false="--output-mode combine" true="" truthIsSitesOnlyVcf} \
        --decompose -t rtg_ref \
        ~{"--threads "+threads} -o output_dir

    for f in output_dir/*; do
        mv $f ~{outputPre}_"$(basename "$f")";
    done

    /bin/rtg-tools/rtg rocplot --precision-sensitivity --title="~{outputPre} SNP"   --svg=~{outputPre}.snp.svg   ~{outputPre}_snp_roc.tsv.gz
    /bin/rtg-tools/rtg rocplot --precision-sensitivity --title="~{outputPre} INDEL" --svg=~{outputPre}.indel.svg ~{outputPre}_non_snp_roc.tsv.gz

    python3 -<<"EOF" ~{outputPre}_snp_roc.tsv.gz ~{outputPre}_non_snp_roc.tsv.gz ~{outputPre}_summary.csv
    import gzip
    import sys

    indel_sensitivity = 0
    indel_precision = 0
    indel_fscore = 0
    indel_TP_Base = 0
    indel_TP_Eval = 0
    indel_FP = 0
    indel_FN = 0

    snp_sensitivity = 0
    snp_precision = 0
    snp_fscore = 0
    snp_TP_Base = 0
    snp_TP_Eval = 0
    snp_FP = 0
    snp_FN = 0

    with gzip.open(sys.argv[1],"rt") as f_snp:
        for line in f_snp:
            try:
                snp_sensitivity = float(line.split()[6])
                snp_precision = float(line.split()[5])
                snp_fscore = float(line.split()[7])
                snp_TP_Eval = float(line.split()[3])
                snp_TP_Base = float(line.split()[1])
                snp_FP = float(line.split()[2])
                snp_FN = float(line.split()[4])
            except ValueError:
                continue
            except IndexError:
                continue
        f_snp.close()
    with gzip.open(sys.argv[2],"rt") as f_indel:
        for line in f_indel:
            try:
                indel_sensitivity = float(line.split()[6])
                indel_precision = float(line.split()[5])
                indel_fscore = float(line.split()[7])
                indel_TP_Eval = float(line.split()[3])
                indel_TP_Base = float(line.split()[1])
                indel_FP = float(line.split()[2])
                indel_FN = float(line.split()[4])
            except ValueError:
                continue
            except IndexError:
                continue
        f_indel.close()

    str_indel_sensitivity = str(indel_sensitivity)
    str_indel_precision = str(indel_precision)
    str_indel_fscore = str(indel_fscore)
    str_snp_sensitivity = str(snp_sensitivity)
    str_snp_precision = str(snp_precision)
    str_snp_fscore = str(snp_fscore)


    if indel_TP_Eval+indel_FP==0:
        str_indel_precision="NA"
    if indel_TP_Base+indel_FN==0:
        str_indel_sensitivity="NA"
    if str_indel_sensitivity=="NA" or str_indel_precision=="NA":
        str_indel_fscore="NA"

    if snp_TP_Eval+snp_FP==0:
        str_snp_precision="NA"
    if snp_TP_Base+snp_FN==0:
        str_snp_sensitivity="NA"
    if str_snp_sensitivity=="NA" or str_snp_precision=="NA":
        str_snp_fscore="NA"



    with open(sys.argv[3],"wt") as f_out:
        f_out.write(",".join(["Type","Precision","Recall","F1_Score","TP_Eval","TP_Base","FP","FN"])+"\n")
        f_out.write(",".join(["SNP",str_snp_precision,str_snp_sensitivity,str_snp_fscore,str(snp_TP_Eval),str(snp_TP_Base),str(snp_FP),str(snp_FN)])+"\n")
        f_out.write(",".join(["INDEL",str_indel_precision,str_indel_sensitivity,str_indel_fscore,str(indel_TP_Eval),str(indel_TP_Base),str(indel_FP),str(indel_FN)])+"\n")
        f_out.close()
    EOF
    >>>

    runtime {
        docker: "ckachulis/rtg-tools:0.1"
        preemptible: select_first([preemptible,0])
        memory: mem
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        Array[File] outs = glob("${outputPre}_*")
        File outSummary="${outputPre}_summary.csv"
        File outVcf="${outputPre}_output.vcf.gz"
        File outVcfIndex="${outputPre}_output.vcf.gz.tbi"
        File outSnpRocPlot="~{outputPre}.snp.svg"
        File outNonRocPlot="~{outputPre}.indel.svg"
        File outSnpRoc="${outputPre}_snp_roc.tsv.gz"
        File outNonSnpRoc="${outputPre}_non_snp_roc.tsv.gz"
    }
}

#Evaluate evalVCF against truthVCF using hap.py
task EvalHappy {
    input{
        File truthVCF
        File truthVCFIndex
        File evalVCF
        File confidenceBed
        File? stratBed
        File ref
        File refDict
        File refIndex
        String outputPre
        String? memUser
        Int? preemptible
        Int? threads
        String happyTag
    }
    String memDefault="16 GB"
    String mem = select_first([memUser,memDefault])

    Int cpu = select_first([threads,1])
    Int disk_size = 10 + ceil(size(truthVCF, "GB") + size(truthVCFIndex, "GB") + 2.2 * size(evalVCF, "GB") + size(confidenceBed, "GB") + size(stratBed, "GB") + size(ref, "GB") + size(refDict, "GB") + size(refIndex, "GB"))


    command <<<
        /opt/hap.py/bin/hap.py ~{truthVCF} ~{evalVCF} -f ~{confidenceBed} -r ~{ref} -V ~{"-T " + stratBed} -L --preprocess-truth ~{"--threads "+threads} -o ~{outputPre}
    >>>

    runtime {
        docker: "pkrusche/hap.py:"+happyTag
        memory: mem
        preemptible: select_first([preemptible,0])
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        Array[File] outs = glob("${outputPre}*")
        File outSummary="${outputPre}.summary.csv"
        File outVcf="${outputPre}.vcf.gz"
        File outVcfIndex="${outputPre}.vcf.gz.tbi"
    }
}

#Evaluate evalVCF against truthVCF using picard GenotypeConcordance
task EvalGATKGC {
    input{
        File truthVCF
        File truthVCFIndex
        File evalVCF
        File evalVCFIndex
        File intervalList
        File? stratIL
        File ref
        File refDict
        File refIndex
        String outputPre
        Int? preemptible
        Int? memoryMaybe
        String gatkTag
    }
    Int memoryDefault = 16
    Int memoryJava = select_first([memoryMaybe,memoryDefault])
    Int memoryRam = memoryJava+2
    Int disk_size = 10 + ceil(size(truthVCF, "GB") + size(truthVCFIndex, "GB") + 2.2 * size(evalVCF, "GB") + size(intervalList, "GB") + size(stratIL, "GB") + size(evalVCFIndex, "GB") + size(ref, "GB") + size(refDict, "GB") + size(refIndex, "GB"))


    command <<<
        gatk --java-options "-Xmx~{memoryJava}G" GenotypeConcordance -TV ~{truthVCF} -CV ~{evalVCF} -R ~{ref} -INTERVALS ~{intervalList} ~{"-INTERVALS " + stratIL} -USE_VCF_INDEX -OUTPUT_VCF -O ~{outputPre}
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:"+gatkTag
        preemptible: select_first([preemptible,0])
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: "16"
        memory: memoryRam + " GB"
    }

    output {
        Array[File] out = glob("${outputPre}*")
        File outSummary="${outputPre}.genotype_concordance_summary_metrics"
        File outCounts="${outputPre}.genotype_concordance_contingency_metrics"
        File outVcf="${outputPre}.genotype_concordance.vcf.gz"
        File outVcfIndex="${outputPre}.genotype_concordance.vcf.gz.tbi"
    }
}

#takes in either a .bed or .intervallist and returns both a .bed and .intervallist version of the input
task ConvertIntervals {
    input {
        File inputIntervals
        File? subset_interval
        Int? preemptible
        Int? memoryMaybe
        File refDict
        String gatkTag
        String? dummyInputForTerraCallCaching
    }

    Int memoryDefault = 16
    Int memoryJava = select_first([memoryMaybe,memoryDefault])
    Int memoryRam = memoryJava+2
    Int disk_size = 10 + ceil(3 * size(inputIntervals, "GB") + size(refDict, "GB"))

    command <<<
        set -xeuo pipefail

        # convert bed to interval_list, or copy interval_list
        if [[ ~{inputIntervals} == *.bed || ~{inputIntervals} == *.bed.gz ]]; then
                gatk --java-options "-Xmx~{memoryJava}G" \
                    BedToIntervalList \
                    -I ~{inputIntervals} \
                    -O initial_intervals.interval_list \
                    -SD ~{refDict}
        else
            cp ~{inputIntervals} initial_intervals.interval_list
        fi

        # optionally intersect interval_list with subset_interval

        if [ ! -z ~{subset_interval} ];  then
            gatk --java-options "-Xmx~{memoryJava}G" \
                IntervalListTools \
                -I initial_intervals.interval_list \
                -I ~{subset_interval} \
                -ACTION INTERSECT \
                -O intervals.interval_list
        else
            mv initial_intervals.interval_list intervals.interval_list
        fi

        # convert result to BED
        gatk --java-options "-Xmx~{memoryJava}G" \
            IntervalListToBed \
            -I intervals.interval_list \
            -O intervals.bed
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:"+gatkTag
        preemptible: select_first([preemptible,0])
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: "16"
        memory: memoryRam + " GB"
    }

    output {
        File bed="intervals.bed"
        File intervalList="intervals.interval_list"
    }
}

#For now, due to a bug, the ouput annotated VCF from hap.py has a hardcoded HG19 header, so the VCF header must be fixed to represent the correct reference.
#Additionally, due to a separate bug, UpdateVCFSequenceDictionary crashes in this case when trying to index a bgzipped vcf.  So for now have to perform indexing in a separate command.
task FixVcfHeader {
    input {
        File vcf
        File vcfIndex
        File ref
        File refDict
        File refIndex
        Int? preemptible
        Int? memoryMaybe
        String gatkTag
    }
    Int memoryDefault = 16
    Int memoryJava = select_first([memoryMaybe,memoryDefault])
    Int memoryRam = memoryJava+2
    Int disk_size = 10 + ceil(2.2 * size(vcf, "GB") + 2.2 * size(vcfIndex, "GB") + size(ref, "GB") + size(refDict, "GB") + size(refIndex, "GB"))

    command <<<
        set -xeuo pipefail

        gatk --java-options "-Xmx~{memoryJava}G" UpdateVCFSequenceDictionary -V ~{vcf} -O fixed.vcf.gz --source-dictionary ~{refDict} --replace --create-output-variant-index false
        gatk --java-options "-Xmx~{memoryJava}G" IndexFeatureFile -F fixed.vcf.gz
    >>>

    runtime {
            docker: "us.gcr.io/broad-gatk/gatk:"+gatkTag
            preemptible: select_first([preemptible,0])
            disks: "local-disk " + disk_size + " HDD"
            bootDiskSizeGb: "16"
            memory: memoryRam + " GB"
    }

    output {
        File outVcf="fixed.vcf.gz"
        File outVcfIndex="fixed.vcf.gz.tbi"
    }

}

#Count number of variants which were outside confidence region based on vcfeval annotated vcf
task CountUNKVcfEval {
    input {
        File? vcf=""
        File? vcfIndex=""
        Int? preemptible
        Int? memoryMaybe
        String gatkTag
    }
    Int memoryDefault = 16
    Int memoryJava = select_first([memoryMaybe,memoryDefault])
    Int memoryRam = memoryJava+2
    Int disk_size = 10 + ceil(size(vcf, "GB") + size(vcfIndex, "GB"))

    command <<<
        set -xeuo pipefail

        gatk --java-options "-Xmx~{memoryJava}G" SelectVariants -V ~{vcf} -O selected.unk.snp.vcf -select "(CALL == 'OUT')" --select-type-to-include SNP
        gatk --java-options "-Xmx~{memoryJava}G" SelectVariants -V ~{vcf} -O selected.unk.indel.vcf -select "(CALL == 'OUT')" --select-type-to-include INDEL

        UNK_SNP="$(gatk --java-options "-Xmx~{memoryJava}G" CountVariants -V selected.unk.snp.vcf | tail -1)"
        UNK_INDEL="$(gatk --java-options "-Xmx~{memoryJava}G" CountVariants -V selected.unk.indel.vcf | tail -1)"

        echo "$UNK_SNP" > unk_snp.txt
        echo "$UNK_INDEL" > unk_indel.txt
    >>>
    runtime {
                docker: "us.gcr.io/broad-gatk/gatk:"+gatkTag
                preemptible: select_first([preemptible,0])
                disks: "local-disk " + disk_size + " HDD"
                bootDiskSizeGb: "16"
                memory: memoryRam + " GB"
    }
    output {
        Int UNK_SNP = read_int("unk_snp.txt")
        Int UNK_INDEL = read_int("unk_indel.txt")
    }
}

#Count number of variants which were outside confidence region based on GenotypeConcordance annotated vcf
task CountUNKGC {
        input {
            File? vcfAnnotated=""
            File? vcfIndexAnnotated=""
            File vcfOrig
            File vcfIndexOrig
            Int? preemptible
            Int? memoryMaybe
            File? stratIL
            String gatkTag
        }
        Int memoryDefault = 16
        Int memoryJava = select_first([memoryMaybe,memoryDefault])
        Int memoryRam = memoryJava+2
        Int disk_size = 10 + ceil(size(vcfAnnotated, "GB") + size(vcfIndexAnnotated, "GB") + 2 * size(vcfOrig, "GB") + size(vcfIndexOrig, "GB") + size(stratIL, "GB"))

        command <<<
                set -xeuo pipefail

                gatk --java-options "-Xmx~{memoryJava}G" SelectVariants -V ~{vcfOrig} -O selected.unk.snp.vcf.gz ~{"-L "+stratIL} --select-type-to-include SNP --discordance ~{vcfAnnotated}
                gatk --java-options "-Xmx~{memoryJava}G" SelectVariants -V ~{vcfOrig} -O selected.unk.indel.vcf.gz ~{"-L "+stratIL} --select-type-to-include INDEL --discordance ~{vcfAnnotated}

                UNK_SNP="$(gatk --java-options "-Xmx~{memoryJava}G" CountVariants -V selected.unk.snp.vcf.gz | tail -1)"
                UNK_INDEL="$(gatk --java-options "-Xmx~{memoryJava}G" CountVariants -V selected.unk.indel.vcf.gz | tail -1)"

                echo "$UNK_SNP" > unk_snp.txt
                echo "$UNK_INDEL" > unk_indel.txt
            >>>


        runtime {
                    docker: "us.gcr.io/broad-gatk/gatk:"+gatkTag
                    preemptible: select_first([preemptible,0])
                    disks: "local-disk " + disk_size + " HDD"
                    bootDiskSizeGb: "16"
                    memory: memoryRam + " GB"
        }

        output {
            Int UNK_SNP = read_int("unk_snp.txt")
            Int UNK_INDEL = read_int("unk_indel.txt")
        }
}

#Count TP,FP,FN for a particular selection of variants given by jexl
task EvalForVariantSelection {
    input {
        File? vcf=""
        File? vcfIndex=""
        String jexl
        Int? preemptible
        Int? memoryMaybe
        String engine
        String selectTPCall
        String selectTPBase
        String selectFN
        String selectFP
        String sampleCall
        String sampleBase
        String gatkTag

        Array[String]? annotationNames=[]
        File? gatkJarForAnnotation
        File reference
        File refDict
        File refIndex
    }

    Int memoryDefault = 16
    Int memoryJava = select_first([memoryMaybe,memoryDefault])
    Int memoryRam = memoryJava+2

    String selectionTPCall = jexl + " && " + selectTPCall
    String selectionTPBase = jexl + " && " + selectTPBase
    String selectionFN = jexl + " && " + selectFN
    String selectionFP = jexl + " && " + selectFP

    Int disk_size = 10 + ceil(4.2 * size(vcf, "GB") + 2.2 * size(vcfIndex, "GB") + size(reference, "GB"))

    command <<<
        set -xeuo pipefail

        VCF=~{vcf}
        if [[ ! -z "~{gatkJarForAnnotation}" ]]; then
            java -jar ~{gatkJarForAnnotation} VariantAnnotator -V ~{vcf} -O annotated.vcf.gz ~{true="-A" false="" length(annotationNames)>0} ~{sep=" -A " annotationNames} -R ~{reference}
            VCF = annotated.vcf.gz
        else
            touch annotated.vcf.gz
        fi

        gatk --java-options "-Xmx~{memoryJava}G" SelectVariants -V $VCF -O selected.TP_CALL.vcf.gz -select "~{selectionTPCall}" -sn ~{sampleCall}
        gatk --java-options "-Xmx~{memoryJava}G" SelectVariants -V $VCF -O selected.TP_BASE.vcf.gz -select "~{selectionTPBase}" -sn ~{sampleBase}
        gatk --java-options "-Xmx~{memoryJava}G" SelectVariants -V $VCF -O selected.FN.vcf.gz -select "~{selectionFN}" -sn ~{sampleBase}
        gatk --java-options "-Xmx~{memoryJava}G" SelectVariants -V $VCF -O selected.FP.vcf.gz -select "~{selectionFP}" -sn ~{sampleCall}

        TP_CALL="$(gatk --java-options "-Xmx~{memoryJava}G" CountVariants -V selected.TP_CALL.vcf.gz | tail -1)"
        TP_BASE="$(gatk --java-options "-Xmx~{memoryJava}G" CountVariants -V selected.TP_BASE.vcf.gz | tail -1)"
        FN="$(gatk --java-options "-Xmx~{memoryJava}G" CountVariants -V selected.FN.vcf.gz | tail -1)"
        FP="$(gatk --java-options "-Xmx~{memoryJava}G" CountVariants -V selected.FP.vcf.gz | tail -1)"

        echo "$TP_CALL" > tp_call.txt
        echo "$TP_BASE" > tp_base.txt
        echo "$FN" > fn.txt
        echo "$FP" > fp.txt
    >>>

    runtime {
            docker: "us.gcr.io/broad-gatk/gatk:"+gatkTag
            preemptible: select_first([preemptible,0])
            disks: "local-disk " + disk_size + " HDD"
            bootDiskSizeGb: "16"
            memory: memoryRam + " GB"
    }

    output {
        Int TP_CALL = read_int("tp_call.txt")
        Int TP_BASE = read_int("tp_base.txt")
        Int FP = read_int("fp.txt")
        Int FN = read_int("fn.txt")
        File selectedTPCall="selected.TP_CALL.vcf.gz"
        File selectedTPBase="selected.TP_BASE.vcf.gz"
        File selectedFP="selected.FP.vcf.gz"
        File selectedFN="selected.FN.vcf.gz"

        File annotated="annotated.vcf.gz"
        File selectedTPCallIndex="selected.TP_CALL.vcf.gz.tbi"
        File selectedTPBaseIndex="selected.TP_BASE.vcf.gz.tbi"
        File selectedFPIndex="selected.FP.vcf.gz.tbi"
        File selectedFNIndex="selected.FN.vcf.gz.tbi"
    }

}

#create csv file of statistics based on TP,FP,FN
task SummariseForIndelSelection {
    input {
        String evalLabel
        String truthLabel
        String? stratLabel
        String indelLabel
        String engine
        String igvSession
        Int TP_CALL
        Int TP_BASE
        Int FP
        Int FN
        Int? preemptible

    }

    command <<<
        set -xeuo pipefail

        Rscript -<<"EOF" ~{TP_CALL} ~{TP_BASE} ~{FN} ~{FP} ~{evalLabel} ~{truthLabel} ~{indelLabel} ~{engine} ~{default="" stratLabel} ~{igvSession}
        GetSelectionValue<-function(name, target) {

                  if(target=="insertion" || target=="deletion") {
                    return(NA)
                  }
                  pos_start<-regexpr(target,name)
                  sub = substring(name,pos_start+attr(pos_start,"match.length")+1,nchar(name))
                  split_sub = strsplit(sub,"_")
                  val = if(grepl("^m",split_sub[[1]][[1]])) -as.double(gsub("m","",split_sub[[1]][[1]])) else as.double(split_sub[[1]][[1]])
                }

        args <-commandArgs(trailingOnly = TRUE)
        indel_options <-c("deletion","insertion","indel_fine","indel_coarse")
        indel_type <- mapply(grepl,indel_options,args[7])
        indel_type <- indel_options[indel_type[indel_options]]
        indel_length <- GetSelectionValue(args[7],indel_type)
        if (length(args)<10) {
          stratifier <- NA
        } else {
          stratifier <- args[9]
        }
        table <- data.frame("Name"=args[5], "Truth_Set"=args[6],"Comparison_Engine"=args[8],"Stratifier"=stratifier,
                            "IndelLength"= indel_length,
                            "Recall"=as.numeric(args[2])/(as.numeric(args[2])+as.numeric(args[3])),"Precision"=as.numeric(args[1])/(as.numeric(args[1])+as.numeric(args[4])),"TP_Base"=as.numeric(args[2]),"TP_Eval"=as.numeric(args[1]),
                            "FP"=as.numeric(args[4]),"FN"=as.numeric(args[3]),"IGV_Session"=args[length(args)],"Summary_Type"=indel_type)
        table$F1_Score <- 2*table$Precision*table$Recall/(table$Precision+table$Recall)
        write.csv(table,paste(args[8],".",indel_type,".summary.csv",sep=""),row.names = FALSE)
        EOF
    >>>

    runtime {
            docker: "rocker/tidyverse"
            preemptible: select_first([preemptible,0])
            disks: "local-disk 10 HDD"
        }

    output {
        File summaryOut = glob("*.summary.csv")[0]
    }

}

#create csv file of statistics based on TP,FP,FN
task SummariseForVariantSelection {
    input {
        String evalLabel
        String truthLabel
        String? stratLabel
        String variantLabel
        String engine
        String igvSession
        Int TP_CALL
        Int TP_BASE
        Int FP
        Int FN
        Int? preemptible

    }

    command <<<
        set -xeuo pipefail

        Rscript -<<"EOF" ~{TP_CALL} ~{TP_BASE} ~{FN} ~{FP} ~{evalLabel} ~{truthLabel} ~{variantLabel} ~{engine} ~{default="" stratLabel} ~{igvSession}
        args <-commandArgs(trailingOnly = TRUE)
        if (length(args)<10) {
          stratifier <- NA
        } else {
          stratifier <- args[9]
        }
        table <- data.frame("Name"=args[5], "Truth_Set"=args[6],"Comparison_Engine"=args[8],"Stratifier"=stratifier,
                            "IndelLength"= NA,
                            "Recall"=as.numeric(args[2])/(as.numeric(args[2])+as.numeric(args[3])),"Precision"=as.numeric(args[1])/(as.numeric(args[1])+as.numeric(args[4])),"TP_Base"=as.numeric(args[2]),"TP_Eval"=as.numeric(args[1]),
                            "FP"=as.numeric(args[4]),"FN"=as.numeric(args[3]),"IGV_Session"=args[length(args)],"Summary_Type"="summary","Type"=args[7])
        table$F1_Score <- 2*table$Precision*table$Recall/(table$Precision+table$Recall)
        write.csv(table,paste(args[8],".",args[7],".summary.csv",sep=""),row.names = FALSE)
        EOF
    >>>

    runtime {
            docker: "rocker/tidyverse"
            preemptible: select_first([preemptible,0])
            disks: "local-disk 10 HDD"
        }

    output {
        File summaryOut = glob("*.summary.csv")[0]
    }

}

#Convert vcfeval output statistics to final output format
task SummariseVcfEval {
    input {
        String evalLabel
        String truthLabel
        String areVariants
        String? igvSession
        String? stratLabel
        File? summaryFile
        Int? unkSNP
        Int? unkINDEL
        Int? preemptible
    }
    Int disk_size = 10 + ceil(2.2 * size(summaryFile, "GB"))

    command <<<
        set -xeuo pipefail

        Rscript -<<"EOF" ~{evalLabel} ~{truthLabel} ~{default="" summaryFile} ~{default="" stratLabel} ~{default="" igvSession} ~{default="" unkSNP} ~{default="" unkINDEL} ~{areVariants}
        args <- commandArgs(trailingOnly = TRUE)
        if (args[length(args)]=="yes") {
            table_vcfeval <- read.csv(args[3])
            if (length(args)==7) {
                table_vcfeval$Stratifier <- NA
                table_vcfeval$IGV_Session <- args[4]
                table_vcfeval$UNK[table_vcfeval$Type=="SNP"]=args[5]
                table_vcfeval$UNK[table_vcfeval$Type=="INDEL"]=args[6]
            } else {
                table_vcfeval$Stratifier <- args[4]
                table_vcfeval$IGV_Session <- args[5]
                table_vcfeval$UNK[table_vcfeval$Type=="SNP"]=args[6]
                table_vcfeval$UNK[table_vcfeval$Type=="INDEL"]=args[7]
            }
            write(ifelse(is.na(table_vcfeval[table_vcfeval$Type=="SNP",  ]$Precision), 0, table_vcfeval[table_vcfeval$Type=="SNP",  ]$Precision), "snpPrecision.txt")
            write(ifelse(is.na(table_vcfeval[table_vcfeval$Type=="INDEL",]$Precision), 0, table_vcfeval[table_vcfeval$Type=="INDEL",]$Precision), "indelPrecision.txt")
            write(ifelse(is.na(table_vcfeval[table_vcfeval$Type=="SNP",  ]$Recall),    0, table_vcfeval[table_vcfeval$Type=="SNP",  ]$Recall),    "snpRecall.txt")
            write(ifelse(is.na(table_vcfeval[table_vcfeval$Type=="INDEL",]$Recall),    0, table_vcfeval[table_vcfeval$Type=="INDEL",]$Recall),    "indelRecall.txt")
            write(ifelse(is.na(table_vcfeval[table_vcfeval$Type=="SNP",  ]$F1_Score),  0, table_vcfeval[table_vcfeval$Type=="SNP",  ]$F1_Score),  "snpF1Score.txt")
            write(ifelse(is.na(table_vcfeval[table_vcfeval$Type=="INDEL",]$F1_Score),  0, table_vcfeval[table_vcfeval$Type=="INDEL",]$F1_Score),  "indelF1Score.txt")
        } else {
            types <- c("INDEL","SNP")
            recall <- c(NA,NA)
            precision <- c(NA,NA)
            f1_score <- c(NA,NA)
            tp <- c(0,0)
            fp <- c(0,0)
            fn <- c(0,0)
            unk <- c(0,0)
            igv_session <- c(NA,NA)
            table_vcfeval <- data.frame("Type"=types,"Recall"=recall,"Precision"=precision,"F1_Score"=f1_score,"TP_Base"=tp,"TP_Eval"=tp,"FP"=fp,"FN"=fn,"UNK"=unk,"IGV_Session"=igv_session)
            if (length(args)==3) {
                table_vcfeval$Stratifier <- NA
            } else {
                table_vcfeval$Stratifier <- args[3]
            }
            write(0,"snpPrecision.txt")
            write(0,"indelPrecision.txt")
            write(0,"snpRecall.txt")
            write(0,"indelRecall.txt")
            write(0,"snpF1Score.txt")
            write(0,"indelF1Score.txt")
        }
        table_vcfeval$Name <- args[1]
        table_vcfeval$Truth_Set <- args[2]
        table_vcfeval$Summary_Type <- "summary"
        table_vcfeval$Comparison_Engine <-"VcfEval"
        write.csv(table_vcfeval,"vcfeval.summary.csv",row.names = FALSE)
        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse"
        preemptible: select_first([preemptible,0])
        disks: "local-disk " + disk_size + " HDD"
    }

    output{
        File summaryOut="vcfeval.summary.csv"
        Float snpPrecision = read_float("snpPrecision.txt")
        Float indelPrecision = read_float("indelPrecision.txt")
        Float snpRecall = read_float("snpRecall.txt")
        Float indelRecall = read_float("indelRecall.txt")
        Float snpF1Score = read_float("snpF1Score.txt")
        Float indelF1Score = read_float("indelF1Score.txt")
    }
}

#Convert hap.py output statistics to final output format
task SummariseHappy {
     input {
        String evalLabel
        String truthLabel
        String areVariants
        String? stratLabel
        File? summaryFile
        Int? preemptible
        String? igvSession
     }

     Int disk_size = 10 + ceil(2.2 * size(summaryFile, "GB"))

    command <<<
        Rscript -<<"EOF" ~{evalLabel} ~{truthLabel} ~{default="" summaryFile} ~{default="" stratLabel} ~{default="" igvSession} ~{areVariants}
        args <-commandArgs(trailingOnly = TRUE)
        if (args[length(args)]=="yes") {
          table_happy <- read.csv(args[3])
          colnames(table_happy)[10]<-"Recall"
          colnames(table_happy)[11]<-"Precision"
          colnames(table_happy)[13]<-"F1_Score"
          colnames(table_happy)[4]<-"TP_Base"
          table_happy$TP_Eval = table_happy$TP_Base
          colnames(table_happy)[5]<-"FN"
          colnames(table_happy)[7]<-"FP"
          colnames(table_happy)[8]<-"UNK"
          if (length(args)==5) {
            table_happy$Stratifier <- NA
            table_happy$IGV_Session <- args[4]
          } else {
            table_happy$Stratifier <-args[4]
            table_happy$IGV_Session <- args[5]
          }


        } else {
          types <- c("INDEL","SNP")
          recall <- c(NA,NA)
          precision <- c(NA,NA)
          f1_score <- c(NA,NA)
          igv_session <- c(NA,NA)
          tp <- c(0,0)
          fp <- c(0,0)
          fn <- c(0,0)
          unk <- c(0,0)
          filters <- c("PASS","PASS")
          table_happy <- data.frame("Type"=types,"Recall"=recall,"Precision"=precision,"F1_Score"=f1_score,"TP_Base"=tp,"TP_Eval"=tp,"FP"=fp,"FN"=fn,"UNK"=unk,"Filter"=filters,"IGV_Session"=igv_session)
          if (length(args)==3) {
            table_happy$Stratifier <- NA
          } else {
            table_happy$Stratifier <-args[3]
          }

        }
        table_happy$Name <- args[1]
        table_happy$Truth_Set <- args[2]
        table_happy$Comparison_Engine <-"Happy"
        table_happy$Summary_Type <- "summary"
        table_happy<-table_happy[table_happy$Filter=="PASS",c("Type","Recall","Precision","Name","Truth_Set","Comparison_Engine","F1_Score","TP_Eval","TP_Base","FP","FN","UNK","Stratifier","IGV_Session","Summary_Type")]
        if (nrow(subset(table_happy,Type=="INDEL"))==0) {
          table_add <- data.frame("Type"="INDEL","Recall"=NA,"Precision"=NA,"F1_Score"=NA,"TP_Eval"=0,"TP_Base"=0,"FP"=0,"FN"=0,"UNK"=0,"Name"=args[1],"Truth_Set"=args[2],"Comparison_Engine"="Happy","Stratifier"=table_happy$Stratifier[1],"IGV_Session"=NA,"Summary_Type"="summary")
          table_happy <- rbind(table_happy,table_add)
        }
        write.csv(table_happy,"happy.summary.csv",row.names = FALSE)
        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse"
        preemptible: select_first([preemptible,0])
        disks: "local-disk " + disk_size + " HDD"
    }

    output{
        File summaryOut="happy.summary.csv"
    }
}

#Convert GenotypeConcordance output statistics to final output format
task SummariseGATKGC {
    input {
        String evalLabel
        String truthLabel
        String areVariants
        String? stratLabel
        File? summaryFile
        File? summaryCounts
        Int? unkSNP
        Int? unkINDEL
        Int? preemptible
        String? igvSession
    }

    Int disk_size = 10 + ceil(2.2 * size(summaryFile, "GB") + 2.2 * size(summaryCounts, "GB"))

    command <<<
        set -xeuo pipefail

        Rscript -<<"EOF" ~{evalLabel} ~{truthLabel} ~{default="" summaryFile} ~{default="" summaryCounts} ~{default="" stratLabel} ~{default="" igvSession} ~{default="" unkSNP} ~{default="" unkINDEL} ~{areVariants}
        args <- commandArgs(trailingOnly = TRUE)
        if (args[length(args)]=="yes") {
            table_GC <- read.table(args[3],skip = 6, header = TRUE,sep="\t",na.strings="?")
            table_counts_GC <- read.table(args[4],skip = 6, header = TRUE,sep="\t",na.strings="?")
            table_GC$F1_Score <- 2*(table_GC$VAR_PPV*table_GC$VAR_SENSITIVITY)/(table_GC$VAR_PPV+table_GC$VAR_SENSITIVITY)
            table_GC$TP_Eval <- table_counts_GC$TP_COUNT
            table_GC$TP_Base <- table_counts_GC$TP_COUNT
            table_GC$FP <- table_counts_GC$FP_COUNT
            table_GC$FN <- table_counts_GC$FN_COUNT
            colnames(table_GC)[10]<-"Recall"
            colnames(table_GC)[11]<-"Precision"
            colnames(table_GC)[1]<-"Type"
            if (length(args)==8) {
                table_GC$Stratifier <- NA
                table_GC$IGV_Session <- args[5]
                table_GC$UNK[table_GC$Type=="SNP"]=args[6]
                table_GC$UNK[table_GC$Type=="INDEL"]=args[7]

            } else {
                table_GC$Stratifier <- args[5]
                table_GC$IGV_Session <- args[6]
                table_GC$UNK[table_GC$Type=="SNP"]=args[7]
                table_GC$UNK[table_GC$Type=="INDEL"]=args[8]
            }
        } else {
            types <- c("INDEL","SNP")
            recall <- c(NA,NA)
            precision <- c(NA,NA)
            f1_score <- c(NA,NA)
            tp <- c(0,0)
            fp <- c(0,0)
            fn <- c(0,0)
            unk <- c(0,0)
            igv_session <- c(NA,NA)
            table_GC <- data.frame("Type"=types,"Recall"=recall,"Precision"=precision,"F1_Score"=f1_score,"TP_Eval"=tp,"TP_Base"=tp,"FP"=fp,"FN"=fn,"UNK"=unk,"IGV_Session"=igv_session)
            if (length(args)==3) {
                table_GC$Stratifier <- NA
            } else {
                table_GC$Stratifier <- args[3]
            }
        }
        table_GC$Name <- args[1]
        table_GC$Truth_Set <- args[2]
        table_GC$Comparison_Engine <-"GATK_GC"
        table_GC$Summary_Type <- "summary"
        table_GC <- table_GC[,c("Type","Recall","Precision","Name","Truth_Set","Comparison_Engine","F1_Score","TP_Eval","TP_Base","FP","FN","UNK","Stratifier","IGV_Session","Summary_Type")]
        write.csv(table_GC,"gatkgc.summary.csv",row.names = FALSE)
        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse"
        preemptible: select_first([preemptible,0])
        disks: "local-disk " + disk_size + " HDD"
    }

    output{
        File summaryOut="gatkgc.summary.csv"
    }
}

#Combine summaries from multiple csv into a single csv
task CombineSummaries {
    input {
        Array[File] summaries
        Int? preemptible
    }
    String dollar="$"

    Int disk_size = 10 + ceil(2 * size(summaries, "GB"))
    command <<<
        set -xeuo pipefail

        Rscript -<<"EOF"
        library(readr)
        library(dplyr)
        library(purrr)
        summary_files <- read_csv("~{write_lines(summaries)}", col_names = FALSE)
        merged<- as.list(summary_files$X1) %>% map(read_csv) %>% reduce(bind_rows)
        write.csv(merged,"summary.csv",row.names = FALSE)
        EOF
    >>>

    runtime {
            docker: "rocker/tidyverse"
            preemptible: select_first([preemptible,0])
            disks: "local-disk " + disk_size + " HDD"
        }

    output{
        File summaryOut="summary.csv"
    }
}

#Use CrosscheckFingerprints to match evaluation vcfs to appropriate truth vcfs
task MatchEvalTruth {
    input{
        File evalVcf
        File truthVcf
        File evalVcfIndex
        File truthVcfIndex
        File hapMap
        Int? preemptible
        Int? memoryMaybe
        String gatkTag
    }
    Int memoryDefault = 16
    Int memoryJava = select_first([memoryMaybe,memoryDefault])
    Int memoryRam = memoryJava+2
    Int disk_size = 10 + ceil(size(hapMap, "GB") + size(evalVcf, "GB") + size(evalVcfIndex, "GB") + size(truthVcf, "GB") + size(truthVcfIndex, "GB"))

    command <<<
        gatk --java-options "-Xmx~{memoryJava}G" CrosscheckFingerprints -I ~{evalVcf} -SI ~{truthVcf} -H ~{hapMap} --CROSSCHECK_MODE CHECK_ALL_OTHERS --CROSSCHECK_BY FILE --EXPECT_ALL_GROUPS_TO_MATCH
    >>>

    runtime {
            docker: "us.gcr.io/broad-gatk/gatk:"+gatkTag
            preemptible: select_first([preemptible,0])
            disks: "local-disk " + disk_size + " HDD"
            bootDiskSizeGb: "16"
            memory: memoryRam + " GB"
        }
}

# creates an IGV session
# given a list of IGV compatible file paths
task WriteXMLfile {
    input {
        Array[String] input_files
        String reference_version
        String file_name

        Array[String]? input_names
    }

    Array[String] input_names_prefix = if defined(input_names) then prefix('-n ', select_first([input_names])) else []

    command <<<
        set -euxo pipefail

        # because of some nonsense above, we need to play some bash tricks here to make sure that
        # the inputs are labeled correctly:
        echo '~{sep=" " input_names_prefix}' | sed -e 's#-n[ \t]*-n#-n#g' -e 's#-n[ \t]*$##' > labels.txt

        bash /usr/writeIGV.sh ~{reference_version} ~{sep=" " input_files} $(cat labels.txt) > "~{file_name}.xml"
    >>>
    runtime {
        docker: "quay.io/mduran/generate-igv-session_2:v1.0"
    }
    output {
        File igv_session = "${file_name}.xml"
    }
}

task CreateIntervalList{
    input {
        File reference
        File reference_index
        File reference_dict
        String interval_string
        String gatkTag
        String? dummyInputForTerraCallCaching
    }
    command {
        gatk PreprocessIntervals \
        -R ~{reference} \
        -L ~{interval_string} \
        -O output.interval_list \
        --bin-length 0 \
        -imr OVERLAPPING_ONLY \
        -padding 0
    }
    output {
        File interval_list = "output.interval_list"
    }
    runtime {
        preemptible: 3
        docker: "us.gcr.io/broad-gatk/gatk:"+gatkTag
        disks: "local-disk 100 HDD"
        memory: "4 GB"
    }
}

#Print given message to stderr and return an error
task ErrorWithMessage{
    input {
        String message
    }
    command <<<
    >&2 echo "Error: ~{message}"
    exit 1
    >>>

    runtime {
        docker: "ubuntu"
    }
}
