version 1.0

import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow RemoveRedundantReblocksInGVCF {
    meta {
        description:
        "A posthoc fix for some gVCF where the concatenation introduces redundant ref blocks"
    }

    parameter_meta {
        original_shard_number:
        "If the original gVCF were produced by our standard pipeline with standard config, then this value is 21 for GRCh38 and 13 for T2T."
    }

    input {
        File input_gvcf
        File input_gvcf_tbi
        Int original_shard_number
        File ref_map_file

        String gcs_out_root_dir
    }

    output {
        File fixed_gvcf     = SaveGVCF.gcs_path
        File fixed_gvcf_tbi = SaveGVCFtbi.gcs_path
        File nonref_vcf     = SaveVCF.gcs_path
        File nonref_vcf_tbi = SaveVCFtbi.gcs_path
    }

    Map[String, String] ref_map = read_map(ref_map_file)
    String output_dir = sub(gcs_out_root_dir, "/$", "") + "/RemoveRedundantReblocksInGVCF"

    Array[String] canonical_chromosomes = [
                                            'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
                                            'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                                            'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
                                            'chrX', 'chrY', 'chrM'
                                          ]
    scatter (chromosome in canonical_chromosomes){
        call VU.SubsetVCF { input:
            vcf_gz  = input_gvcf,
            vcf_tbi = input_gvcf_tbi,
            locus   = chromosome,
            prefix  = basename(input_gvcf, ".g.vcf.gz") + ".~{chromosome}.g",
            omit_cmd_in_header = true
        }
        call FixAChromosome { input:
            input_gvcf = SubsetVCF.subset_vcf,
            ref_fasta = ref_map['fasta'],
            ref_fasta_fai = ref_map['fai'],
            original_shard_number = original_shard_number,
        }
    }
    call VU.MergeHumanPrimaryChromosomesVCFs as MergeGVCF { input:
        vcfs = FixAChromosome.fixed_gvcf,
        out_prefix = basename(input_gvcf, ".g.vcf.gz") + ".refblock_dedup.g"
    }
    call FF.FinalizeToFile as SaveGVCF { input:
        file = MergeGVCF.vcf,
        outdir = output_dir
    }
    call FF.FinalizeToFile as SaveGVCFtbi { input:
        file = MergeGVCF.tbi,
        outdir = output_dir
    }
    call VU.MergeHumanPrimaryChromosomesVCFs as MergeVCF { input:
        vcfs = FixAChromosome.extracted_vcf,
        out_prefix = basename(input_gvcf, ".g.vcf.gz")
    }
    call FF.FinalizeToFile as SaveVCF { input:
        file = MergeVCF.vcf,
        outdir = output_dir
    }
    call FF.FinalizeToFile as SaveVCFtbi { input:
        file = MergeVCF.tbi,
        outdir = output_dir
    }
}

task FixAChromosome {
    meta {
        description:
        "A posthoc fix for some gVCF where the concatenation introduces redundant ref blocks"
    }

    parameter_meta {
        original_shard_number:
        "If the original gVCF were produced by our standard pipeline with standard config, then this value is 21 for GRCh38 and 13 for T2T."
    }

    input {
        File input_gvcf
        File ref_fasta
        File ref_fasta_fai
        Int original_shard_number
    }

    output {
        File fixed_gvcf    = "~{prefix}.cleanedup.g.vcf.gz"
        File extracted_vcf = "~{prefix}.cleanedup.converted.non_ref.vcf.gz"

        # File? original_gt_stats   = "~{prefix}.original.GT.stats"
        # File? fixed_gvcf_gt_stats = "~{prefix}.fixedgVCF.GT.stats"
        # File? vcf_gt_stats        = "~{prefix}.nonRefVCF.GT.stats"
    }

    String prefix = basename(input_gvcf, ".g.vcf.gz")

    command <<<
    set -euxo pipefail

        ####################################################################
        # sometimes, there are no variants in the input gVCF (i.e. chrY on female samples)
        ####################################################################
        var_cnt=$(bcftools view -H ~{input_gvcf} | wc -l | awk '{print $1}')
        if [[ "${var_cnt}" -eq 0 ]]; then
            if [[ "~{prefix}" != *'chrY' && "~{prefix}" != *'chrM' ]]; then
                echo "No variants detected in the whole chromosome (except chrY/chrM), something is wrong."
                exit 1
            fi
            echo "No variants in the input gVCF, hence nothing to fix"
            cp ~{input_gvcf} \
               "~{prefix}.cleanedup.g.vcf.gz"

            echo "convert to VCF..."
            time \
            bcftools convert \
                --fasta-ref "~{ref_fasta}" \
                --gvcf2vcf \
                --write-index -O 'z' \
                -o "~{prefix}.cleanedup.converted.vcf.gz" \
                "~{prefix}.cleanedup.g.vcf.gz"
            echo "done"

            echo "filter away REF sites..."
            time \
            bcftools view -O 'z' \
                -i 'GT[*]="alt"' \
                -o "~{prefix}.cleanedup.converted.non_ref.vcf.gz" \
                "~{prefix}.cleanedup.converted.vcf.gz"
            echo "done"

            exit 0
        fi

        ####################################################################
        # for GRCh38, we shard by 22 shards, including the alt shard
        # but since DV 1.6 has issues with the alt shard, we ran with 21 shards
        # therefore,
        # if a line occurs 21 times (21-mers), it means that block is a true ref ref lock,
        #                                      so we want to keep just one copy
        # if a line occurs 20 times (20-mers), it means that block is a ref block on shards that don't have reads in that block,
        #                                      so we want to drop such lines
        # if a line occurs  1 time,  we want to keep it
        # if a line occurs other times, error out

        # for T2T, we shard by 13 shards, and no shard is skipped.
        # so change them accrodingly
        ####################################################################
        KEEP_ONE_COPY_OCCUR=~{original_shard_number}
        KEEP_ZERO_COPY_OCCUR=$(( KEEP_ONE_COPY_OCCUR - 1 ))

        ####################################################################
        # first gather the lines that hold the signature, also count appearances
        ####################################################################
        echo "extracting lines matching that pattern...."

        time \
        bcftools view -H "~{input_gvcf}" \
            | grep $'GT:GQ:MIN_DP:PL\t0/0:1:0:0,0,0$' \
            | sort \
            | uniq -c \
            | sed -E "s/^\s+//" \
        > tmp.suspect.duplicate.blocks.count

        echo "done extracting"

        awk -F ' ' '{print $1}' tmp.suspect.duplicate.blocks.count \
            | sort \
            | uniq -c \
            | awk '{print $2}' \
            | sort \
        > tmp.occurances.txt

        set +e
        should_be_zero=$(grep -vw "1" tmp.occurances.txt | grep -vw "${KEEP_ZERO_COPY_OCCUR}" | grep -vw "${KEEP_ONE_COPY_OCCUR}" | wc -l | awk '{print $1}')
        set -e
        if [[ "${should_be_zero}" -ne 0 ]]; then
            echo "num. of occurances of that pattern should be one of (1, ${KEEP_ZERO_COPY_OCCUR}, ${KEEP_ONE_COPY_OCCUR})."
            cat tmp.occurances.txt
            exit 1
        fi

        ####################################################################
        # extract the 20-mers and 21-mers separately, and sort them
        ####################################################################
        touch "tmp.suspect.duplicate.${KEEP_ONE_COPY_OCCUR}mers.txt"
        if grep -qE "^${KEEP_ONE_COPY_OCCUR}\s+" tmp.suspect.duplicate.blocks.count; then
            grep -E "^${KEEP_ONE_COPY_OCCUR}\s+" tmp.suspect.duplicate.blocks.count \
                | sed -E "s/^${KEEP_ONE_COPY_OCCUR}\s+//" \
                | sort \
            > "tmp.suspect.duplicate.${KEEP_ONE_COPY_OCCUR}mers.txt"
        fi

        touch "tmp.suspect.duplicate.${KEEP_ZERO_COPY_OCCUR}mers.txt"
        if grep -qE "^${KEEP_ZERO_COPY_OCCUR}\s+" tmp.suspect.duplicate.blocks.count; then
            grep -E "^${KEEP_ZERO_COPY_OCCUR}\s+" tmp.suspect.duplicate.blocks.count \
                | sed -E "s/^${KEEP_ZERO_COPY_OCCUR}\s+//" \
                | sort \
            > "tmp.suspect.duplicate.${KEEP_ZERO_COPY_OCCUR}mers.txt"
        fi

        ####################################################################
        # sort the VCF body and the file holding the 20-mers, then comm to remove those 20-mers lines
        ####################################################################
        cat "tmp.suspect.duplicate.${KEEP_ONE_COPY_OCCUR}mers.txt" \
            "tmp.suspect.duplicate.${KEEP_ZERO_COPY_OCCUR}mers.txt" \
            | sort \
        > tmp.remove.txt

        echo "extracting gVCF body and lexico sorting, just wait...."
        time \
        bcftools view -H "~{input_gvcf}" \
            | sort \
            | uniq \
        > tmp.sorted.body

        echo "done extration & lexico sorting, remove false and redundant Ref blocks...."
        time \
        comm -23 tmp.sorted.body tmp.remove.txt \
        > "tmp.sorted.sans-${KEEP_ZERO_COPY_OCCUR}-${KEEP_ONE_COPY_OCCUR}.body"
        rm tmp.sorted.body
        echo "done removal"

        ####################################################################
        # add back just one copy of the 21-mers, then output indexed gVCF
        ####################################################################
        zgrep "^#" "~{input_gvcf}" > header.txt
        echo "final step: coordinate sort the deduplicated gVCF...."
        time \
        cat header.txt \
            "tmp.suspect.duplicate.${KEEP_ONE_COPY_OCCUR}mers.txt" \
            "tmp.sorted.sans-${KEEP_ZERO_COPY_OCCUR}-${KEEP_ONE_COPY_OCCUR}.body" \
        | bcftools sort \
            -m 10G -O 'z' \
            -o "~{prefix}.cleanedup.g.vcf.gz" \
            -
        echo "done"

        echo "tabix...."
        time
        tabix -p vcf \
            "~{prefix}.cleanedup.g.vcf.gz"
        echo "done"

        ls -lh tmp*

        ####################################################################
        # finally, convert to VCF and filter away the 0/0 sites
        ####################################################################
        echo "convert to VCF"
        time \
        bcftools convert \
            --fasta-ref "~{ref_fasta}" \
            --gvcf2vcf \
            --write-index -O 'z' \
            -o "~{prefix}.cleanedup.converted.vcf.gz" \
            "~{prefix}.cleanedup.g.vcf.gz"

        echo "filter REF sites"
        time \
        bcftools view -O 'z' \
            -i 'GT[*]="alt"' \
            -o "~{prefix}.cleanedup.converted.non_ref.vcf.gz" \
            "~{prefix}.cleanedup.converted.vcf.gz"

        ####################################################################
        # QC on GT
        ####################################################################
        bcftools view -H "~{input_gvcf}" \
            | awk -F '\t' '{print $10}' \
            | awk -F':' '{print $1}' \
            | sort | uniq -c \
        > ~{prefix}.original.GT.stats &

        bcftools view -H "~{prefix}.cleanedup.g.vcf.gz" \
            | awk -F '\t' '{print $10}' \
            | awk -F':' '{print $1}' \
            | sort | uniq -c \
        > ~{prefix}.fixedgVCF.GT.stats &

        bcftools view -H "~{prefix}.cleanedup.converted.non_ref.vcf.gz" \
            | awk -F '\t' '{print $10}' \
            | awk -F':' '{print $1}' \
            | sort | uniq -c \
        > ~{prefix}.nonRefVCF.GT.stats &

        wait
        tail -n20 *GT.stats

        touch should_be_empty.txt
        diff <(grep -vF "0/0" ~{prefix}.original.GT.stats) \
                <(grep -vF "0/0" ~{prefix}.fixedgVCF.GT.stats) \
        > should_be_empty.txt

        touch should_also_be_empty.txt
        diff ~{prefix}.nonRefVCF.GT.stats \
                <(grep -vF "0/0" ~{prefix}.fixedgVCF.GT.stats | grep -vF "./.") \
        > should_also_be_empty.txt

        if [[ (-s should_be_empty.txt)  || (-s should_also_be_empty.txt) ]]; then
            echo "The fixing procedure messed up GT!!!"
            exit 1
        fi
    >>>

    runtime {
        cpu:            2
        memory:         "10 GiB"
        disks:          "local-disk " + (10 + ceil(size([input_gvcf, ref_fasta], "GiB"))) + " SSD"
        preemptible:    0
        maxRetries:     0
        docker:         "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}
