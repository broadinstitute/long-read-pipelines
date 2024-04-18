version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Utility/Finalize.wdl" as FF

import "../Utility/ShardWholeGenome.wdl"


workflow SplitCohortBAMsAndVCFsForPhasing {
    meta {
        desciption:
        "Preprocess the WGS BAMs and joint-call small variant VCF and SV VCF by splitting into a matrix of [sample, chromosome]"
    }
    input {
        Array[File] wholegenome_bams_from_all_samples
        Array[File] wholegenome_bais_from_all_samples

        File wholegenome_joint_vcf
        File wholegenome_joint_vcf_tbi
        File wholegenome_joint_sv
        File wholegenome_joint_sv_tbi

        File chrs_of_interest
        File ref_dict

        String cohort_name
        String gcs_out_root_dir
    }
    output {
        File sharded_bams_manifest = SaveBamManifest.gcs_path
        File sharded_bais_manifest = SaveBaiManifest.gcs_path
        # File sharded_smallvar_vcf_manifest
        # File sharded_strucvar_vcf_manifest
    }

    Array[String] chromosomes = read_lines(chrs_of_interest)
    String bam_out_dir = sub(gcs_out_root_dir, "/+$", "") + "/~{cohort_name}/sharded_bams"
    ############################################################################
    # handle BAMs
    scatter (pair in zip(wholegenome_bams_from_all_samples, wholegenome_bais_from_all_samples)) { # each sample
        call BU.InferSampleName { input: bam = pair.left }
        call ShardWholeGenome.Split as SplitBamByChr { input:
            ref_dict = ref_dict,
            bam = pair.left, bai = pair.right
        }
        Array[Pair[String, Pair[File, File]]] chr_bam_bai = SplitBamByChr.id_bam_bai_of_shards
    }

    # save data
    scatter (sample_chr_bam_bai in zip(InferSampleName.sample_name, chr_bam_bai)) { # each sample
        String sample_name = sample_chr_bam_bai.left
        Array[Pair[String, Pair[File, File]]] chrName_bam_bai = sample_chr_bam_bai.right
        scatter (pp in chrName_bam_bai) { # each chromosome
            String chr = pp.left
            Pair[File, File] bam_bai = pp.right
            call FF.FinalizeToFile as SaveChoppedSampleChrBam { input:
                file = bam_bai.left,
                outdir = bam_out_dir + "/~{sample_name}",
                name = "~{sample_name}.~{chr}.bam"
            }
            call FF.FinalizeToFile as SaveChoppedSampleChrBai { input:
                file = bam_bai.right,
                outdir = bam_out_dir + "/~{sample_name}",
                name = "~{sample_name}.~{chr}.bam.bai"
            }
        }
        Array[String] ordered_chromosomes = chr
        call WriteOneLine as BamManifestLine { input:
            index = sample_name,
            line = SaveChoppedSampleChrBam.gcs_path
        }
        call WriteOneLine as BaiManifestLine { input:
            index = sample_name,
            line = SaveChoppedSampleChrBai.gcs_path
        }
    }
    Array[Array[String]] dummy = ordered_chromosomes
    call WriteOneLine as WriteHeaderLine { input:
        index = "sample",
        line = dummy[0]
    }

    # make and save manifest
    call GU.ConcatenateFiles as AssembleBamManifest { input:
        af = flatten([[WriteHeaderLine.one_liner], BamManifestLine.one_liner]),
        out_name = "~{cohort_name}.sharded_bams_manifest.tsv"
    }
    call GU.ConcatenateFiles as AssembleBaiManifest { input:
        af = flatten([[WriteHeaderLine.one_liner], BaiManifestLine.one_liner]),
        out_name = "~{cohort_name}.sharded_bais_manifest.tsv"
    }

    call FF.FinalizeToFile as SaveBamManifest { input:
        file = AssembleBamManifest.merged,
        outdir = bam_out_dir
    }

    call FF.FinalizeToFile as SaveBaiManifest { input:
        file = AssembleBaiManifest.merged,
        outdir = bam_out_dir
    }
    ############################################################################
    # # split VCF by chromosome, then by sample

}

task WriteOneLine {
    input {
        String index
        Array[String] line
    }

    output {
        File one_liner = "one_line.txt"
    }

    Array[String] padded = flatten([[index], line])
    command <<<
    set -euxo pipefail
        echo -e ~{sep = "\t" padded} > one_line.txt
    >>>
    runtime {
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
