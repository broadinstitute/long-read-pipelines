version 1.0

import "tasks/CCSPepper.wdl"
import "tasks/Finalize.wdl" as FF

workflow PhaseMergedSNVVCFs {
    input {
        File snv_vcf
        File? snv_vcf_tbi

        File matching_ccs_bam
        File matching_ccs_bai

        String sample_name

        Int phasing_memory

        File ref_map_file

        String gcs_out_root_dir
    }

    String prefix = basename(basename(snv_vcf, ".gz"), ".vcf") + ".phased_with." + basename(matching_ccs_bam, ".bam")
    String final_dir = sub(gcs_out_root_dir, "/$", "") + "/finalized_variants/~{sample_name}/CCS/small"

    Map[String, String] ref_map = read_map(ref_map_file)

    call CCSPepper.MarginPhase {
        input:
            bam = matching_ccs_bam,
            bai = matching_ccs_bai,

            unphased_vcf = snv_vcf,
            unphased_vcf_tbi = snv_vcf_tbi,

            ref_fasta = ref_map['fasta'],
            ref_fasta_fai = ref_map['fai'],

            memory = phasing_memory,

            out_prefix = prefix + '.margin_phase'
    }
    call FF.FinalizeToFile as FinalizeMarginVCF { input: outdir = final_dir, file = MarginPhase.phasedVCF }
    call FF.FinalizeToFile as FinalizeMarginTbi { input: outdir = final_dir, file = MarginPhase.phasedtbi }
    call FF.FinalizeToFile as FinalizeMarginBED { input: outdir = final_dir, file = MarginPhase.phaseset_bed }

    call WhatsHap {
        input :
            bam = matching_ccs_bam,
            bai = matching_ccs_bai,

            unphased_vcf = snv_vcf,
            unphased_vcf_tbi = snv_vcf_tbi,

            ref_fasta = ref_map['fasta'],
            ref_fasta_fai = ref_map['fai'],

            memory = phasing_memory,

            out_prefix = prefix + '.whatshap'
    }

    call FF.FinalizeToFile as FinalizeWhatsHapVCF { input: outdir = final_dir, file = WhatsHap.phasedVCF }
    call FF.FinalizeToFile as FinalizeWhatsHapTbi { input: outdir = final_dir, file = WhatsHap.phasedTbi }
    call FF.FinalizeToFile as FinalizeWhatsHapTar { input: outdir = final_dir, file = WhatsHap.stats_targz }

    output {
        File margin_phased_merged_vcf = FinalizeMarginVCF.gcs_path
        File margin_phased_merged_vcf_tbi = FinalizeMarginTbi.gcs_path
        File margin_phased_merged_bed = FinalizeMarginBED.gcs_path
        File whatshap_phased_merged_vcf = FinalizeWhatsHapVCF.gcs_path
        File whatshap_phased_merged_vcf_tbi = FinalizeWhatsHapTbi.gcs_path
        File whatshap_phased_merged_vcf_stats = FinalizeWhatsHapTar.gcs_path
    }
}

task WhatsHap {
    input {
        File bam
        File bai

        File unphased_vcf
        File? unphased_vcf_tbi

        File ref_fasta
        File ref_fasta_fai

        Int memory

        String out_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = if bam_sz > 200 then 2*bam_sz else bam_sz + 200

    Int cores = 64

    command <<<
        set -euxo

        whatshap phase \
            --indels \
            -o ~{out_prefix}.vcf \
            --reference=~{ref_fasta} \
            ~{unphased_vcf} \
            ~{bam}
        ls -lh .

        whatshap stats \
            --tsv ~{out_prefix}.stats.tsv \
            --block-list ~{out_prefix}.stats.txt \
            --gtf ~{out_prefix}.stats.gtf \
            ~{out_prefix}.vcf
        ls -lh .

        bgzip -c ~{out_prefix}.vcf > ~{out_prefix}.vcf.gz
        tabix -p vcf ~{out_prefix}.vcf.gz
        tar -czvf \
            ~{out_prefix}.stats.tar.gz \
            ~{out_prefix}.stats.tsv \
            ~{out_prefix}.stats.txt \
            ~{out_prefix}.stats.gtf
    >>>

    output {
        File phasedVCF = "~{out_prefix}.vcf.gz"
        File phasedTbi = "~{out_prefix}.vcf.gz.tbi"
        File stats_targz = "~{out_prefix}.stats.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cores,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       20,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-whatshap:1.2.1"
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