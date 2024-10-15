version 1.0

workflow BedtoolsIntersection {
    input {
        String sample_id
        File eval_vcf
        File eval_vcf_index
        File bed_file
    }
    call Intersection { input:
        joint_vcf = eval_vcf,
        joint_vcf_index = eval_vcf_index,
        bedfile = bed_file,
        samplename = sample_id
    }


     output {
        File result = Intersection.outputfile
        
    }
}


task Intersection {
    input{       
        File joint_vcf
        File joint_vcf_index
        File bedfile
        String samplename
    }
    
    command <<<
        set -x pipefail

        time bcftools view --threads 4 -s ~{samplename} ~{joint_vcf} -o ~{samplename}.subset.g.vcf.gz

        bcftools view -e 'GT="0/0"' ~{samplename}.subset.g.vcf.gz | bcftools query -f'%CHROM\t%POS0\t%POS\n' > filtered_variants.txt

        bedtools intersect -c -a ~{bedfile} -b filtered_variants.txt > variant_counts.txt

        awk 'FNR==NR {count[$1"\t"$2"\t"$3]=$4; next} {print $0"\t"count[$1"\t"$2"\t"$3]}' variant_counts.txt "$BED_FILE" > bed_with_counts.bed




    >>>
    
    output {
		File outputfile = "bed_with_counts.bed"
    }


    Int disk_size = 1 + ceil(2 * (size(joint_vcf, "GiB")))

    runtime {
        cpu: 4
        memory: "24 GiB"
        disks: "local-disk " + disk_size + " SSD" #"local-disk 100 HDD"
        preemptible: 1
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}
