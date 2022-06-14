# Mitochondria Variant Calling Documentation

## Sample: HG00514.2
12,167 reads after subsetting 
<br />
<br />
<br />
<br />


## Genome Assembly:


1) Hifiasm + unfiltered reads w/o genome size parameter:
    - Quast report:
        - Genome fraction (%): 100
        - Duplication ratio: 84.7
        - Largest alignment: 17,534
        - Total aligned length: 1,404,365
        - GC content: 44.49
        - number of contigs: 59

2) Hifiasm + unfiltered reads w/ genome size (ran on unfiltered reads): 
- workflow failed: <br />
    ```./hifiasm -o "hifi_test" -t 4 --hg-size 16k ~{reads.fq}```

    ```awk '/^S/{print ">"$2;print $3}' hifi_test.bp.p_ctg.gfa > hifi_test.fa```

    ``` samtools faidx hifi_test.fa ``` 
    ``` wc -l hifi_test.fa.fai ```

3) Hifiasm + filtered reads w/ genome size
    - output: 59 contigs 
    - mean length: 18k
        - Alighment using minimap2: 
           -  ```./minimap2 -aYL --MD -t 8 ../chrM_ref.fa ../test_1.p_ctg.fa | samtools sort -o "chrM_test.bam" ``` 
                - potential further steps: 
                    - polishing: collapse alignment into a single contig (CD-HIT) http://weizhong-lab.ucsd.edu/cd-hit/
                        1) CONSENT:
                            - free(): invalid pointer. software not maintained.
                        2) abPOA: 
                            - segmentation fault
                        3) CD-HIT:
                            - downloaded latest release and uploaded to VM
                            - ```awk '/^>/ {print; next; } { seqlen = length($0); print seqlen}' file.fa```
                            - with wrapped fasta:
                                - ```awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}'```
                            - 99% identity ng threshold)
                            - output: 9 clusters with 1 with similar length
                        4) circlator
                        ```sudo docker run --rm -it -v /home/ewan/data:/data sangerpathogens/circlator circlator all /data/test_1.p_ctg.fa /data/classified_reads.for_assembly.bam /data/output_directory```
                            
                        - OR trim soft-clipped region and align separately or throw away 
                        - filter out contigs that don't stretch the full length 16k
                        - abPOA: collapse contig (talk to John about the python version) 
                        - OR: run assembly twice

            - Contig selection and trimming:
             ```samtools view  chrM_test.bam | grep -F "h1tg000049" | less -S > select_contig.txt```
             - minimap2 align again: 
             ```./minimap2 -aYL -R'@RG\tID:blah\tSM:HG00514' --MD -t 8 ../chrM_ref.fa contig_trimmed.fa | samtools sort -o "second_alignment_chrM.bam"```

            - run paftools (paftools_mito.wdl):
                - convert PAF to VCF: <br />
                ```awk 'BEGIN {OFS="\t"} {print "chrM", $1-1, $1, $2"->"$3}' paf_fin_out.txt > mito_paf.bed```

    - align reads to assembly using Clair3: <br />
        - align filtered reads to contig (check if contig deviates from reads) <br />
        ```./minimap2 -aYL -R'@RG\tID:ID\tSM:HG00514' --MD -t 8 contig_trimmed.fa chrM_raw_reads.fa | samtools sort -o "readstocontig_chrM.bam"```
    - Clair3: evaluate assembly - compare raw reads to selected contig
        - output: pileup.vcf; full_alignment.vcf; merged.vcf
        - visualize vcf output on IGV with readstocontig_chrM.bam (change reference to contig_trimmed.fa on IGV)
        - IGV shows:
            - variants identified (shown in vcf) are minor (not seen in most reads aka low allele frequency)
            - if selected contig is the representative one, there should be almost no variants?? 
        


 
                        
<br />
<br />


3) Canu 2.2 with filtered reads:
    -  canu-correct issues:
        1) error due to less than 100 reads being corrected
        2) reads without overlapping

        ```
        --                             original      original
        --                            raw reads     raw reads
        --   category                w/overlaps  w/o/overlaps
        --   -------------------- ------------- -------------
        --   Number of Reads                282          7189
        --   Number of Bases            3066432        136616
        --   Coverage                   191.652         8.539
        --   Median                       10382             0
        --   Mean                         10873            19
        --   N50                          10577          8563
        --   Minimum                       8277             0
        --   Maximum                      16615          9539
        ```


    - Two approaches for correcting error caused by < 100 reads: 
    1) change source code 
    2) use parameter ```corOutCoverage```

    -  ```canu -correct corOutCoverage=100```
        - Quast stats:
            - 1 contig
            - largest alignment: 15,689
            - duplication ratio: 1.887
            - total length: 31,266
            - GC %: 44.42



4) MitoHiFi
    - container setup:
        ```
        sudo docker login
        sudo docker build .
        sudo docker images
        sudo docker tag <image name> <dockerID/imageID>
        sudo docker push <dockerID/imageID>
        sudo docker pull <dockerID/imageID>
        sudo docker run -it <imageID>
        ```
        convert fastq to fa:
        ```
        seqtk/seqtk seq -a classified_reads.for_assembly.fastq.gz > chrM_reads.fa
        ```
    - mount necessary files: 
        ```
        sudo docker run -it -v /home/ewan/MitoHiFi/reads_ref:/powerhouse_reads --entrypoint /bin/bash f84b6087a290
        ```
    - get GenBank file:
        ```
        findMitoReference.py --species "Homo sapiens" --email ewan@broadinstitute.org --outfolder /bin --min_length 16000
        ```

    - assemble using mitohifi.py

        ```
         mitohifi.py -r ../powerhouse_reads/chrM_reads.fa -f ../powerhouse_reads/chrM_ref.fa -g ../powerhouse_reads/LC657585.1.gb -t 10
        ```

        1. First we map your Pacbio HiFi reads to the close-related mitogenome
        2. Now we filter out any mapped reads that are larger than the reference mitogenome to avoid NUMTS
            2.1 First we convert the mapped reads from BAM to FASTA format: 
            ```samtools fasta reads.HiFiMapped.bam > gbk.HiFiMapped.bam.fasta```
            Total number of mapped reads: 7471
            2.2 Then we filter reads that are larger than 16569 bp. Number of filtered reads: 7286
        3. Now let's run hifiasm to assemble the mapped and filtered reads!
            ```hifiasm --primary -t 10 -f 0 -o gbk.HiFiMapped.bam.filtered.assembled gbk.HiFiMapped.bam.filtered.fasta```



### CD-HIT: 
- clustering and comparing protein or sequences 
- compares 2 datasets and identifies the sequences that are similar 
- identifies duplicates from single or paired Illumina reads
- identifies overlapping reads

BED file and visualize 

future work: 
- clustering variants (should see cluster around ethnicity)




create truth dataset:
- supporting reads





Hifiasm notes: 
- prefix`.r_utg.gfa: haplotype-resolved raw unitig graph. This graph keeps all haplotype information.
- `prefix`.p_utg.gfa: haplotype-resolved processed unitig graph without small bubbles. Small bubbles might
be caused by somatic mutations or noise in data, which are not the real haplotype information. Hifiasm automatically pops such small bubbles based on coverage
- `prefix`.p_ctg.gfa: assembly graph of primary contigs. This graph includes a complete assembly with long
stretches of phased blocks.
- `prefix`.a_ctg.gfa: assembly graph of alternate contigs. This graph consists of all contigs that are discarded
in primary contig graph.
- `prefix`.*hap*.p_ctg.gfa: phased contig graph. This graph keeps the phased contigs.



hifi parameters: 
- N 20 didn't change contig numbers but looks diff on IGV


./hifiasm -o "523_test" -t 4 -k 51 --pri-range 100 16700 --n-hap 1 --hg-size 16k chrM_filtered_rawreads.fastq