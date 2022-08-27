version 1.0

import "tasks/Fusilli.wdl" as Fusilli

workflow FusilliCall {
    input {
        String fusilli_db_gcs

        File ref_meta
        Array[String] ref_ids
        Array[File] ref_fastas
        Array[File] ref_gffs

        String fusilli_run_gcs
        String fusilli_call_id
        Array[String] sample_ids
        File? hmm_config
    }

    # Perform Tesserae2 alignment across multiple shards (by default 1000 per worker)
    scatter(sample_id in sample_ids) {
        String sample_contigs = fusilli_run_gcs + "/" + sample_id + "/" + sample_id + ".contigs.fasta"

        call Fusilli.ChunkSampleContigs as ChunkSampleContigs {
            input:
                sample_contigs = sample_contigs
        }

        scatter(chunk in ChunkSampleContigs.chunks) {
            if(defined(hmm_config) && hmm_config != "") {
                call Fusilli.TesseraeAlign as TesseraeAlignWithConfig {
                    input:
                        sample_id = sample_id,
                        sample_contigs = chunk,
                        fusilli_run_gcs = fusilli_run_gcs,
                        hmm_config = hmm_config
                }
            }

            if(!defined(hmm_config) || hmm_config == "") {
               call Fusilli.TesseraeAlign as TesseraeAlign {
                   input:
                       sample_id = sample_id,
                       sample_contigs = chunk,
                       fusilli_run_gcs = fusilli_run_gcs
               }
            }

            Array[File] alignments = select_first([TesseraeAlignWithConfig.alignments, TesseraeAlign.alignments])
            File unaligned_fasta = select_first([TesseraeAlignWithConfig.unaligned, TesseraeAlign.unaligned])
        }

        call Fusilli.MakeCalls as MakeCalls {
            input:
                fusilli_db_gcs = fusilli_db_gcs,
                ref_ids = ref_ids,
                ref_fastas = ref_fastas,

                fusilli_run_gcs = fusilli_run_gcs,
                sample_id = sample_id,
                aligned_sample_contigs = flatten(alignments)
        }

        call Fusilli.FinalizeCalls as FinalizeCalls {
            input:
                fusilli_run_gcs = fusilli_run_gcs,
                call_id = fusilli_call_id,
                sample_id = sample_id,

                alignment_stats = MakeCalls.alignment_stats,
                unaligned_fastas = unaligned_fasta,
                aligned_sample_contigs = flatten(alignments),
                aligned_ref_contigs = MakeCalls.aligned_ref_contigs,
                aligned_ref_contigs_bai = MakeCalls.aligned_ref_contigs_bai,
                per_ref_vcfs = MakeCalls.per_ref_vcfs,
                per_ref_bed = MakeCalls.per_ref_bed
        }
    }

    call Fusilli.CreatePseudoRef as CreatePseudoRef {
        input:
            ref_meta = ref_meta,
            ref_ids = ref_ids,
            ref_fastas = ref_fastas,
            ref_gffs = ref_gffs,

            alignment_stats = MakeCalls.alignment_stats
    }

    scatter(i in range(length(sample_ids))) {
        call Fusilli.TranslateCoordinatesToPseudoRef as TranslateCoordinatesToPseudoRef {
            input:
                pseudo_ref = CreatePseudoRef.pseudo_ref,
                sample_id = sample_ids[i],
                per_ref_vcfs = MakeCalls.per_ref_vcfs[i],
                per_ref_bed = MakeCalls.per_ref_bed[i]
        }
    }

    call Fusilli.FinalizePseudoRef as FinalizePseudoRef {
        input:
            fusilli_run_gcs = fusilli_run_gcs,
            call_id = fusilli_call_id,

            pseudo_ref = CreatePseudoRef.pseudo_ref,
            pseudo_ref_fai = CreatePseudoRef.pseudo_ref_fai,
            pseudo_gff = CreatePseudoRef.pseudo_gff,
            pseudo_vcfs = TranslateCoordinatesToPseudoRef.pseudo_vcf,
            pseudo_vcf_csis = TranslateCoordinatesToPseudoRef.pseudo_vcf_csi,
            pseudo_beds = TranslateCoordinatesToPseudoRef.pseudo_bed,
    }

    output {
        String pseudo_ref = FinalizePseudoRef.pseudo_ref
        String pseudo_gff = FinalizePseudoRef.pseudo_gff
        String pseudo_sample_data_dir = FinalizePseudoRef.sample_data_dir
    }
}
