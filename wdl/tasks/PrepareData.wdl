version 1.0

import "Structs.wdl"
import "Utils.wdl" as Utils

workflow PrepareData {
    input {
        String gcs_dir
        String sample_name
        Int? num_reads_per_split
    }

    call Utils.DetectRunInfo as DetectRunInfo {
        input:
            gcs_dir = gcs_dir,
            sample_name = sample_name,
    }
    String platform = DetectRunInfo.run_info['PL']
    String id = DetectRunInfo.run_info['ID']

    call Utils.ShardLongReads as ShardLongReads { input: unmapped_files = DetectRunInfo.files, num_reads_per_split = num_reads_per_split }
    call Utils.PrepareManifest as PrepareFileManifest { input: files = DetectRunInfo.files }
    call Utils.PrepareManifest as PrepareShardManifest { input: files = ShardLongReads.unmapped_shards }

    File manifest = if (DetectRunInfo.run_info['TY'] == "BAM" && length(DetectRunInfo.files) == 1) then PrepareShardManifest.manifest else PrepareFileManifest.manifest
    Int manifest_lines_per_chunk = if (DetectRunInfo.run_info['TY'] == "BAM" && length(DetectRunInfo.files) == 1) then 1 else 100

    call Utils.ChunkManifest as ChunkManifest { input: manifest = manifest, manifest_lines_per_chunk = manifest_lines_per_chunk }

    output {
        Array[File] manifest_chunks = ChunkManifest.manifest_chunks
        Map[String, String] run_info = DetectRunInfo.run_info
        Array[String] files = DetectRunInfo.files
    }
}
