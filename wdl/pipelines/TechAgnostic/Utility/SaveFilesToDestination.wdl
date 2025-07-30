version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/Finalize.wdl" as FF

struct FinalizationManifestLine {
    Array[File]+ files_to_save
    Array[String]? file_names       # must be same length as files_to_save, if provided; ignored if pack_name is also provided
    String? pack_name               # if provided, file_names will be ignored: single file will be (block-)gzipped and multiple files will be tar-gzed.
    Boolean is_singleton_file       # if this is just a single file [when files_to_save is intended to be an array, even though it might be length 1, this should be false]
    String destination              # where to save the files;
                                    #  if saving a single file, or if packing multiple files, then the result will be saved under this destination
                                    #  if saving multiple files but not packing them, then the files will be saved under this destination with their default or custom names
    Boolean? individual_compress    # has effect only for multiple files; if true, will
    Boolean? use_bgzip              # block-gzip files or not
    String output_attribute_name    # mostly for terra usage
}

workflow SaveFilestoDestination {
    meta {
        desciption:
        "If your workflow needs to save outputs from various tasks into slightly different locations, this is the sub-workflow for you."
    }
    parameter_meta {
        instructions:
        "A workflow, when calling this sub-workflow, should construct this array by hand"
        result:
        "The gs path for each line in your instructions"
        key_file:
        "Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
    }
    input {
        Array[FinalizationManifestLine] instructions
        Array[Map[String, String]]? already_finalized
        File? key_file
    }
    output {
        Map[String, String] result = select_first([PackAllSavings.merged, collect.output_map])
    }

    scatter (line in instructions) {

        if (line.is_singleton_file) {
            if (!defined(line.pack_name)) {  # as-is
                call FF.FinalizeToFile as SF { input:
                    file = line.files_to_save[0],
                    name = if (defined(line.file_names)) then select_first([line.file_names])[0] else basename(line.files_to_save[0]),
                    outdir = line.destination
                }
            }
            if (defined(line.pack_name)) {  # (block-)gzip single file
                call FF.CompressAndFinalize as SF_gz { input:
                    file = line.files_to_save[0],
                    name = select_first([line.pack_name]),
                    block_gzip = if (defined(line.use_bgzip)) then select_first([line.use_bgzip]) else false,
                    outdir = line.destination
                }
            }
        }

        if (!line.is_singleton_file) {
            if (!defined(line.pack_name)) {  # as-is
                call FF.FinalizeToDir as MF_default { input:
                    files = line.files_to_save,
                    file_names = line.file_names,
                    outdir = line.destination
                }
            }

            if (defined(line.individual_compress)) {  # (b)gz individual files
                call FF.FinalizeAndCompress as MF_gz { input:
                    files = line.files_to_save,
                    outdir = line.destination,
                    block_gzip = if (defined(line.use_bgzip)) then select_first([line.use_bgzip]) else false
                }
            }

            if (defined(line.pack_name)) { # tar.gz to a single file
                call FF.TarGZFilesAndSave as MF_pack { input:
                    files = line.files_to_save,
                    name = select_first([line.pack_name]),
                    outdir = line.destination
                }
            }
        }

        String attr = line.output_attribute_name
        String final_path = select_first([SF_gz.gcs_path, SF.gcs_path,
                                          MF_default.gcs_dir, MF_gz.gcs_path, MF_pack.gcs_path])
    }

    call GU.CoerceArrayOfPairsToMap as collect { input: keys = attr, values = final_path }

    Boolean concatenate = (defined(already_finalized)) && (length(select_first([already_finalized])) > 0)
    if (concatenate) {
        scatter(mm in select_first([already_finalized])) {
            call GU.MapToTsv { input: m = mm}
        }
        call GU.ConcatenateFiles { input: af = MapToTsv.tsv, out_name = "does_not_matter.tsv" }

        Map[String, String] pack_one = read_map(ConcatenateFiles.merged)

        call GU.MergeMaps as PackAllSavings { input:
            one = pack_one, two = collect.output_map
        }
    }
}
