version 1.0

import "ONTMethylation.wdl" as M
import "tasks/Finalize.wdl" as FF

workflow ONTMergeMethylation {
  input {
    Array[File] mod_dbs
    Array[File] mod_txts
    String prefix

    String gcs_out_root_dir
  }

  String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTMethylation/~{prefix}"

  if (length(mod_dbs) > 1) {
    call M.MergeModifiedBaseCallDBs as MergeDBs { input: dbs = mod_dbs }
  }

  File merged_txt = select_first([MergeDBs.per_read_modified_base_calls_txt, mod_txts[0]])

  call FF.FinalizeToFile as FinalizeModifiedTxt { input: file=merged_txt, outdir=outdir, name="~{prefix}.per_read_modified_base_calls.txt" }

  output {
    File per_read_modified_base_calls_txt = FinalizeModifiedTxt.gcs_path
  }

}
