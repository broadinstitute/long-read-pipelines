version 1.0

import "../../../tasks/Utility/Finalize.wdl" as FF
import "PfalciparumTypeDrugResistanceMarkers.wdl" as DRUG_RES

workflow PfalciparumDrugResistanceSummary {
    meta {
        desciption: "Create a drug resistance report based on the given raw drug resistance loci report."
    }

    input {
        File raw_drug_resistance_info

        String participant_name

        String? gcs_out_root_dir
    }

    parameter_meta {
        raw_drug_resistance_report: "File containing a raw drug resistance report to use to determine drug resistance."
        participant_name:    "Participant (or sample) name for the given bam file."
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs.  If no bucket is given, outputs will not be finalized and instead will remain in their native execution location."
    }

    call DRUG_RES.CreateDrugResistanceSummary as CreateDrugResistanceSummary {
        input:
            raw_drug_resistance_info = raw_drug_resistance_info,
            prefix = participant_name
    }

    if (defined(gcs_out_root_dir)) {
        String concrete_gcs_out_root_dir = select_first([gcs_out_root_dir])
        String outdir = sub(concrete_gcs_out_root_dir, "/$", "") + "/PfalciparumDrugResistanceSummary/~{participant_name}"
        call FF.FinalizeToFile as FinalizeDrugResistanceSummary { input: outdir = outdir, file = CreateDrugResistanceSummary.resistance_summary }
    }

    output {
        File drug_resistance_summary = select_first([FinalizeDrugResistanceSummary.gcs_path, CreateDrugResistanceSummary.resistance_summary])

        String drug_status_chloroquine   = CreateDrugResistanceSummary.chloroquine_status
        String drug_status_pyrimethamine = CreateDrugResistanceSummary.pyrimethamine_status
        String drug_status_sulfadoxine   = CreateDrugResistanceSummary.sulfadoxine_status
        String drug_status_mefloquine    = CreateDrugResistanceSummary.mefloquine_status
        String drug_status_artemisinin   = CreateDrugResistanceSummary.artemisinin_status
        String drug_status_piperaquine   = CreateDrugResistanceSummary.piperaquine_status
    }
}
