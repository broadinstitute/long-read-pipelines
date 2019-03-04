import argparse
from google.cloud import bigquery
import os
import peewee as pw
import tqdm

from vcf_utils import parse_vcf_file
from bigquery_utils import export_to_bigquery


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--description", help="(BigQuery only) description of table contents. This will be saved in table metadata")
    p.add_argument("--dataset-name", help="(BigQuery only) dataset name", choices=["weisburd_grch37", "weisburd_grch38"])
    p.add_argument("--project-id", help="(BigQuery only) gcloud project id", default="broad-dsp-spec-ops")
    p.add_argument("--limit", type=int, help="Number of records to read from the .dat file. Useful for debugging.")

    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--sqlite", action="store_true", help="export to sqlite")
    g.add_argument("-bq", "--bigquery", action="store_true", help="export to BigQuery")

    p.add_argument("-t", "--table-name", help="Table name in BigQuery or filename prefix in sqlite. Defaults to the vcf filename")
    p.add_argument("vcf_path")

    args = p.parse_args()

    if ".vcf" not in args.vcf_path:
        p.error(f"Input file {args.vcf_path} must have a .vcf* extension")

    return args


VCF_SCHEMA = [
    # type: STRING, BYTES, INTEGER, FLOAT, BOOLEAN, BOOL TIMESTAMP, DATE, TIME, DATETIME, RECORD or STRUCT (same as RECORD) (indicates that the field contains a nested schema).
    # mode: NULLABLE, REQUIRED and REPEATED. The default value is NULLABLE.
    bigquery.SchemaField("chrom", "STRING", mode="REQUIRED"),
    bigquery.SchemaField("start", "INTEGER", mode="REQUIRED", description="0-based, like in .bed format"),
    bigquery.SchemaField("end", "INTEGER", mode="REQUIRED", description="1-based, like in .bed format"),
    bigquery.SchemaField("ref", "STRING", mode="REQUIRED"),
    bigquery.SchemaField("alt", "STRING", mode="REQUIRED"),

    bigquery.SchemaField("repeat_unit", "STRING", mode="REQUIRED"),
    bigquery.SchemaField("repeat_count", "FLOAT", mode="REQUIRED"),
    bigquery.SchemaField("repeat_unit_length", "INTEGER", mode="REQUIRED"),

    bigquery.SchemaField("variant_length", "INTEGER", mode="REQUIRED"),
    bigquery.SchemaField("variant_type", "STRING", mode="REQUIRED"),
    bigquery.SchemaField("repeat_class", "STRING", mode="REQUIRED"),
    bigquery.SchemaField("repeat_count_in_variant", "FLOAT", mode="REQUIRED"),
    bigquery.SchemaField("repeat_count_in_variant_and_ref_left", "FLOAT", mode="REQUIRED"),
    bigquery.SchemaField("repeat_count_in_variant_and_ref_right", "FLOAT", mode="REQUIRED"),

    bigquery.SchemaField("full_repeat_sequence", "STRING", mode="NULLABLE"),
    bigquery.SchemaField("full_repeat_sequence_length", "INTEGER", mode="NULLABLE"),

    bigquery.SchemaField("lacks_matching_repeats_in_ref", "BOOL", mode="REQUIRED"),
    bigquery.SchemaField("left_entropy", "FLOAT", mode="NULLABLE"),
    bigquery.SchemaField("right_entropy", "FLOAT", mode="NULLABLE"),
    bigquery.SchemaField("entropy", "FLOAT", mode="NULLABLE"),
    #bigquery.SchemaField("full_repeat_sequence", "STRING", mode="NULLABLE"),
    #bigquery.SchemaField("full_repeat_sequence_length", "INTEGER", mode="NULLABLE"),
    bigquery.SchemaField("left_flank", "STRING", mode="NULLABLE"),
    bigquery.SchemaField("right_flank", "STRING", mode="NULLABLE"),
]


def generate_db_record(vcf_record):
    r = vcf_record

    repeat_count = max(
        float(r["RepeatCountInVariant"]),
        float(r["RepeatCountInVariantAndRefLeft"]),
        float(r["RepeatCountInVariantAndRefRight"]),
        float(r["RepeatCountInVariantAndRefLeft"]) + float(r["RepeatCountInVariantAndRefRight"]) - float(r["RepeatCountInVariant"]))

    repeat_unit_length = len(r["RU"])
    full_repeat_sequence = r["RU"]*int(repeat_count)
    if int(repeat_count) < repeat_count:
        full_repeat_sequence += r["RU"][:int(repeat_unit_length * repeat_count) - repeat_unit_length * int(repeat_count)]

    assert "," not in r["alt"]

    return {
        "chrom": r["chrom"],
        "start": int(r["pos"]) - 1,
        "end": int(r["pos"]) + len(r["ref"]) - 1,
        "ref": r["ref"],
        "alt": r["alt"],

        "repeat_unit": r["RU"],
        "repeat_count": repeat_count,

        "variant_length": max(len(r["ref"]), len(r["alt"])) - 1,
        "variant_type": "DEL" if len(r["ref"]) > len(r["alt"]) else ("INS" if len(r["ref"]) < len(r["alt"]) else "SNP"),
        "repeat_unit_length": len(r["RU"]),

        "repeat_class": r["RepeatClass"],
        "repeat_count_in_variant": float(r["RepeatCountInVariant"]),

        "repeat_count_in_variant_and_ref_left": float(r["RepeatCountInVariantAndRefLeft"]),
        "repeat_count_in_variant_and_ref_right": float(r["RepeatCountInVariantAndRefRight"]),

        "full_repeat_sequence": full_repeat_sequence,
        "full_repeat_sequence_length": len(full_repeat_sequence),

        "lacks_matching_repeats_in_ref": r["LacksMatchingRepeatsInRef"].lower() == "true",
        "left_entropy": float(r["LeftEntropy"]),
        "right_entropy": float(r["RightEntropy"]),
        "entropy": float(r["RightEntropy"] if float(r["RepeatCountInVariantAndRefRight"]) > 0 else r["LeftEntropy"]),

        "left_flank": r["LeftFlank"],
        "right_flank": r["RightFlank"],
    }


def main():

    args = parse_args()

    if not args.table_name:
        args.table_name = os.path.basename(args.vcf_path).replace(".vcf", "").replace(".gz", "")

    print(f"Parsing and loading table {args.table_name} \n      from {args.vcf_path}")
    records = tqdm.tqdm(parse_vcf_file(args.vcf_path, limit=args.limit), unit=" records")

    if args.bigquery:
        if args.dataset_name is None:
            raise ValueError("--dataset-name must be specified")

        client = bigquery.Client(project=args.project_id)
        table_ref = client.dataset(args.dataset_name).table(args.table_name)
        table = bigquery.Table(table_ref, schema=VCF_SCHEMA)

        if any([t._properties['tableReference']['tableId'] == args.table_name for t in client.list_tables(args.dataset_name)]):
            print(f"Deleting table: {table}")
            client.delete_table(table)

        print(f"Creating table: {table}")
        client.create_table(table)

        export_to_bigquery((generate_db_record(r) for r in records),
            schema=VCF_SCHEMA,
            table_name=args.table_name,
            dataset_name=args.dataset_name,
            project_id=args.project_id,
            table_description=f"source: {os.path.abspath(os.path.expanduser(args.vcf_path))}\n{args.description or ''}",
            delete_table_if_exists=True,
            verbose=True)

    if args.sqlite:
        print(f"Creating sqlite database {args.table_name}.db")
        db = pw.SqliteDatabase(f"{args.table_name}.db")

        class VCFRecord(pw.Model):
            chrom = pw.CharField(max_length=5)
            start = pw.IntegerField(index=True)
            end = pw.IntegerField(index=True)
            ref = pw.TextField()
            alt = pw.TextField()

            repeat_unit = pw.TextField(index=True)
            repeat_count = pw.FloatField(index=True)
            repeat_unit_length = pw.FloatField(index=True)
            full_repeat_sequence = pw.TextField(null=True)
            full_repeat_sequence_length = pw.FloatField(null=True)

            #percent_matches = pw.FloatField(null=True)
            #percent_indels = pw.FloatField(null=True)
            #alignment_score = pw.FloatField(null=True)
            entropy = pw.FloatField(null=True)
            #full_repeat_sequence = pw.TextField(null=True)
            #full_repeat_sequence_length = pw.FloatField(null=True)
            left_flank = pw.TextField(null=True)
            right_flank = pw.TextField(null=True)
            repeat_class = pw.TextField(null=True)
            repeat_count_in_variant = pw.FloatField(index=True, null=True)
            repeat_count_in_variant_and_ref_left = pw.FloatField(index=True, null=True)
            repeat_count_in_variant_and_ref_right = pw.FloatField(index=True, null=True)
            lacks_matching_repeats_in_ref = pw.BooleanField(index=True)
            left_entropy = pw.FloatField(null=True)
            right_entropy = pw.FloatField(null=True)

            class Meta:
                database = db
                indexes = (
                    (('chrom', 'start', 'end', 'repeat_unit'), False),
            )

        Record = VCFRecord

        if Record.table_exists():
            Record.delete().execute()

        Record.create_table()
        with db.atomic():
            Record.bulk_create((
                Record(**generate_db_record(r)) for r in records), batch_size=1000)

if __name__ == "__main__":
    main()
