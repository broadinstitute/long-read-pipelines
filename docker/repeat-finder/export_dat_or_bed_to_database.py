import argparse
from collections import defaultdict
from google.cloud import bigquery
import os
import peewee as pw
import tqdm
from bed_utils import parse_bed_file
from dat_utils import parse_dat_file
from bigquery_utils import export_to_bigquery


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--dat-in-original-format", action="store_true", help="For .dat files generated without the trf -ngs flag")
    p.add_argument("--description", help="(BigQuery only) description of table contents. This will be saved in table metadata")
    p.add_argument("--dataset-name", help="(BigQuery only) dataset name", choices=["weisburd_grch37", "weisburd_grch38"])
    p.add_argument("--project-id", help="(BigQuery only) gcloud project id", default="broad-dsp-spec-ops")
    p.add_argument("--limit", type=int, help="Number of records to read from the .dat file. Useful for debugging.")
    #p.add_argument("-m", "--merge", action="store_true", help="merge with existing intervals in database")

    p.add_argument("--stretch-tsv", action="store_true", help=".bed file is from STRetch output")

    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--sqlite", action="store_true", help="export to sqlite")
    g.add_argument("-bq", "--bigquery", action="store_true", help="export to BigQuery")

    p.add_argument("-t", "--table-name", help="Table name in BigQuery or filename prefix in sqlite. Defaults to the dat_or_bed_file_path filename")
    p.add_argument("dat_or_bed_file_path")

    args = p.parse_args()

    if args.dat_or_bed_file_path.endswith(".dat"):
        args.is_bed_file = False
        args.is_dat_file = True
        args.is_stretch_tsv = False
    elif args.dat_or_bed_file_path.endswith('.bed'):
        args.is_bed_file = True
        args.is_dat_file = False
        args.is_stretch_tsv = False
    elif args.dat_or_bed_file_path.endswith('.tsv') and args.stretch_tsv:
        args.is_bed_file = False
        args.is_dat_file = False
        args.is_stretch_tsv = True
    else:
        p.error(f"Input file {args.dat_or_bed_file_path} must have .dat or .bed extension or be a .tsv file and --stretch-tsv must be specified")

    return args

BED_SCHEMA = [
    # type: STRING, BYTES, INTEGER, FLOAT, BOOLEAN, BOOL TIMESTAMP, DATE, TIME, DATETIME, RECORD or STRUCT (same as RECORD) (indicates that the field contains a nested schema).
    # mode: NULLABLE, REQUIRED and REPEATED. The default value is NULLABLE.
    bigquery.SchemaField("chrom", "STRING", mode="REQUIRED"),
    bigquery.SchemaField("start", "INTEGER", mode="REQUIRED", description="0-based, like in .bed format"),
    bigquery.SchemaField("end", "INTEGER", mode="REQUIRED", description="1-based, like in .bed format"),
    bigquery.SchemaField("repeat_unit", "STRING", mode="REQUIRED"),
    bigquery.SchemaField("repeat_count", "FLOAT", mode="REQUIRED"),
    bigquery.SchemaField("repeat_unit_length", "INTEGER", mode="REQUIRED"),
    bigquery.SchemaField("full_repeat_sequence", "STRING", mode="NULLABLE"),
    bigquery.SchemaField("full_repeat_sequence_length", "INTEGER", mode="NULLABLE"),
]

STRETCH_SCHEMA = list(BED_SCHEMA) + [
    bigquery.SchemaField("p_adj", "FLOAT", mode="REQUIRED"),
    bigquery.SchemaField("z_score", "FLOAT", mode="REQUIRED"),
    bigquery.SchemaField("coverage", "FLOAT", mode="REQUIRED"),

    bigquery.SchemaField("variant_type", "STRING", mode="REQUIRED"),
    bigquery.SchemaField("variant_length", "INTEGER", mode="REQUIRED"),
    bigquery.SchemaField("repeat_count_in_variant", "FLOAT", mode="REQUIRED"),
    bigquery.SchemaField("lacks_matching_repeats_in_ref", "BOOL", mode="REQUIRED"),
]

DAT_SCHEMA = list(BED_SCHEMA) + [
    bigquery.SchemaField("percent_matches", "FLOAT", mode="NULLABLE"),
    bigquery.SchemaField("percent_indels", "FLOAT", mode="NULLABLE"),
    bigquery.SchemaField("alignment_score", "FLOAT", mode="NULLABLE"),
    bigquery.SchemaField("entropy", "FLOAT", mode="NULLABLE"),
    #bigquery.SchemaField("full_repeat_sequence", "STRING", mode="NULLABLE"),
    #bigquery.SchemaField("full_repeat_sequence_length", "INTEGER", mode="NULLABLE"),
    bigquery.SchemaField("left_flank", "STRING", mode="NULLABLE"),
    bigquery.SchemaField("right_flank", "STRING", mode="NULLABLE"),
]


def generate_db_record(dat_or_bed_record):
    r = dat_or_bed_record
    if "full_repeat_sequence_length" not in r:
        r["full_repeat_sequence_length"] = r["end"] - r["start"]

    if "repeat_unit_length" not in r:
        r["repeat_unit_length"] = len(r["repeat_unit"])

    return r


def main():

    args = parse_args()

    if not args.table_name:
        args.table_name = os.path.splitext(os.path.basename(args.dat_or_bed_file_path))[0]

    print(f"Parsing and loading table {args.table_name} \n      from {args.dat_or_bed_file_path}")
    if args.is_dat_file:
        records = tqdm.tqdm(parse_dat_file(args.dat_or_bed_file_path, dat_in_original_format=args.dat_in_original_format, limit=args.limit), unit=" records")
        schema = DAT_SCHEMA
    elif args.is_bed_file:
        records = tqdm.tqdm(parse_bed_file(args.dat_or_bed_file_path, limit=args.limit), unit=" records")
        schema = BED_SCHEMA
    elif args.is_stretch_tsv:
        schema = STRETCH_SCHEMA
        records = []
        with open(args.dat_or_bed_file_path) as tsv_file:
            for i, line in enumerate(tsv_file):
                line = line.rstrip("\n")
                fields = line.split("\t")
                if i == 0:
                    header = fields
                    assert header[0] == "chrom" # header should be: chrom   start   end     sample  repeatunit      reflen  locuscoverage   outlier p_adj   bpInsertion     repeatUnit
                    continue

                fields_dict = dict(zip(header, fields))

                record = {
                    "chrom": fields_dict["chrom"],
                    "start": int(fields_dict["start"]) - 1,  # convert to 0-based coordinates in .bed
                    "end": int(fields_dict["end"]),
                    "repeat_unit":  fields_dict["repeatunit"],
                    "repeat_count": float(fields_dict["repeatUnits"]),
                    "repeat_count_in_variant": float(fields_dict["repeatUnits"]) - float(fields_dict["reflen"]),  # how many repeat units are in the variant = total repeat units - repeat units in reference
                    "repeat_unit_length": len(fields_dict["repeatunit"]),
                    "p_adj": float(fields_dict["p_adj"]),
                    "z_score": float(fields_dict["outlier"]),
                    "coverage": float(fields_dict["locuscoverage"]),
                    "variant_type": "INS",
                    "lacks_matching_repeats_in_ref": False,
                }

                record["variant_length"] = int(record["repeat_unit_length"] * record["repeat_count_in_variant"])
                record["full_repeat_sequence_length"] = int(record["repeat_unit_length"] * record["repeat_count"])
                records.append(record)


    else:
        raise ValueError("Internal error")


    if args.bigquery:
        if args.dataset_name is None:
            raise ValueError("--dataset-name must be specified")

        #if args.merge:
        #    raise ValueError("--merge not supported with --bigquery")

        client = bigquery.Client(project=args.project_id)
        table_ref = client.dataset(args.dataset_name).table(args.table_name)
        table = bigquery.Table(table_ref, schema=schema)

        if any([t._properties['tableReference']['tableId'] == args.table_name for t in client.list_tables(args.dataset_name)]):
            print(f"Deleting table: {table}")
            client.delete_table(table)

        print(f"Creating table: {table}")
        client.create_table(table)

        export_to_bigquery((generate_db_record(r if isinstance(r, dict) else r._asdict()) for r in records),
            schema=schema,
            table_name=args.table_name,
            dataset_name=args.dataset_name,
            project_id=args.project_id,
            table_description=f"source: {os.path.abspath(os.path.expanduser(args.dat_or_bed_file_path))}\n{args.description or ''}",
            delete_table_if_exists=True,
            verbose=True)

    if args.sqlite:
        print(f"Creating sqlite database {args.table_name}.db")
        db = pw.SqliteDatabase(f"{args.table_name}.db")

        class BedRecord(pw.Model):
            chrom = pw.CharField(max_length=5)
            start = pw.IntegerField(index=True)
            end = pw.IntegerField(index=True)
            repeat_unit = pw.TextField(index=True)
            repeat_count = pw.FloatField(index=True)
            repeat_unit_length = pw.FloatField(index=True)
            full_repeat_sequence = pw.TextField(null=True)
            full_repeat_sequence_length = pw.FloatField(null=True)
            class Meta:
                database = db
                indexes = (
                    (('chrom', 'start', 'end', 'repeat_unit'), False),
                )

        class DatRecord(BedRecord):
            percent_matches = pw.FloatField(null=True)
            percent_indels = pw.FloatField(null=True)
            alignment_score = pw.FloatField(null=True)
            entropy = pw.FloatField(null=True)
            #full_repeat_sequence = pw.TextField(null=True)
            #full_repeat_sequence_length = pw.FloatField(null=True)
            left_flank = pw.TextField(null=True)
            right_flank = pw.TextField(null=True)


        if args.is_dat_file:
            Record = DatRecord
        elif args.is_bed_file:
            Record = BedRecord
        else:
            raise ValueError("Internal error")


        if Record.table_exists():
            Record.delete().execute()
        Record.create_table()
        with db.atomic():
            Record.bulk_create((
                Record(**generate_db_record(r._asdict())) for r in records
            ), batch_size=1000)

if __name__ == "__main__":
    main()











"""
        if args.merge:
            counters = defaultdict(int)
            for r in records:
                # for all overlapping records, select the one with the highest alignment score (for .dat records), or the longest one (for .bed)
                counters['total'] += 1
                current_best = r
                current_best__sorted_repeat_unit = "".join(sorted(r.repeat_unit))
                for r2 in Record.select().where(Record.chrom == r.chrom & Record.start <= r.end & Record.end >= r.start):
                    r2__sorted_repeat_unit = "".join(sorted(r2.repeat_unit))
                    has_alignment_scores = r2.alignment_score is not None and current_best.alignment_score is not None
                    if has_alignment_scores and r2.alignment_score > current_best.alignment_score:
                        current_best = r2
                        current_best__sorted_repeat_unit = r2__sorted_repeat_unit
                    elif not has_alignment_scores and r2.repeat_unit_length * r2.repeat_count > current_best.repeat_unit_length * current_best.repeat_count:
                        current_best = r2
                        current_best__sorted_repeat_unit = r2__sorted_repeat_unit
                    else:
                        if r2__sorted_repeat_unit == current_best__sorted_repeat_unit:
                            # replace r2 with current_best
                            r2.delete()
                            counters['deleted previous'] += 1
                        else:
                            counters['kept previous'] += 1

                if current_best is r:
                    Record.create(**r._asdict())
                    counters['added new'] += 1
                else:
                    counters['selected previous'] += 1

            print("Counters: ")
            for key, value in sorted(counters.items(), reverse=True, key=lambda x: x[1]):
                print("   ==> {} {}".format(value, key))
        else:
"""