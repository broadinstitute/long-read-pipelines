import json
from google.cloud import bigquery
import os
import sys
import tempfile


def export_to_bigquery(
    dict_records,
    schema,
    table_name,
    table_description=None,
    dataset_name="weisburd",
    project_id="broad-dsp-spec-ops",
    delete_table_if_exists=True,
    verbose=True):
    """Takes sequence (or iterator) of dictionaries and exports them to a BigQuery table.
    Works by writing all the dictionaries to a temp .json file (in NEWLINE_DELIMITED_JSON format), and then uploading
    this to BigQuery using the load_table_from_file API.

    Args:
        dict_records (iter): sequence or iterator over dictionaries
        schema (list): list of BigQuery fields to define in the table (each entry in the list should be a bigquery.SchemaField)
        table_name (str): BigQuery table name
        dataset_name (str): BigQuery dataset that will contain this table
        project_id (str): gcloud project that contains the dataset
        delete_table_if_exists (bool): whether to delete the table before creating it.
        verbose (bool): print info
    """

    client = bigquery.Client(project=project_id)
    table_ref = client.dataset(dataset_name).table(table_name)
    table = bigquery.Table(table_ref, schema=schema)

    if delete_table_if_exists and any([t._properties['tableReference']['tableId'] == table_name for t in client.list_tables(dataset_name)]):
        if verbose:
            print(f"Deleting table: {table}")
        client.delete_table(table)

    if verbose:
        print(f"Creating table: {table}")

    table.description = "\n\n".join([table_description or "", f"cwd:\n{os.getcwd()}", "cli:\n" + " ".join(sys.argv)])

    #table = client.update_table(table, ["description"])  # API request

    client.create_table(table)

    json_file = tempfile.NamedTemporaryFile('w', delete=False, suffix=".json")
    if verbose:
        print(f"Writing to {json_file.name}")

    for d in dict_records:
        json.dump(d, json_file)
        json_file.write("\n")

    json_file.close()

    if verbose:
        print(f"Wrote out {json_file.name}")

    if verbose:
        print(f"Loading .json data into {dataset_name}:{table_name}.")

    with open(json_file.name, 'rb') as f:
        job_config = bigquery.LoadJobConfig()
        job_config.source_format = bigquery.SourceFormat.NEWLINE_DELIMITED_JSON
        job_config.schema = schema
        job = client.load_table_from_file(f, destination=table_ref, job_config=job_config)

    try:
        job.result()  # Waits for table load to complete.
    except:
        print(f"ERRORS: {job.errors}")
    else:
        if verbose:
            print(f"Loaded {job.output_rows} rows into {dataset_name}:{table_name}.")

        os.remove(json_file.name)



# from itertools import islice, chain
# def batch(iterable, size):
#     sourceiter = iter(iterable)
#     while True:
#         batchiter = islice(sourceiter, size)
#         yield chain([next(batchiter)], batchiter)
#
# for records_batch in batch(dat_records, 10**4):
#     errors = client.insert_rows_json(table, (convert_record_to_dict(d) for d in tqdm.tqdm(records_batch, unit=" records"))
#     )
#
#     if errors:
#         for error in errors[:99]:
#             print(f"Error in record #{error['index']}: " + ", ".join([e['message'] for e in error['errors']]))
#         print("Exiting..")
#         break
# else:
#     print("Done")
#     #for key, value in sorted(counters.items(), key=lambda x: x[1], reverse=True):
#     #    print("%10d  %s" % (value, key))


def update_table_description(
    table_name,
    table_description,
    dataset_name="weisburd",
    project_id="broad-dsp-spec-ops"):

    client = bigquery.Client(project=project_id)
    table_ref = client.dataset(dataset_name).table(table_name)
    table = bigquery.Table(table_ref)
    table.description = table_description

    table = client.update_table(table, ["description"])  # API request
