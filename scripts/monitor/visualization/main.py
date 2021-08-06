#!
import pickle

from template import *
from week_one import *
from week_two import *

ds_meta_runtime = pickle.load(open("monitoring_dataset_50e68afe-a342-4357-9397-4de7b9300b10.metadata_runtime.pickle", "rb"))
ds_metrics = pickle.load(open("monitoring_dataset_50e68afe-a342-4357-9397-4de7b9300b10.metrics.pickle", "rb"))
# above is very temporary measure


def generate_individual_plots():
    pass


def generate_report_homepage(meta_runtime: pd.DataFrame, metrics: pd.DataFrame):

    html_template = HTML_STYLE+HTML_HEADING+generate_html_tables_section(table_html(meta_runtime, metrics))+HTML_VIOLIN+HTML_SCATTER
    with open('home.html', 'w') as home:
        #home.write('') #alternate
        home.write(html_template)


def table_html(meta_runtime: pd.DataFrame, metrics: pd.DataFrame):
    dicts = num_task_shards(meta_runtime)
    shard_count_table = pd.DataFrame(data=zip(dicts.keys(), dicts.values()), columns=["Task", "Shard Count"]).to_html(justify="center")
    shardless_table = shardless_task_table(meta_runtime, metrics).to_html(justify="center")
    return [shardless_table.replace("<td", '''<td style="text-align:center"'''), shard_count_table.replace("<td", '''<td style="text-align:center"''')]


def main(meta_runtime: pd.DataFrame, metrics: pd.DataFrame):
    generate_individual_plots()
    generate_report_homepage(meta_runtime, metrics)


if __name__ == 'main':
    main(ds_meta_runtime, ds_metrics)


def open_():  # for jupyter notebook debugging
    main(ds_meta_runtime, ds_metrics)
    os.system("open home.html")
