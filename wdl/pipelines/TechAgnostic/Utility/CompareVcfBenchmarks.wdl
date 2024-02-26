version 1.0

workflow CompareVcfBenchmarks {

    meta {
        description: "The purpose of this workflow is to create a table to facilitate a comparison of precision and sensitivity between different configurations (e.g. pipeline versions, chemistry changes) of samples vs. truth that have been obtained with the BenchmarkVCFs workflow.  After the CompareBenchmarks workflow is run, you can export the generated CSV table into a Google Sheets spreadsheet with automatic formatting using the ExportToGoogleSheets Colab Notebook (https://github.com/broadinstitute/palantir-workflows/blob/mg_benchmark_compare/BenchmarkVCFs/README_CompareBenchmarks.md)."
    }
    parameter_meta {
        sample_ids: "The names of one or more different samples. The comparisons will be made between configurations for each sample individually."
        configurations: "The labels for the different configurations that should be compared to each other."
        benchmark_summaries: "The output summary.tsv files from the BenchmarkVCFs workflow."
        stratifiers: "This input requires the same labels as the stratLabels input that has been passed to BenchmarkVCFs."
        include_counts: "If set to false, the resulting metrics will be Sensitivity, Precision and F-Measure. If set to true, the output table will also include the number of TP, FP and FN variants for each stratifier."
        generate_gc_plots: "If set to true, will generate plots of GC."
        order_of_samples: "If multiple different sample names are provided you can specify the order of those samples in the resulting table. Here, each sample name only has to be specified once, not once for each input VCF. If not specified, the order will be determined by the input files."
        order_of_configurations: "This input determines the order of the configurations in the resulting table. Just as above, each configuration only has to be specified once, not once for each input VCF. If not specified, the order will be determined by the input files."
        deltas: "If specified, the output table will contain columns that compares the results of different configurations to each other. Each delta column is defined by two entries in the array (referenced by zero-based index). If, for example there are three configurations A, B and C and you want to compare configurations B to A and C to A, provide the following data: [1, 0, 2, 0]. Note that you will want to define the order_of_configurations in this case to make sure that the indices refer to the correct configurations."
        mem_gb: "Optional input overriding the default memory."
        preemptible: "Optional input overriding the default number of preemptible attempts."
    }

    input {
        Array[String] sample_ids
        Array[String] configurations
        Array[File] benchmark_summaries
        Array[String]? stratifiers

        Boolean include_counts = true
        Boolean generate_gc_plots = false

        Array[String]? order_of_samples
        Array[String]? order_of_configurations
        Array[Int]? deltas

        Int? mem_gb
        Int? preemptible
    }

    call CompareBenchmarksTask {
        input:
            sample_ids = sample_ids,
            configurations = configurations,
            benchmark_summaries = benchmark_summaries,
            stratifiers = stratifiers,
            include_counts = include_counts,
            order_of_samples = order_of_samples,
            order_of_configurations = order_of_configurations,
            deltas = deltas,
            mem_gb = mem_gb,
            preemptible = preemptible
    }

    if (generate_gc_plots) {
        call CreateGCPlotsTask {
            input:
                sample_ids = sample_ids,
                configurations = configurations,
                benchmark_summaries = benchmark_summaries,
                order_of_samples = order_of_samples,
                order_of_configurations = order_of_configurations,
                mem_gb = mem_gb,
                preemptible = preemptible
        }
    }

    output {
        File comparison_csv = CompareBenchmarksTask.comparison_csv
        File raw_data = CompareBenchmarksTask.raw_data
        Array[File]? gc_plots = CreateGCPlotsTask.gc_plots
    }
}

task CreateGCPlotsTask {
    input {
        Array[String] sample_ids
        Array[String] configurations
        Array[File] benchmark_summaries

        Array[String]? order_of_samples
        Array[String]? order_of_configurations

        Int mem_gb = 4
        Int preemptible = 0
    }

    String order_of_samples_arg = if !defined(order_of_samples) then "" else "--order-of-samples"
    Array[String] order_of_samples_or_empty = select_first([order_of_samples, []])
    String order_of_configurations_arg = if !defined(order_of_configurations) then "" else "--order-of-configurations"
    Array[String] order_of_configurations_or_empty = select_first([order_of_configurations, []])

    command <<<
        set -xeuo pipefail

        source activate compare_benchmarks

        cat <<'EOF' > script.py
import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.family'] = 'serif'

def get_value_from_table(data, sample_id, configuration, stratifier, var_type, column):
    return data.query('sample_id == @sample_id and configuration == @configuration and Stratifier == @stratifier and Type == @var_type').iloc[0][column]

def calculate_metrics(data, unique_sample_ids, unique_configurations, stratifiers):
    recalculated_data = pd.DataFrame(columns=['sample_id', 'configuration', 'Stratifier', 'Type', 'TP', 'FP', 'FN', 'Precision', 'Sensitivity', 'F-Measure'])
    for sample_id in unique_sample_ids:
        for configuration in unique_configurations:
            for stratifier in stratifiers:
                tp, fp, fn = dict(), dict(), dict()
                tp['SNP'] = get_value_from_table(data, sample_id, configuration, stratifier, 'SNP', 'TP_Base')
                fp['SNP'] = get_value_from_table(data, sample_id, configuration, stratifier, 'SNP', 'FP')
                fn['SNP'] = get_value_from_table(data, sample_id, configuration, stratifier, 'SNP', 'FN')
                tp['INDEL'] = get_value_from_table(data, sample_id, configuration, stratifier, 'INDEL', 'TP_Base')
                fp['INDEL'] = get_value_from_table(data, sample_id, configuration, stratifier, 'INDEL', 'FP')
                fn['INDEL'] = get_value_from_table(data, sample_id, configuration, stratifier, 'INDEL', 'FN')
                tp['all'] = tp['SNP'] + tp['INDEL']
                fp['all'] = fp['SNP'] + fp['INDEL']
                fn['all'] = fn['SNP'] + fn['INDEL']

                for var_type in ['SNP', 'INDEL', 'all']:
                    recalculated_data = recalculated_data.append({
                        'sample_id': sample_id,
                        'configuration': configuration,
                        'Stratifier': stratifier,
                        'Type': var_type,
                        'TP': tp[var_type],
                        'FP': fp[var_type],
                        'FN': fn[var_type],
                        'Precision': tp[var_type]/(tp[var_type] + fp[var_type]) if tp[var_type] + fp[var_type] > 0 else np.nan,
                        'Sensitivity': tp[var_type]/(tp[var_type] + fn[var_type]) if tp[var_type] + fn[var_type] > 0 else np.nan,
                        'F-Measure': tp[var_type]/(tp[var_type] + 0.5*(fp[var_type] + fn[var_type])) if tp[var_type] + fp[var_type] + fn[var_type] > 0 else np.nan
                        }, ignore_index=True)
    return recalculated_data

def plot_sample(data, i_sample, sample_id, unique_configurations, stratifiers):
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for column, var_type in enumerate(['SNP', 'INDEL']):
        for row, metric in enumerate(['Sensitivity', 'Precision']):
            ax = axes[row, column]
            X = np.arange(len(stratifiers))
            for configuration in unique_configurations:
                ax.plot(X, [get_value_from_table(data, sample_id, configuration, stratifier, var_type, metric) for stratifier in stratifiers], label=configuration)

            ax.set_xticks(X)
            ax.set_xticklabels([stratifier.replace('gc', '') for stratifier in stratifiers], rotation=45)
            ax.set_xlabel('GC bin')
            ax.set_ylabel(metric)
            ax.set_title(var_type, zorder=0)
    axes[0, 0].legend()
    fig.suptitle(sample_id)
    plt.tight_layout()

    fig.savefig(f'gc_plot_{i_sample}_{sample_id}.png', dpi=100)


def main(sample_ids, configurations, summaries, order_of_samples, order_of_configurations):
    if len(sample_ids) != len(configurations) or len(sample_ids) != len(summaries):
        raise RuntimeError('The number of sample_ids, configurations, and summaries tables must be equal.')

    samples_data = []
    for i in range(len(sample_ids)):
        sample_data = pd.read_csv(summaries[i])

        # Filter out everything other than SNP or INDEL rows, and stratifiers starting with "gc"
        sample_data = sample_data.loc[((sample_data['Type'] == 'SNP') | (sample_data['Type'] == 'INDEL')) & (sample_data['Stratifier'].str.startswith("gc"))]

        # Add sample_id and configuration names
        sample_data['sample_id'] = sample_ids[i]
        sample_data['configuration'] = configurations[i]
        samples_data.append(sample_data)
    data = pd.concat(samples_data)

    if order_of_samples is None:
        unique_sample_ids = data['sample_id'].unique()
    else:
        unique_sample_ids = order_of_samples

    if order_of_configurations is None:
        unique_configurations = data['configuration'].unique()
    else:
        unique_configurations = order_of_configurations

    gc_stratifiers = data['Stratifier'].unique().astype(str)
    gc_stratifiers.sort()

    data = calculate_metrics(data, unique_sample_ids, unique_configurations, gc_stratifiers)

    for i_sample_id, sample_id in enumerate(unique_sample_ids):
        plot_sample(data, i_sample_id, sample_id, unique_configurations, gc_stratifiers)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a table to compare output of BenchmarkVCFs.')
    parser.add_argument('--order-of-samples', type=str, nargs='+', help='Order of samples. If not specified, the order will be the same as the supplied inputs.')
    parser.add_argument('--order-of-configurations', type=str, nargs='+', help='Order of configurations. If not specified, the order will be the same as the supplied inputs.')
    required_named = parser.add_argument_group('Required named arguments')
    required_named.add_argument('--sample-ids', required=True, type=str, nargs='+')
    required_named.add_argument('--configurations', required=True, type=str, nargs='+')
    required_named.add_argument('--summaries', required=True, type=str, nargs='+')
    args = parser.parse_args()
    main(args.sample_ids, args.configurations, args.summaries, args.order_of_samples, args.order_of_configurations)
EOF
        python script.py --sample-ids ~{sep=' ' sample_ids} --configurations ~{sep=' ' configurations} --summaries ~{sep=' ' benchmark_summaries} ~{order_of_samples_arg} ~{sep=' ' order_of_samples_or_empty} ~{order_of_configurations_arg} ~{sep=' ' order_of_configurations_or_empty}
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/benchmark_vcfs/compare_benchmarks:1.1.0"
        preemptible: preemptible
        memory: mem_gb + " GB"
        disks: "local-disk 20 HDD"
    }

    output {
        Array[File] gc_plots = glob("gc_plot_*.png")
    }
}



task CompareBenchmarksTask {
    input {
        Array[String] sample_ids
        Array[String] configurations
        Array[File] benchmark_summaries
        Array[String]? stratifiers

        Boolean include_counts

        Array[String]? order_of_samples
        Array[String]? order_of_configurations
        Array[Int]? deltas

        Int mem_gb = 4
        Int preemptible = 0
    }

    String stratifiers_arg = if !defined(stratifiers) then "" else "--stratifiers"
    Array[String] stratifiers_or_empty = select_first([stratifiers, []])

    String order_of_samples_arg = if !defined(order_of_samples) then "" else "--order-of-samples"
    Array[String] order_of_samples_or_empty = select_first([order_of_samples, []])
    String order_of_configurations_arg = if !defined(order_of_configurations) then "" else "--order-of-configurations"
    Array[String] order_of_configurations_or_empty = select_first([order_of_configurations, []])

    String deltas_arg = if !defined(deltas) then "" else "--deltas"
    Array[Int] deltas_or_empty = select_first([deltas, []])

    command <<<
        set -xeuo pipefail

        source activate compare_benchmarks

        cat <<'EOF' > script.py
import argparse
import numpy as np
import pandas as pd

# The purpose of this class is to significantly simplify writing a TSV (or otherwise separated) table.
# For example, when passing an array of strings to cells(arr), the elements will be written separated
# by the separator sep, which is equivalent to file.write(sep.join(arr)). However, the methods of this class allow
# for chaining of these write commands to easily concatenate multiple calls, which can be helpful for writing
# custom-formatted tables, e.g. output.cells(arr1).sep().cells(arr2).new_line() instead of having to issue multiple
# file.write() calls.
class ChainableOutput:
    def __init__(self, file, sep:str):
        self.file = file
        self.separator = sep

    def cells(self, cells:list):
        self.file.write(self.separator.join([str(cell) for cell in cells]))
        return self
    def sep(self):
        self.file.write(self.separator)
        return self
    def new_line(self):
        self.file.write('\n')
        return self

def write_header(output:ChainableOutput, unique_sample_ids, unique_configurations, deltas, include_counts):
    # This function will write the header for the comparison table. In general, the layout looks like this:
    # | | | |               Precision               |              Sensitivity              | ...
    # | | | |      sample1      |      sample2      |      sample1      |      sample2      | ...
    # | | | | config1 | config2 | config1 | config2 | config1 | config2 | config1 | config2 | ...

    # If include_counts is True then TP, FP and FN will be written before Precision, Sensitivity and F-Measure

    # Line 1
    output.cells([''] * 4).sep()
    if include_counts:
        output.cells(['TP'] + [''] * (len(unique_sample_ids) * (len(unique_configurations) + len(deltas)) - 1)).sep()
        output.cells(['FP'] + [''] * (len(unique_sample_ids) * (len(unique_configurations) + len(deltas)) - 1)).sep()
        output.cells(['FN'] + [''] * (len(unique_sample_ids) * (len(unique_configurations) + len(deltas)) - 1)).sep()
    output.cells(['Precision'] + [''] * (len(unique_sample_ids) * (len(unique_configurations) + len(deltas)) - 1)).sep()
    output.cells(['Sensitivity'] + [''] * (len(unique_sample_ids) * (len(unique_configurations) + len(deltas)) - 1)).sep()
    output.cells(['F-Measure'] + [''] * (len(unique_sample_ids) * (len(unique_configurations) + len(deltas)) - 1)).sep()
    output.cells([''] * 2).new_line()

    # Line 2
    output.cells([''] * 4)
    for _ in (['TP', 'FP', 'FN'] if include_counts else []) + ['Precision', 'Sensitivity', 'F-Measure']:
        for sample_id in unique_sample_ids:
            output.sep().cells([sample_id] + [''] * (len(unique_configurations) + len(deltas) - 1))
    output.sep().cells([''] * 2).new_line()

    # Line 3
    output.cells([''] * 4)
    for metric in (['TP', 'FP', 'FN'] if include_counts else []) + ['Precision', 'Sensitivity', 'F-Measure']:
        for sample_id in unique_sample_ids:
            # Write configurations themselves
            for configuration in unique_configurations:
                output.sep().cells([configuration])
            # Write deltas
            for delta_pair in deltas:
                if metric in ['TP', 'FP', 'FN']:
                    output.sep().cells([f'delta({unique_configurations[delta_pair[1]]}-{unique_configurations[delta_pair[0]]})'])
                else:
                    output.sep().cells([f'delta%({unique_configurations[delta_pair[1]]}-{unique_configurations[delta_pair[0]]})'])
    output.sep().cells([''] * 2).new_line()

def get_value_from_table(data, sample_id, configuration, stratifier, var_type, column):
    try:
        return data.query('sample_id == @sample_id and configuration == @configuration and Stratifier == @stratifier and Type == @var_type').iloc[0][column]
    except IndexError as e:
        raise RuntimeError(f'Failed querying table for sample_id: {sample_id}, configuration: {configuration}, stratifier {stratifier}, var_type: {var_type}, column: {column}. Make sure that the set of stratifiers and samples in the order_of_ arguments exactly matches the set of stratifiers and samples used for BenchmarkVCFs.')

def write_stratifier(output:ChainableOutput, stratifier:str, data:pd.DataFrame, unique_sample_ids:list, unique_configurations:list, deltas:list, include_counts):
    for var_type in ['SNP', 'INDEL', 'all']:
        # Placeholder to print percentage of genome for each stratification
        output.cells(['', '']).sep()
        # Only print stratifier name in the first row
        output.cells([stratifier if var_type == 'SNP' else '', var_type])
        for metric in (['TP', 'FP', 'FN'] if include_counts else []) + ['Precision', 'Sensitivity', 'F-Measure']:
            for sample_id in unique_sample_ids:
                for configuration in unique_configurations:
                    output.sep().cells(['{:.5f}'.format(get_value_from_table(data, sample_id, configuration, stratifier, var_type, metric))])
                for delta_pair in deltas:
                    base_value = get_value_from_table(data, sample_id, unique_configurations[delta_pair[0]], stratifier, var_type, metric)
                    current_value = get_value_from_table(data, sample_id, unique_configurations[delta_pair[1]], stratifier, var_type, metric)
                    if metric in ['TP', 'FP', 'FN']:
                        delta = int(current_value) - int(base_value)
                        output.sep().cells(['{}'.format(delta)])
                    else:
                        delta_pct = (current_value - base_value) / base_value
                        output.sep().cells(['{:.2%}'.format(delta_pct)])
        output.sep().cells([var_type, stratifier if var_type == 'SNP' else '']).new_line()

def calculate_metrics(data, unique_sample_ids, unique_configurations, stratifiers):
    recalculated_data = pd.DataFrame(columns=['sample_id', 'configuration', 'Stratifier', 'Type', 'TP', 'FP', 'FN', 'Precision', 'Sensitivity', 'F-Measure'])
    for sample_id in unique_sample_ids:
        for configuration in unique_configurations:
            for stratifier in stratifiers:
                tp, fp, fn = dict(), dict(), dict()
                tp['SNP'] = get_value_from_table(data, sample_id, configuration, stratifier, 'SNP', 'TP_Base')
                fp['SNP'] = get_value_from_table(data, sample_id, configuration, stratifier, 'SNP', 'FP')
                fn['SNP'] = get_value_from_table(data, sample_id, configuration, stratifier, 'SNP', 'FN')
                tp['INDEL'] = get_value_from_table(data, sample_id, configuration, stratifier, 'INDEL', 'TP_Base')
                fp['INDEL'] = get_value_from_table(data, sample_id, configuration, stratifier, 'INDEL', 'FP')
                fn['INDEL'] = get_value_from_table(data, sample_id, configuration, stratifier, 'INDEL', 'FN')
                tp['all'] = tp['SNP'] + tp['INDEL']
                fp['all'] = fp['SNP'] + fp['INDEL']
                fn['all'] = fn['SNP'] + fn['INDEL']

                for var_type in ['SNP', 'INDEL', 'all']:
                    recalculated_data = recalculated_data.append({
                        'sample_id': sample_id,
                        'configuration': configuration,
                        'Stratifier': stratifier,
                        'Type': var_type,
                        'TP': tp[var_type],
                        'FP': fp[var_type],
                        'FN': fn[var_type],
                        'Precision': tp[var_type]/(tp[var_type] + fp[var_type]) if tp[var_type] + fp[var_type] > 0 else np.nan,
                        'Sensitivity': tp[var_type]/(tp[var_type] + fn[var_type]) if tp[var_type] + fn[var_type] > 0 else np.nan,
                        'F-Measure': tp[var_type]/(tp[var_type] + 0.5*(fp[var_type] + fn[var_type])) if tp[var_type] + fp[var_type] + fn[var_type] > 0 else np.nan
                        }, ignore_index=True)
    return recalculated_data



def main(sample_ids, configurations, summaries, stratifiers, order_of_samples, order_of_configurations, deltas_array, include_counts):
    if len(sample_ids) != len(configurations) or len(sample_ids) != len(summaries):
        raise RuntimeError('The number of sample_id, configurations, and summary tables must be equal.')
    if deltas_array is None:
        deltas_array = []
    if len(deltas_array) % 2 != 0:
        raise RuntimeError('The number of indices in the delta argument must be even. Please use --help or check the documentation on how to use this argument.')

    deltas = [(int(deltas_array[i]), int(deltas_array[i+1])) for i in range(0, len(deltas_array), 2)]

    samples_data = []
    for i in range(len(sample_ids)):
        sample_data = pd.read_csv(summaries[i])

        # Filter out everything other than SNP or INDEL rows
        sample_data = sample_data.loc[(sample_data['Type'] == 'SNP') | (sample_data['Type'] == 'INDEL')]

        # Add sample_id and configuration names
        sample_data['sample_id'] = sample_ids[i]
        sample_data['configuration'] = configurations[i]
        samples_data.append(sample_data)
    data = pd.concat(samples_data)
    data = data.fillna({'Stratifier': 'all'})

    if order_of_samples is None:
        unique_sample_ids = data['sample_id'].unique()
    else:
        unique_sample_ids = order_of_samples

    if order_of_configurations is None:
        unique_configurations = data['configuration'].unique()
    else:
        unique_configurations = order_of_configurations

    if stratifiers is None:
        stratifiers = data['Stratifier'].unique()
    else:
        stratifiers = ['all'] + stratifiers

    data = calculate_metrics(data, unique_sample_ids, unique_configurations, stratifiers)

    data.to_csv('raw_data.csv')

    with open('comparison.csv', 'w') as output_file:
        chainable_output = ChainableOutput(output_file, ',')
        write_header(chainable_output, unique_sample_ids, unique_configurations, deltas, include_counts)
        for stratifier in stratifiers:
            write_stratifier(chainable_output, stratifier, data, unique_sample_ids, unique_configurations, deltas, include_counts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a table to compare output of BenchmarkVCFs.')
    parser.add_argument('--stratifiers', type=str, nargs='*', help='Explicitly specify the stratifiers that have to be present in all samples. "all" will automatically be added. If not specified, the stratifiers will be inferred. If specified, this argument also defines the order of the stratifiers.')
    parser.add_argument('--order-of-samples', type=str, nargs='+', help='Order of samples. If not specified, the order will be the same as the supplied inputs.')
    parser.add_argument('--order-of-configurations', type=str, nargs='+', help='Order of configurations. If not specified, the order will be the same as the supplied inputs.')
    parser.add_argument('--deltas', type=str, nargs='+', help='A list of configuration (zero-based) indices to compare. E.g. for comparing configurations 0 to 1 and 0 to 2, pass the values 0 1 0 2.')
    parser.add_argument('--include-counts', action='store_true', help='If set, include the TP/FP/FN counts in the output table.')
    required_named = parser.add_argument_group('Required named arguments')
    required_named.add_argument('--sample-ids', required=True, type=str, nargs='+')
    required_named.add_argument('--configurations', required=True, type=str, nargs='+')
    required_named.add_argument('--summaries', required=True, type=str, nargs='+')
    args = parser.parse_args()
    main(args.sample_ids, args.configurations, args.summaries, args.stratifiers, args.order_of_samples, args.order_of_configurations, args.deltas, args.include_counts)
EOF
        python script.py --sample-ids ~{sep=' ' sample_ids} --configurations ~{sep=' ' configurations} --summaries ~{sep=' ' benchmark_summaries} ~{stratifiers_arg} ~{sep=' ' stratifiers_or_empty} ~{order_of_samples_arg} ~{sep=' ' order_of_samples_or_empty} ~{order_of_configurations_arg} ~{sep=' ' order_of_configurations_or_empty} ~{deltas_arg} ~{sep=' ' deltas_or_empty} ~{true="--include-counts" false="" include_counts}
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/benchmark_vcfs/compare_benchmarks:1.1.0"
        preemptible: preemptible
        memory: mem_gb + " GB"
        disks: "local-disk 20 HDD"
    }

    output {
        File comparison_csv = "comparison.csv"
        File raw_data = "raw_data.csv"
    }
}
