from typing import Dict
import matplotlib.pyplot as plt
from bokeh.models import Title
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import pandas as pd
import pandas_bokeh
from operator import truediv
import copy
# 2 main functions are run and setup.
# Once combine is run, class instance becomes solely dedicated to analyzing the multiple datasets. All unrelated methods are deleted
#


class ResourceUsage:
    instances = 0  # count number of class instances (can't find a way to update this when an instance is deleted)
    _level = 0
    # _level will be used later to determine what sort of suggestions can be automatically made based on the given data
    # the more datasets the higher the level.
    # For example, with only one dataset there is no way to determine how much variation there is.

    def __init__(self, meta_runtime: pd.DataFrame, metrics: pd.DataFrame, name: str) -> None:
        """
        initialization step: bring meta, runtime, and metrics data into scope
        :param meta_runtime: meta_runtime table
        :param metrics: metrics table
        :param name: name for this instance (used to identify which files are created by this specific instance)
        """
        self.name = name
        self._id = ResourceUsage.instances
        self.meta_runtime = meta_runtime
        self.metrics = metrics
        tasks_list = self.num_task_shards().keys()
        self.avg_max_cpu_usage = pd.DataFrame(index=["avg_max_value", "requested"], columns=tasks_list)  # set in violin plot functions
        self.max_mem_usage = pd.DataFrame(index=["max_value", "requested"], columns=tasks_list)
        self.max_disk_usage = pd.DataFrame(index=["max_value", "requested"], columns=tasks_list)
        ResourceUsage.instances += 1
        self.setup()

    def setup(self) -> None:  # maybe do this as part of initialization
        """
        Load all files and information for the webpage
        :return: None
        """
        self.generate_plot_files()
        self.generate_basic_suggestions()

    def run(self) -> None:
        """
        display the webpage
        :return: None
        """
        self.generate_report_homepage()

    def generate_report_homepage(self) -> None:  # finish
        """
        Set up webpage
        :return: None
        """

        pass
        # jekyll stuff

    def generate_plot_files(self) -> None:
        """
        Create all files and put them where other functions can find them
        :return: None
        """
        tasks = self.num_task_shards()
        unsharded_tasks = dict([(k, v) for k, v in tasks.items() if v == 0])
        sharded_tasks = dict([(k, v) for k, v in tasks.items() if v != 0])
        for task in unsharded_tasks:
            print(f"Working on Task: {task}")
            self.plot_aggregate(task, -1)
        for task in sharded_tasks:
            self.save_bokeh(task)

    def generate_basic_suggestions(self):
        tasks = self.num_task_shards().keys()
        with open('suggestions.txt', 'a+') as suggestion_text:
            suggestion_text.truncate(0)
            for task in tasks:
                # thresholds can be tweaked later
                # cpu is trickier, since it can have multiple cores and only requested # of cores is given
                suggestion_text.write(f"{task}'s allocated CPU can be {'reduced' if self.avg_max_cpu_usage.loc['avg_max_value'][task] > 90 else 'increased'}. \n")
                suggestion_text.write(f"{task}'s allocated Memory can be {'reduced' if self.max_mem_usage.loc['max_value'][task] > 90 else 'increased'}. \n")
                suggestion_text.write(f"{task}'s allocated Disk can be {'reduced' if self.max_disk_usage.loc['max_value'][task] > 90 else 'increased'}. \n")
        pass

    def combine(self, *args):  # UNTESTED: NO IDEA IF THIS WORKS
        """
        Combines multiple class instances (with their workflow data)
        Meant for class instances with the same workflow
        :param args: other class instances to combine- for data on multiple workflows
        :return: None
        """
        Frames = list(self.meta_runtime)
        for class_instance in args:
            if class_instance.meta_runtime.loc[0, "meta_workflow_name"] == self.meta_runtime.loc[0, "meta_workflow_name"]:
                #
                # NOTE TO SELF: it's easier to just check that the workflow names are the same
                #
                Frames.append(class_instance.meta_runtime)  # this is ok as long as I check the meta_workflow_id later
            else:
                raise Exception("Not the same workflow, data cannot be combined")
        self.meta_runtime = pd.concat(Frames)
        # Some methods will no longer be compatible after the datasets are combined

    def num_task_shards(self) -> Dict[str, int]:
        """
        Returns dict with task names as key and shard count as value
        :return: dict of task:shard_count
        """

        task_and_shard_idx = self.meta_runtime.loc[:, ['runtime_task_call_name', 'runtime_shard']]
        res = dict()
        for task in task_and_shard_idx['runtime_task_call_name'].unique():
            res[task] = 1 + max(task_and_shard_idx
                                [task_and_shard_idx['runtime_task_call_name'] == task]['runtime_shard'])
        return res

    def disk_requested_for_task(self, task: str) -> pd.DataFrame:
        """
        Disk requested for a task (only successful runs, if sharded display all shards)
        :param task: task name
        :return: dataframe of disk requested for task
        """
        columns_to_select = ["runtime_task_call_name", "runtime_shard", "runtime_disk_total_gb",
                             "meta_disk_total_gb", "meta_execution_status"]
        renaming_scheme = {"runtime_task_call_name": "task",
                           "runtime_shard": "shard_idx"}
        df = self.meta_runtime.loc[:, columns_to_select].rename(columns=renaming_scheme)
        df = df[df["meta_execution_status"] == "Done"]
        if task in df['task'].unique():
            return df[df["task"] == task].sort_values(by=['shard_idx'])
        else:
            raise KeyError("Task name does not exist")

    def mem_requested_for_task(self, task: str) -> pd.DataFrame:
        """
        Memory requested for a task (only successful runs, if sharded display all shards)
        :param task: task name
        :return: dataframe of memory meta for task
        """
        columns_to_select = ["runtime_task_call_name", "runtime_shard", "runtime_mem_total_gb",
                             "meta_mem_total_gb", "meta_execution_status"]
        renaming_scheme = {"runtime_task_call_name": "task", "runtime_shard": "shard"}
        df = self.meta_runtime.loc[:, columns_to_select].rename(columns=renaming_scheme)
        df = df[df["meta_execution_status"] == "Done"]
        if task in df['task'].unique():
            return df[df["task"] == task].sort_values(by=['shard'])
        else:
            raise KeyError("Task name does not exist")

    def plot_task_disk_requested(self, task: str) -> plt.Figure:
        """
        Useful for tasks split into many shards
        Will only plot the first disk in array, use in conjunction with disk_requested_for_task
        :param task: task name
        :return: plot
        """
        x = self.disk_requested_for_task(task).loc[:, "shard"]
        y1 = self.disk_requested_for_task(task).loc[:, "runtime_disk_total_gb"]\
            .map(lambda arr: arr[0]).astype("float64", copy=True)
        y2 = self.disk_requested_for_task(task).loc[:, "meta_disk_total_gb"]\
            .map(lambda arr: arr[0]).astype("float64", copy=True)
        plt.plot(x, y1)
        plt.plot(x, y2)
        plt.legend(["Runtime Disk Usage(gb)", "Meta Disk Usage(gb)"])
        plt.show()
        # plt.savefig(f"{__name__}.pdf")
        return plt.figure()

    def plot_task_mem_requested(self, task: str) -> plt.Figure:
        """
        Useful for with tasks split into many shards
        Use in conjunction with mem_requested_for_task
        :param task: task name
        :return: plot
        """
        x = self.mem_requested_for_task(task).loc[:, "shard"]
        y1 = self.mem_requested_for_task(task).loc[:, "runtime_mem_total_gb"].astype("float64", copy=True)
        y2 = self.mem_requested_for_task(task).loc[:, "meta_mem_total_gb"].astype("float64", copy=True)
        plt.plot(x, y1)
        plt.plot(x, y2)
        plt.legend(["Runtime Memory Usage(gb)", "Meta Memory Usage(gb)"])
        plt.show()
        # plt.savefig(f"{__name__}.pdf")
        return plt.figure()

    def task_input_size(self, task_name_or_index, shard_index=None) -> float:
        """
        Summation of a task's input file sizes
        Use:
        task_input_size(row index of task/shard), or task_input_size(task name, shard index)
        :param task_name_or_index: task name or index (details above)
        :param shard_index: shard index if task name is given
        :return:
        """
        pos = 0
        if shard_index is not None:
            tasks = self.meta_runtime.loc[:, "runtime_task_call_name"]  # find row index
            shards = self.meta_runtime.loc[:, "runtime_shard"]
            for i in range(0, tasks.size):
                if tasks[i] == task_name_or_index and shards[i] == shard_index:
                    pos = i
        else:
            if isinstance(task_name_or_index, int):
                pos = task_name_or_index  # if only task_name is supplied it's treated as the row index
            else:
                raise Exception(f"Need an integer row index but received {task_name_or_index}")
        input_ds = pd.DataFrame(data=self.meta_runtime.loc[pos, "meta_inputs"])
        return input_ds[input_ds['type'].astype(str) == "file"].loc[:, "value"].astype(float).sum()

    def cpu_use_violin(self, instance_id: int) -> plt.Figure:
        """
        Violin plot for CPU use
        :param instance_id: id of vm a shard is running on
        :return: plot
        """
        by_id = self.metrics[self.metrics["metrics_instance_id"] == instance_id]\
            .sort_values("metrics_timestamp").loc[:, "metrics_cpu_used_percent"]
        data = list()
        if len(by_id) == 0:
            print(f"Instance id: {instance_id} for task: {self.get_shard_from_instance_id(instance_id)[0]} shard: {self.get_shard_from_instance_id(instance_id)[1]} was not found in metrics")
            fig = plt.figure()
            axis = fig.add_axes([0, 0, 1, 1])
            name = self.get_shard_from_instance_id(instance_id)
            axis.set_title(f"FAILED: Empty CPU Usage Plot for Failed Shard {name[0]}, shard {name[1]}")
            plt.show()
            return fig

        for i in range(len(by_id.iloc[0])):
                data.append(pd.Series(list(e[i] for e in by_id)))
        fig = plt.figure()
        axis = fig.add_axes([0, 0, 1, 1])
        name = self.get_shard_from_instance_id(instance_id)

        if len(by_id.iloc[0]) > 10:
            avg_data = list()

            for i in range(len(by_id)):
                avg_data.append(sum(list(e for e in by_id.iloc[i]))/len(by_id.iloc[0]))

            _ = axis.violinplot(avg_data, showmeans=True, showextrema=True, showmedians=True)
            axis.set_xlabel("Average CPU Core")
            axis.set_ylabel("CPU % used")
            axis.set_title(f"Average CPU Usage Across all {len(by_id.iloc[0])} Cores for {name[0]}, shard {name[1]}")
        else:
            _ = axis.violinplot(data, showmeans=True, showextrema=True, showmedians=True)
            axis.set_xlabel("CPU core #")
            axis.set_ylabel("CPU % used")
            axis.set_title(f"CPU Usage for {name[0]}, shard {name[1]}")
        plt.show()
        try:
            averaged_data = [sum(each)/len(each) for each in data]
            self.avg_max_cpu_usage.loc["avg_max_value"][name[0]] = max(averaged_data)
        except ValueError:
            print(name)
            print(self.avg_max_cpu_usage.columns)
            # print(data[0])
            raise
        return fig

    def get_instance_id(self, task: str, shard_idx: int) -> numpy.int64:
        """
        Gets instance ID given a shard- if task is unsharded, use -1 as shard index
        :param task: task name
        :param shard_idx: shard index #
        :return: instance ID
        """
        tasks = self.meta_runtime.loc[:, "runtime_task_call_name"]  # find row index
        shards = self.meta_runtime.loc[:, "runtime_shard"]
        pos = None
        for i in range(len(shards)):
            if tasks[i] == task and shards[i] == shard_idx:
                pos = i
        try:
            return self.meta_runtime.loc[pos, "runtime_instance_id"]
        except KeyError:
            print(f"task {task} shard {shard_idx} was not found")

    def cpu_plot(self, task: str, shard_idx: int) -> None:
        """
        Combining cpu_use_violin and get_instance_id (is this considered a wrapper function? maybe not)
        :param task: task name
        :param shard_idx: shard index #
        :return: plot
        """
        self.cpu_use_violin(self.get_instance_id(task, shard_idx))

    def mem_use_violin(self, instance_id: int) -> plt.Figure:
        """
        Violin plot for memory use
        :param instance_id: instance ID
        :return: plot
        """
        by_id = self.metrics[self.metrics["metrics_instance_id"] == instance_id]\
            .sort_values("metrics_timestamp").loc[:, "metrics_mem_used_gb"]
        row = self.get_runtime_row_from_instance_id(instance_id)
        requested_mem = float(row["meta_mem_total_gb"])
        data = list()
        data.append(pd.Series((e / requested_mem) * 100 for e in by_id))
        fig = plt.figure()
        axis = fig.add_axes([0, 0, 1, 1])
        if len(data[0]) == 0:
            # print(f"Instance id: {instance_id} for task: {self.get_shard_from_instance_id(instance_id)[0]} shard: {self.get_shard_from_instance_id(instance_id)[1]} was not found in metrics")
            fig = plt.figure()
            axis = fig.add_axes([0, 0, 1, 1])
            name = self.get_shard_from_instance_id(instance_id)
            axis.set_title(f"FAILED: Empty Memory Usage Plot for Failed Shard {name[0]}, shard {name[1]}")
            plt.show()
            return fig
        _ = axis.violinplot(data, showmeans=True, showextrema=True, showmedians=True)
        name = self.get_shard_from_instance_id(instance_id)
        axis.set_xlabel("mem")
        axis.set_ylabel("% GB used")
        axis.set_title(f"Memory Usage for {name[0]}, shard {name[1]}")
        plt.show()
        max = data[0].max()
        self.max_mem_usage.loc["max_value"][name[0]] = max  # max(data[0])
        self.max_mem_usage.loc["requested"][name[0]] = requested_mem
        return fig

    # @staticmethod
    def __get_first_val(self, row, key) -> float:
        """
        After getting a list from a row using the key, returns first value in list
        :param row: row
        :param key: key
        :return: float
        """
        return row[key][0]

    def disk_use_violin(self, instance_id: int) -> plt.Figure:
        """
        Violin plot for disk use
        :param instance_id: instance ID table
        :return: plot
        """
        by_id = self.metrics[self.metrics["metrics_instance_id"] == instance_id]\
            .sort_values("metrics_timestamp").loc[:, "metrics_disk_used_gb"]
        row = self.get_runtime_row_from_instance_id(instance_id)
        requested_disk = float(self.__get_first_val(row, 'meta_disk_total_gb'))
        data = list()
        data.append(pd.Series((float(e[0]) / requested_disk) * 100 for e in by_id))
        fig = plt.figure()
        axis = fig.add_axes([0, 0, 1, 1])
        if len(data[0]) == 0:
            # print(f"Instance id: {instance_id} for task: {self.get_shard_from_instance_id(instance_id)[0]} shard: {self.get_shard_from_instance_id(instance_id)[1]} was not found in metrics")
            fig = plt.figure()
            axis = fig.add_axes([0, 0, 1, 1])
            name = self.get_shard_from_instance_id(instance_id)
            axis.set_title(f"FAILED: Empty Disk Usage Plot for Failed Shard {name[0]}, shard {name[1]}")
            plt.show()
            return fig
        _ = axis.violinplot(data, showmeans=True, showextrema=True, showmedians=True)
        name = self.get_shard_from_instance_id(instance_id)
        axis.set_xlabel("disk")
        axis.set_ylabel("% GB used")
        axis.set_title(f"Disk Usage for {name[0]}, shard {name[1]}")
        plt.show()
        max = data[0].max()
        self.max_disk_usage.loc["max_value"][name[0]] = max
        self.max_disk_usage.loc["requested"][name[0]] = requested_disk

        return fig

    def get_shard_from_instance_id(self, instance_id: int) -> list:
        """
        Gets shard (task + shard index) from instance ID
        :param instance_id: instance ID
        :return: list of [task name, shard index #]
        """
        temp = self.meta_runtime[self.meta_runtime["runtime_instance_id"] == instance_id]
        # temp = temp[temp["meta_"]]
        name = temp["runtime_task_call_name"].to_numpy()
        idx = temp["runtime_shard"].to_numpy()
        return [name[0], idx[0]]

    def get_runtime_row_from_instance_id(self, instance_id: int) -> pd.DataFrame:
        """
        Gets the row for an instance ID in meta_runtime
        :param instance_id: instance ID
        :return: row
        """
        return self.meta_runtime[self.meta_runtime["runtime_instance_id"] == instance_id].iloc[0, :]  # Taking first row because the monitoring logger sometimes logs duplicate rows

    def get_shard_instance_id_list(self, task: str) -> list:
        """
        Gets a list of instance IDs for each shard of a given task
        :param task: task name
        :return: list of instance IDs
        """
        return list(self.meta_runtime[self.meta_runtime["runtime_task_call_name"] == task]
                    .sort_values(by=["runtime_shard"]).loc[:, "runtime_instance_id"])

    def avg_mem_across_shards(self, task: str):
        """
        Plots average memory use for each shard of a given task
        :param task: task name
        :return: plot
        """
        # pandas_bokeh.output_file(f"avg_mem_across_shards_{task}")
        averages = pd.Series([], dtype="float64")
        id_list = self.get_shard_instance_id_list(task)
        faster = self.metrics.loc[:, "metrics_instance_id"].astype('str').tolist()
        for instance_id in id_list:

            if str(instance_id) not in faster:  # 3870777233638278011 doesn't exist in metrics table
                print(f"the instance id: {instance_id} in meta_runtime does not exist in metrics")
                continue

            by_id = self.metrics[self.metrics["metrics_instance_id"] == instance_id]\
                        .sort_values("metrics_timestamp").loc[:, "metrics_mem_used_gb"]
            row = self.get_runtime_row_from_instance_id(instance_id)
            requested_mem = float(row["meta_mem_total_gb"])
            data = list()
            data.append([(e / requested_mem) * 100 for e in by_id])
            averages = averages.append(pd.Series([sum(data[0]) / len(by_id)]),
                                       ignore_index=True)
        # averages.plot_bokeh(kind="scatter", title=f"Average Memory for {task}",
        #                    xlabel="Shard #", ylabel="Memory % Usage")
        # return averages.plot.line(title=f"Average Memory for {task}",
        #                          xlabel="Shard #", ylabel="Memory % Usage").get_figure()
        return averages

    def avg_cpu_across_shards(self, task: str):
        """
        Plots average CPU use for each shard of a given task
        :param task: task name
        :return: plot
        """
        # pandas_bokeh.output_file(f"avg_cpu_across_shards_{task}")
        id_list = self.get_shard_instance_id_list(task)
        averages = pd.Series([], dtype="float64")
        for instance_id in id_list:
            by_id = self.metrics[self.metrics["metrics_instance_id"] == instance_id]\
                        .sort_values("metrics_timestamp").loc[:, "metrics_cpu_used_percent"]
            data = list()
            data.append([sum(i) / len(i) for i in by_id])
            try:
                averages = averages.append(pd.Series([sum(data[0]) / len(data[0])]),
                                           ignore_index=True)  # 3870777233638278011 doesn't exist in metrics table
            except ZeroDivisionError:
                print(f"the instance id: {instance_id} in meta_runtime does not exist in metrics")
        # averages.plot_bokeh(kind="scatter", title=f"Average CPU for {task}", xlabel="Shard #", ylabel="CPU % Usage")
        # return averages.plot\
        # .line(title=f"Average CPU for {task}", xlabel="Shard #", ylabel="CPU % Usage").get_figure()
        return averages

    def avg_disk_across_shards(self, task: str):
        """
        Plots average disk use for each shard of a given task
        :param task: task name
        :return: plot
        """
        # pandas_bokeh.output_file(f"avg_disk_across_shards_{task}")
        id_list = self.get_shard_instance_id_list(task)
        averages = pd.Series([], dtype="float64")
        faster = self.metrics.loc[:, "metrics_instance_id"].astype('str').tolist()
        for instance_id in id_list:

            if str(instance_id) not in faster:  # 3870777233638278011 doesn't exist in metrics table
                print(f"the instance id: {instance_id} in meta_runtime does not exist in metrics")
                continue

            by_id = self.metrics[self.metrics["metrics_instance_id"] == instance_id]\
                        .sort_values("metrics_timestamp").loc[:, "metrics_disk_used_gb"]
            row = self.get_runtime_row_from_instance_id(instance_id)
            requested_disk = float(self.__get_first_val(row, 'meta_disk_total_gb'))
            data = list()
            data.append([(float(e[0]) / requested_disk) * 100 for e in by_id])
            #try:
            averages = averages.append(pd.Series([sum(data[0]) / len(by_id)]),
                                       ignore_index=True)  # 3870777233638278011 doesn't exist in metrics table
            #except ZeroDivisionError:
        # averages.plot_bokeh(kind="scatter", title=f"Average Disk for {task}", xlabel="Shard #", ylabel="Disk % Usage")
        # return averages.plot\
        # .line(title=f"Average Disk for {task}", xlabel="Shard #", ylabel="Disk % Usage").get_figure()
        return averages

    def avg_all_across_shards(self, task: str):
        """
        Creates a dataset of the combined data for disk space, memory, and cpu that can be easily made into bokeh plot
        :param task: task name
        :return: dataset
        """
        disk = self.avg_disk_across_shards(task)
        mem = self.avg_mem_across_shards(task)
        cpu = self.avg_cpu_across_shards(task)
        input_ = [self.task_input_size(task, i) for i in range(len(disk))]
        typedisk = ["disk %" for _ in range(len(disk))]
        typemem = ["mem %" for _ in range(len(mem))]
        typecpu = ["cpu %" for _ in range(len(cpu))]
        typeinput = ["input GB" for _ in range(len(disk))]
        x = [i % len(disk) for i in range(len(disk)*4)]
        # display(type(disk))
        # display(type(mem))
        # display(type(cpu))
        # display(type(pd.Series(data=input_)))
        y = disk.append(mem).append(cpu).append(pd.Series(data=input_))
        cat = typedisk + typemem + typecpu + typeinput
        combine = zip(y, cat, x)
        dataset = pd.DataFrame(data=tuple(combine), columns=["Percent Usage", "types", "Shard #"])
        # pandas_bokeh.output_file(f"avg_all_across_shards_{task}home.html")
        # scatter = dataset.plot_bokeh.scatter(x="Shard #", y="Percent Usage", line_alpha=0, fill_alpha=.6,
        #                                       category="types", title=Title(text=f"Combined plot for {task}",
        #                                       align="center"))
        return dataset

    def total_mem_usage(self) -> None:
        """
        Plots global memory usage in the order it's presented in the metrics table
        :return: plot
        """
        fix = self.metrics.loc[:, "metrics_mem_used_gb"].sort_values(by=["metrics_timestamp"]).reset_index(drop=True)
        fix.plot.line()  # x axis only goes to about half length?

    def total_cpu_usage(self) -> None:  # when there are multiple CPUs per shard take average
        """
        Plots global CPU usage in the order it's presented in the metrics table
        :return: plot
        """
        fix = self.metrics.loc[:, "metrics_cpu_used_percent"].reset_index(drop=True)
        fix = pd.Series([sum(fix[i]) / len(fix[i]) for i in range(len(fix))])  # very slow
        fix.plot.line()

    def clean_metrics(self, returned: bool) -> pd.DataFrame:
        """
        Removes failed shards from metrics
        :param returned: whether to return cleaned metrics or failed metrics (true for cleaned false for failed)
        :return: cleaned metrics
        """

        list_of_failures = [self.metrics[self.metrics["metrics_instance_id"] == i] for i in self.get_failed_shards_id()]
        # display(list_of_failures)
        index_list = list()
        for i in range(len(list_of_failures)):
            index_list += list(list_of_failures[i].index.values)  # apparently bad programming practice to do this
        if returned:
            return self.metrics.drop(index_list, axis=0)
        else:
            return self.metrics.loc[index_list, :]

    def get_failed_shards_id(self) -> list:
        """
        Gets list of instance IDs that failed
        :return: list of instance IDs
        """
        return self.meta_runtime[self.meta_runtime["meta_execution_status"] != "Done"]\
            .loc[:, "runtime_instance_id"].tolist()

    def shardless_task_table(self) -> pd.DataFrame:
        """
        Creates table of shardless tasks + available metrics info
        :return: table
        """
        tasks = self.num_task_shards()
        [tasks.pop(i) for i in self.num_task_shards().keys() if tasks[i] != 0]
        tasks = list(tasks.keys())
        ids = [self.get_instance_id(task, -1) for task in tasks]
        meta_disk = [self.meta_runtime[self.meta_runtime["runtime_instance_id"] == id_]["meta_disk_total_gb"]
                         .iloc[0] for id_ in ids]
        meta_mem = [self.meta_runtime[self.meta_runtime["runtime_instance_id"] == id_]
                    .loc[:, "meta_mem_total_gb"].iloc[0] for id_ in ids]
        avg_disk = [sum(map(self.nested_fix, self.metrics[self.metrics["metrics_instance_id"] == id_]
                    ["metrics_disk_used_gb"]))/len(self.metrics[self.metrics["metrics_instance_id"] == id_]
                    ["metrics_disk_used_gb"]) for id_ in ids]
        avg_mem = [sum(self.metrics[self.metrics["metrics_instance_id"] == id_]
                   ["metrics_mem_used_gb"])/len(self.metrics[self.metrics["metrics_instance_id"] == id_]
                   ["metrics_mem_used_gb"]) for id_ in ids]
        avg_cpu_percent = [sum(map(self.nested_fix, self.metrics[self.metrics["metrics_instance_id"] == id_]
                           ["metrics_cpu_used_percent"]))/(len(self.metrics[self.metrics["metrics_instance_id"] == id_]
                                                               ["metrics_cpu_used_percent"])) for id_ in ids]
        avg_disk = list(map(lambda x: x*100, map(truediv, avg_disk, map(self.nested_fix, meta_disk))))
        avg_mem = list(map(lambda x: x*100, map(truediv, avg_mem, meta_mem)))
        inputs = [self.task_input_size(task, -1) for task in tasks]  # pycharm highlights "inputs" strangely
        time = [len(self.metrics[self.metrics["metrics_instance_id"] == id_]) for id_ in ids]
        table = pd.DataFrame(data={"Disk %": avg_disk, "Memory %": avg_mem, "CPU %": avg_cpu_percent,
                                   "Input Size GB": inputs, "Runtime Sec": time},
                             index=tasks)
        # display(table)
        return table

    @staticmethod
    def nested_fix(x):
        """
        Utility function to convert nested list/series data types due to pandas being difficult
        :param x:
        :return:
        """
        return float(sum(x)/len(x))

    def plot_aggregate(self, task: str, shard: int, type_="svg") -> None:
        """
        Saves all matplotlib violin plots for 1 shard to all_plots directory
        :param task: task name
        :param shard: shard #
        :param type_: pdf or svg
        :return:
        """
        a = self.cpu_use_violin(self.get_instance_id(task, shard))
        # d = avg_mem_across_shards(meta_runtime, metrics, task)
        b = self.mem_use_violin(self.get_instance_id(task, shard))
        # e = avg_cpu_across_shards(meta_runtime, metrics, task)
        c = self.disk_use_violin(self.get_instance_id(task, shard))
        # f = avg_disk_across_shards(meta_runtime, metrics, task)

        # d = timed_plot_disk(meta_runtime, metrics, get_instance_id(meta_runtime, task, shard))
        # e = timed_plot_mem(meta_runtime, metrics, get_instance_id(meta_runtime, task, shard))
        # f = timed_plot_cpu(meta_runtime, metrics, get_instance_id(meta_runtime, task, shard))
        if type_ == "pdf":
            pp = PdfPages(f'all_plots/allPlots{task}.pdf')
            pp.savefig(a, bbox_inches='tight')
            pp.savefig(b, bbox_inches='tight')
            pp.savefig(c, bbox_inches='tight')
            pp.close()
        elif type_ == "svg":
            a.savefig(f'all_plots/{self.name}_allPlots{task}.pdf', bbox_inches='tight')
            b.savefig(f'all_plots/{self.name}_allPlots{task}.pdf', bbox_inches='tight')
            c.savefig(f'all_plots/{self.name}_allPlots{task}.pdf', bbox_inches='tight')
        else:
            raise Exception("File type not supported")

    def save_bokeh(self, task: str) -> None:
        """
        Saves all bokeh scatter plots to gridded_plots directory
        :param task: task name
        :return: plot
        """
        p_all = self.avg_all_across_shards(task).plot_bokeh\
            .scatter(x="Shard #", y="Percent Usage", line_alpha=0, fill_alpha=.6, category="types",
                     title=Title(text=f"Combined plot for {task}", align="center"), show_figure=False)
        p_disk = self.avg_disk_across_shards(task)\
            .plot_bokeh(kind="scatter", title=f"Average Disk for {task}", xlabel="Shard #", ylabel="Disk % Usage",
                        line_alpha=0, fill_alpha=.6, show_figure=False)
        p_cpu = self.avg_cpu_across_shards(task)\
            .plot_bokeh(kind="scatter", title=f"Average CPU for {task}", xlabel="Shard #", ylabel="CPU % Usage",
                        line_alpha=0, fill_alpha=.6, show_figure=False)
        p_mem = self.avg_mem_across_shards(task)\
            .plot_bokeh(kind="scatter", title=f"Average Memory for {task}", xlabel="Shard #", ylabel="Memory % Usage",
                        line_alpha=0, fill_alpha=.6, show_figure=False)
        saved = pandas_bokeh.plot_grid([[p_all, p_disk], [p_cpu, p_mem]],
                                       plot_width=450, show_plot=False, return_html=True)
        html = r"""
        <script type="text/x-mathjax-config">
          MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});
        </script>
        <script type="text/javascript"
          src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
        </script>""" + saved
        with open(f'gridded_plots/{self.name}_gridded{task}.html', 'w') as writer:
            writer.write(html)

    def timed_plot_disk(self, instance_id: int):
        """
        Plots disk usage for a process by time
        :param instance_id: instance id of process
        :return: plot
        """
        filtered = self.metrics[self.metrics["metrics_instance_id"] == instance_id]\
            .sort_values(by="metrics_timestamp", ignore_index=True)
        try:
            return pd.DataFrame(map(self.nested_fix, filtered.loc[:, "metrics_disk_used_gb"])).plot_bokeh(kind="line")
        except IndexError:
            print("instance ID not found in metrics")

    def timed_plot_mem(self, instance_id: int):
        """
        Plots memory usage for a process by time
        :param instance_id: instance id of process
        :return: plot
        """
        filtered = self.metrics[self.metrics["metrics_instance_id"] == instance_id]\
            .sort_values(by="metrics_timestamp", ignore_index=True)
        try:
            return pd.DataFrame(filtered.loc[:, "metrics_mem_used_gb"]).plot_bokeh(kind="line")
        except IndexError:
            print("instance ID not found in metrics")

    def timed_plot_cpu(self, instance_id: int):
        """
        Plots cpu usage for a process by time
        :param instance_id: instance id of process
        :return: plot
        """
        filtered = self.metrics[self.metrics["metrics_instance_id"] == instance_id]\
            .sort_values(by="metrics_timestamp", ignore_index=True)
        try:
            return pd.DataFrame(map(self.nested_fix, filtered.loc[:, "metrics_cpu_used_percent"]))\
                .plot_bokeh(kind="line")
        except IndexError:
            print("instance ID not found in metrics")


# Functions for combined instance: measuring variability, predictions, etc.
def variation(*args: ResourceUsage):
    frame_max_cpu = list()
    frame_max_mem = list()
    frame_max_disk = list()
    for resource_usage in args:
        frame_max_cpu.append(resource_usage.avg_max_cpu_usage)
        frame_max_cpu.append(resource_usage.max_mem_usage)
        frame_max_disk.append(resource_usage.max_disk_usage)
    total_max_cpu = pd.concat(frame_max_cpu)
    total_max_mem = pd.concat(frame_max_mem)
    total_max_disk = pd.concat(frame_max_disk)
