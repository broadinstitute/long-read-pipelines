def df_task_shards_dict(df: pd.Dataframe) -> dict: 
    """
    supposed to return dict with task names as key and max shard index # as value
    """
    series_tasks = df.loc[:,'runtime_task_call_name']
    series_shards = df.loc[:,'runtime_shard']
    lowestvalues = [-1] * series_tasks.size
    task_shards_dict = dict(zip(series_tasks, lowestvalues)) #make dict of task:-1 pairs
    list_all = list(zip(series_tasks, series_shards))
    for key in task_shards_dict: #select largest shard index for task
        for element in list_all:
            if element[0] == key and element[1] > task_shards_dict[key]:
                task_shards_dict[key] = element[1]
    return task_shards_dict


def task_disk_meta_runtime(input_df: pd.DataFrame, task: str) -> pd.DataFrame:
    """
    disk usage for a task by shards
    """
    columns_to_select = ["runtime_task_call_name", "runtime_shard", "runtime_disk_total_gb", "meta_disk_total_gb", "meta_execution_status"]
    renaming_scheme = {"runtime_task_call_name":"task","runtime_shard":"shard_idx"}
    df = input_df.loc[:, columns_to_select].rename(columns = renaming_scheme)
    df = df[df["meta_execution_status"] == "Done"]
    if task in df['task'].unique():
        return df[df["task"] == task].sort_values(by=['shard_idx'])
    else:
        raise KeyError("Task name does not exist")
        
        
def task_mem_meta_runtime(input_df: pd.DataFrame, task: str) -> pd.DataFrame:
    """
    mem usage for a task by shards
    """
    columns_to_select = ["runtime_task_call_name", "runtime_shard", "runtime_mem_total_gb", "meta_mem_total_gb", "meta_execution_status"]
    renaming_scheme = {"runtime_task_call_name":"task","runtime_shard":"shard"}
    df = ds_meta_runtime.loc[:,columns_to_select].rename(columns = renaming_scheme)
    df = df[df["meta_execution_status"] == "Done"]
    if task in df['task'].unique():
        return df[df["task"] == task].sort_values(by=['shard'])
    else:
        raise KeyError("Task name does not exist")


def plot_task_disk_usage(df: pd.DataFrame):
    """
    yields more value with tasks split into many shards
    will only plot the first disk in array, use in conjunction with task_disk_meta_runtime
    """
    x1 = df.loc[:,"shard"]
    x2 = df.loc[:,"shard"]
    y1 = df.loc[:,"runtime_disk_total_gb"].map(lambda arr: arr[0]).astype("float64", copy=True)
    y2 = df.loc[:,"meta_disk_total_gb"].map(lambda arr: arr[0]).astype("float64", copy=True)
    plt.plot(x1, y1)
    plt.plot(x2, y2)
    plt.legend(["Runtime Disk Usage(gb)", "Meta Disk Usage(gb)"])
    

def plot_task_mem_usage(df: pd.DataFrame):
    """
    yields more value with tasks split into many shards
    use in conjunction with task_mem_meta_runtime
    """
    x1 = df.loc[:,"shard"]
    x2 = df.loc[:,"shard"]
    y1 = df.loc[:,"runtime_mem_total_gb"].astype("float64", copy=True)
    y2 = df.loc[:,"meta_mem_total_gb"].astype("float64", copy=True)
    plt.plot(x1, y1)
    plt.plot(x2, y2)
    plt.legend(["Runtime Memory Usage(gb)", "Meta Memory Usage(gb)"])


def task_input_size(task_name_or_index, shard_index = None) -> float: 
    """
    summation of a task's input file sizes
    use:
    task_input_size(row index of task/shard), or task_input_size(task name, shard index)
    """
    if shard_index:
        tasks = ds_meta_runtime.loc[:,"runtime_task_call_name"] #find row index
        shards = ds_meta_runtime.loc[:,"runtime_shard"]
        for i in range(0, tasks.size):
            if tasks[i] == task_name_or_index and shards[i] == shard_index: 
                pos = i
    else:
        if isinstance(task_name_or_index, int):
            pos = task_name_or_index #if only task_name is supplied it's treated as the row index
        else:
            raise Exception(f"Need an integer row index but recieved {task_name_or_index}")
    input_ds = pd.DataFrame(data=ds_meta_runtime.loc[pos,"meta_inputs"])
    return input_ds[input_ds['type'].astype(str)=="file"].loc[:,"value"].astype(float).sum()
