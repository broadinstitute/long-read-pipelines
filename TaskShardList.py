import pandas as pd

def df_task_shards_dict(df: pd.DataFrame) -> dict: 
    
    """
    supposed to return dict with task names as key and max shard index # as value
    """

    series_tasks = df.loc[:,'runtime_task_call_name']
    series_shards = df.loc[:,'runtime_shard']
    lowestvalues = [-1] * series_tasks.size
    task_shards_dict = dict(zip(series_tasks, lowestvalues)) #make dict of task:-1 pairs
    list_all = list(zip(series_tasks, series_shards))
    for key in task_shards_dict: 
        for element in list_all:
            if element[0] == key and element[1] > task_shards_dict[key]:
                task_shards_dict[key] = element[1]
    return task_shards_dict
