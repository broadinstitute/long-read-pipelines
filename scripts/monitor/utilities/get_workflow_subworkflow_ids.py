#!/usr/bin/env python3

import json
import sys
import os
import pprint
from datetime import datetime


class GetWorkflowIdTree:

    """
    Given slim-metadata JSON file, parse the root- and sub-workflow ID tree.
    """

    ROOT_ID = 'ROOT_ID'

    def __init__(self, slim_metadata_file):
        with open(slim_metadata_file, 'r') as json_file:
            self.slim_metadata = json.load(json_file)

    @staticmethod
    def _get_time_elapsed(task_metadata, dt_format="%Y-%m-%dT%H:%M:%S.%f%z"):
        """
        Mostly for debugging
        @param task_metadata:
        @param dt_format:
        @return:
        """
        start_time = datetime.strptime(task_metadata['start'], dt_format)
        end_time = datetime.strptime(task_metadata['end'], dt_format)
        delta = end_time - start_time
        msg = f"Started on {start_time} and finished on {end_time} (totaling {delta})."
        return msg

    @staticmethod
    def _deal_with_subworkflow_or_task(this_level_metadata,
                                       parent_level,
                                       parent_level_wid,
                                       parent_name=None,
                                       debug=False):

        """
        Major workpiece, recursively parse (DFS) tree.
        @param this_level_metadata:
        @param parent_level:
        @param parent_level_wid:
        @param debug:
        @return:
        """

        subworkflow_info_for_this_node = dict()
        current_level = 1 + parent_level

        pretty_parent_name = parent_name
        if 0 == parent_level and parent_name is None:
            pretty_parent_name = list(this_level_metadata.keys())[0].split(".")[0]

        for unit_name, metadata in this_level_metadata.items():

            # metadata is a list:
            #   if the unit is a preemptible task, several attempts may be made,
            #     then the elements in 'metadata' will be all the attempts actually made for the task
            #   otherwise it's a subworkflow (which itself may contain its own subworkflows)
            #     and we assume that
            #       1) simple tasks don't have their own sub-workflows
            #       2) (sub-)workflows are not attempted multiple times, only tasks are
            for struct in metadata:

                time_stamps = GetWorkflowIdTree._get_time_elapsed(struct)

                is_subworkflow = 'subWorkflowMetadata' in struct

                if is_subworkflow:
                    # just a rename to better reflect what it is
                    subworkflow_name = unit_name
                    subworkflow = struct

                    pretty_name = pretty_parent_name + "." + subworkflow_name.split(".")[-1]

                    subworkflow_metadata = subworkflow['subWorkflowMetadata']
                    subworkflow_id = subworkflow_metadata['id']

                    subworkflow_info_for_this_node[pretty_name] = {"wid": subworkflow_id,
                                                                   "level": current_level,
                                                                   "parent": parent_level_wid}
                    # DFS
                    child_subworkflow_info = GetWorkflowIdTree._deal_with_subworkflow_or_task(subworkflow_metadata['calls'],
                                                                                              current_level, subworkflow_id,
                                                                                              pretty_name, debug)
                    if 0 < len(child_subworkflow_info):
                        merged_dict = {**subworkflow_info_for_this_node, **child_subworkflow_info}
                        subworkflow_info_for_this_node = merged_dict
                    elif debug:
                        print(f"For subworkflow {subworkflow_name}, its id is {subworkflow_id}, "
                              f"  its status is {subworkflow['executionStatus']} "
                              f"  (time-wise: {time_stamps}).")
                elif debug:
                    task_name = unit_name  # just a rename to better reflect what it is
                    attempt = struct  # just a rename to better reflect what it is
                    numb_attempt = attempt['attempt']
                    backend_status = attempt['backendStatus']
                    print(f"For simple task {task_name} on attempt {numb_attempt}, "
                          f"  its backend status is {backend_status}"
                          f"  (time-wise: {time_stamps}).")

        return subworkflow_info_for_this_node

    def parse(self, debug=False):
        """
        Main interface
        @param debug:
        @return:
        """

        tree = dict()

        root_level_id = self.slim_metadata['id']
        tree[GetWorkflowIdTree.ROOT_ID] = root_level_id

        calls = self.slim_metadata['calls']
        child_tree = GetWorkflowIdTree._deal_with_subworkflow_or_task(calls, 0, root_level_id, debug=debug)

        tree = {**tree, **child_tree}

        pp = pprint.PrettyPrinter(indent=2)
        pp.pprint(tree)


############################################################

if __name__ == "__main__":

    if 1 == len(sys.argv):
        my_name = os.path.basename(__file__)
        print(f"Usage: ./{my_name} <slim_metadata_json_file>")
    else:
        my_parser = GetWorkflowIdTree(sys.argv[1])
        my_parser.parse()
