#!/usr/bin/env python3

import json
import pprint
from datetime import datetime

class GetWorkflowIdTree:

    ROOT_ID = 'ROOT_ID'

    def __init__(self, slim_metadata_file):
        self.subworkflows_and_their_levels = dict()
        with open(slim_metadata_file, 'r') as json_file:
            self.slim_metadata = json.load(json_file)

    @staticmethod
    def _get_time_elapsed(dictionary, DATE_TIME_FMT = "%Y-%m-%dT%H:%M:%S.%f%z"):
        start_time = datetime.strptime(dictionary['start'], DATE_TIME_FMT)
        end_time = datetime.strptime(dictionary['end'], DATE_TIME_FMT)
        tdelta = end_time - start_time
        msg = f"Started on {start_time} and finished on {end_time} (totaling {tdelta})."
        return(msg)

    @staticmethod
    def _deal_with_subworkflow_or_task(this_level_metadata,
                                       parent_level,
                                       parent_level_wid,
                                       parent_level_dict,
                                       debug = False):

        local_dict = dict()
        current_level = 1 + parent_level

        for unit_name, metadata in this_level_metadata.items():

            # metadata is a list:
            #   if the unit is a preemptible task, several attempts may be made,
            #     then the elements in 'metadata' will be all the attempts actually made for the task
            #   otherwise it's a subworkflow (which itself may contain its own subworkflows)
            #     and we assume that subworkflows are not attempted multiple times (hence 1 == len)
            for struct in metadata:

                time_stamps = GetWorkflowIdTree._get_time_elapsed(struct)
                this_is_simple_task = not 'subWorkflowMetadata' in list(struct.keys())

                if ( this_is_simple_task ):
                    if (debug):
                        task_name = unit_name; # just a rename to better reflect what it is
                        attempt = struct
                        numb_attempt = attempt['attempt']
                        backend_status = attempt['backendStatus']
                        print(f"For simple task {task_name} on attempt {numb_attempt}, its backend status is {backend_status} (timewise: {time_stamps}).")
                else:
                    subworkflow_name = unit_name; # just a rename to better reflect what it is
                    subworkflow = struct

                    subworkflow_metadata = subworkflow['subWorkflowMetadata']
                    subworkflow_id = subworkflow_metadata['id']

                    local_dict[ subworkflow_name ] = {"wid": subworkflow_id, \
                                                      "level": current_level, \
                                                      "parent": parent_level_wid}

                    has_further_subworkflow = 'subWorkflowMetadata' in list(subworkflow_metadata['calls'].keys())
                    if (has_further_subworkflow):
                        child_dict = GetWorkflowIdTree.__deal_with_subworkflow_or_task(subworkflow_metadata['calls'], current_level, subworkflow_id, local_dict, debug)
                        merged_dict = {**local_dict, **child_dict}
                        local_dict = merged_dict
                    elif (debug):
                        print(f"For subworkflow {subworkflow_name}, its id is {subworkflow_id}, its status is {subworkflow['executionStatus']} (timewise: {time_stamps}).")

        return(local_dict)

    def parse(self, debug = False):
        root_level_id = self.slim_metadata['id']

        self.subworkflows_and_their_levels = GetWorkflowIdTree._deal_with_subworkflow_or_task(self.slim_metadata['calls'],
                                                                                              0, root_level_id,
                                                                                              dict(), debug)

        self.subworkflows_and_their_levels[GetWorkflowIdTree.ROOT_ID] = root_level_id

        pp = pprint.PrettyPrinter(indent=2)
        pp.pprint(self.subworkflows_and_their_levels)


############################################################

if __name__ == "__main__":
    import sys

    if 1 == len(sys.argv):
        import os
        my_name = os.path.basename(__file__)
        print(f"Usage: ./{my_name} <slim_metadata_json_file>")
    else:
        my_parser = GetWorkflowIdTree(sys.argv[1])
        my_parser.parse()
