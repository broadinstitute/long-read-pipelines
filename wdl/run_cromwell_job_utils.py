from datetime import datetime
import json
from pprint import pprint
import os
import subprocess
import sys
import time


USER_RESPONSE = None
def user_wants_to_continue():
    global USER_RESPONSE
    if USER_RESPONSE is None:
        USER_RESPONSE = input("Continue? [Y/n] ")
    if USER_RESPONSE and USER_RESPONSE.upper().startswith("Y"):
        return True

    return False


def execute(args, tool_name, data, wdl_path, workflow_options, input_json, data_root_path):

    # set up output dir
    output_dir = f"cromwell-executions/{tool_name}/{datetime.now().strftime('%Y%m%d_%H%M%S')}__{args.tool}"
    #if args.label:
    #    output_dir += f"__{args.label[:50].replace('_', ' ')}"
    os.system(f"mkdir -p {output_dir}")


    final_input_json = {}
    for key, value in input_json.items():
        if isinstance(value, str) and (os.path.exists(os.path.join(data_root_path, value)) or ("/" in value and "." in value) or ("GRCh37/" in value) or ("GRCh38/" in value)):
            # make file paths absolute
            value = os.path.join(data_root_path, value)

        key = f"{tool_name}.{key}"

        final_input_json[key] = value

        if (key.endswith("_bam") or key.endswith("_cram")) and (value.endswith(".bam") or value.endswith(".cram")):
            # add .bam.bai path for every .bam
            final_input_json[key.replace("_bam", "_bam_bai").replace("_cram", "_cram_crai")] = value.replace(".bam", ".bam.bai").replace(".cram", ".cram.crai")


    input_json_path = f"{output_dir}/inputs.json"
    with open(input_json_path, "w") as f:
        json.dump(final_input_json, f)

    print("Inputs: ")
    pprint(final_input_json)

    # validate input json
    if not args.skip_input_validation:
        validation_error = os.system(f"python3 validate_json_inputs.py {input_json_path}")
        if validation_error:
            sys.exit(f"Validation failed for {input_json_path}")

    # copy the input wdl
    original_wdl_path = wdl_path
    wdl_name = input_json.get("output_prefix", f"{args.tool}__{data}")
    wdl_path = f"{output_dir}/{wdl_name}.wdl"
    os.system(f"cp {original_wdl_path} {wdl_path}")

    # compute the command string
    if args.local:
        os.system(f"sed -i.bak 's/String[ ]input/File input/g' {wdl_path}")

        cromwell_jar_path = "/Users/weisburd/code/cromwell/cromwell-36.jar"
        command = f"java -jar {cromwell_jar_path} run {wdl_path} -i {input_json_path}"  # -o {workflow_options_path}

    elif args.remote:
        workflow_options_path = f"{output_dir}/workflow_options.json"
        final_workflow_options = dict(workflow_options)
        #if args.tool in ["samtools-index-bam"]:
        #    final_workflow_options["final_workflow_outputs_dir"] = os.path.dirname(final_input_json[f"{WORKFLOWS[args.tool]['name']}.input_bam"])

        with open(workflow_options_path, "w") as f:
            json.dump(final_workflow_options, f)
        command = f"cromshell submit {wdl_path} {input_json_path} {workflow_options_path}"

    # run the command
    print(f"Command: {command}")
    if not args.test:
        os.system(command)

        if args.remote and data != args.data[-1]:
            if not user_wants_to_continue():
                return

            while True:
                time.sleep(15)
                jobs_list = subprocess.check_output("cromshell list -c -u | tail -n 1", shell=True)
                jobs_list = jobs_list.decode('UTF-8')
                print(jobs_list)
                if "Running" in jobs_list:
                    break

                if "Failed" in jobs_list:
                    sys.exit("Job failed. Exiting..")

