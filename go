#!/usr/local/bin/python3
import subprocess
import sys
import time
from os import walk

return_code = subprocess.call(f'java -jar bin/womtool-36.1.jar validate wdl/correct_and_align/correct_and_align.wdl', shell=True)
print(f'womtool validate correct_and_align.wdl return_code={return_code}')

if return_code != 0:
    print("Validation failed.")
else:
    dry_run = True if len(sys.argv) == 1 else False

    print(f'\ndry_run: {dry_run}')

    f = []
    for (dirpath, dirnames, filenames) in walk("data/"):
        f.extend(filenames)
        break

    for file in sorted(f, reverse=True):
        if file.endswith(".json"):
            cmd = f'cromshell submit wdl/correct_and_align/correct_and_align.wdl data/{file} resources/gcs_workflow_options.json'
            if dry_run:
                subprocess.call(f'echo {cmd}', shell=True)
            else:
                return_code = subprocess.call(f'{cmd}', shell=True)
                time.sleep(5)
