#!/bin/bash

set -eu

# For generating the minimal input json for a WDL


wdl=$1

contents=$(womtool inputs "${wdl}" | grep -vF 'runtime_attr_override' | grep -vF 'optional' | grep -v "{" | grep -v "}" | sort);

formatted_contents=${contents%?}

echo -e "{${formatted_contents}\n}"
