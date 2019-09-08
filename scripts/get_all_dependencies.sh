#!/bin/bash

set -eu 

if [[ "$#" -ne 1 ]]; then
	my_scriptname=$(basename "${0}")
	echo "I get all the dependencies you need for running the master workflow:"
	echo "  bash ${my_scriptname} WORKFLOW_WDL_FILE"
	exit 0
fi

function get_dep() {
	ff=$(grep -F 'import' "${1}" | sed 's/import //' | sed -E 's/ as [A-Za-z]+//' | sed 's/"//g')
	echo "${ff}"
}
export -f get_dep

function traverse() {
	for dep in $(get_dep "${1}"); do
		echo " ${dep}"
		traverse "${dep}"
	done
}
export -f traverse

zip "${1}.dependencies.zip" $(traverse "${1}" | sort | uniq)
