#!/bin/bash

# Script to list all workflow and task definitions in wdl/tasks folder
# Outputs a 4-column TSV: name, class, file_path, line_number

set -euo pipefail

WDL_TASKS_DIR="wdl/tasks"

# Function to extract workflows and tasks from a file
extract_definitions() {
    local file="$1"
    
    # Use grep with line numbers to find workflow and task definitions
    # Pattern matches: "workflow Name {" or "task Name {"
    grep -n "^[[:space:]]*\(workflow\|task\)[[:space:]]\+[0-9A-Za-z_-]\+[[:space:]]*{" "$file" 2>/dev/null | while IFS=: read -r line_num match; do
        # Extract the type (workflow or task) and name
        local type name
        type=$(echo "$match" | sed -E 's/^[[:space:]]*(workflow|task)[[:space:]].*/\1/')
        name=$(echo "$match" | sed -E 's/^[[:space:]]*(workflow|task)[[:space:]]+([0-9A-Za-z_-]+)[[:space:]]*\{.*/\2/')
        
        # Output TSV line: name, class, file_path, line_number
        echo -e "${name}\t${type}\t${file}\t${line_num}"
    done
}

# Main execution
main() {
    local output_tsv="$1"
    rm -f "${output_tsv}"

    # Check if WDL tasks directory exists
    if [[ ! -d "$WDL_TASKS_DIR" ]]; then
        echo "Error: WDL tasks directory '$WDL_TASKS_DIR' not found" >&2
        exit 1
    fi
    
    # Print TSV header
    echo -e "name\tclass\tfile_path\tline_number" > "${output_tsv}"
    
    # Find all .wdl files in wdl/tasks and extract definitions
    find "$WDL_TASKS_DIR" -name "*.wdl" -type f | sort | while read -r file; do
        extract_definitions "$file" >> "${output_tsv}"
    done
}

main "$@" ALL_CALLABLES_IN_TASK_FOLDER.tsv