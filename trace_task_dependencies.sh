#!/bin/bash

# Script to trace WDL task dependencies
# Usage: ./trace_task_dependencies.sh <task_name>
# Example: ./trace_task_dependencies.sh GetTodayDate

set -euo pipefail

TASK_NAME="$1"
WDL_DIR="./wdl"
VISITED=()
CHAIN=()

# Function to check if a value exists in an array
contains_element() {
    local element="$1"
    shift
    local array=("$@")
    for item in "${array[@]}"; do
        [[ "$item" == "$element" ]] && return 0
    done
    return 1
}

# Function to add visited file to avoid cycles
add_visited() {
    local file="$1"
    VISITED+=("$file")
}

# Function to check if file was already visited
is_visited() {
    local file="$1"
    if [[ ${#VISITED[@]} -eq 0 ]]; then
        return 1
    fi
    contains_element "$file" "${VISITED[@]}"
}

# Function to find the file containing a task definition
find_task_definition() {
    local task="$1"
    # Search for task definition pattern: "^task TaskName {"
    find "${WDL_DIR}" -name "*.wdl" -exec grep -l "^task ${task}\s*{" {} \; 2>/dev/null || true
}

# Function to find workflows/tasks that call a specific task
find_callers() {
    local task="$1"
    local calling_files=()
    
    # Search for "call TaskName" or "call Namespace.TaskName" patterns
    while IFS= read -r -d '' file; do
        if grep -q "call.*${task}" "$file" 2>/dev/null; then
            calling_files+=("$file")
        fi
    done < <(find "${WDL_DIR}" -name "*.wdl" -print0)
    
    # Only print if array has elements to avoid "unbound variable" error
    if [[ ${#calling_files[@]} -gt 0 ]]; then
        printf '%s\n' "${calling_files[@]}"
    fi
}

# Function to get workflow/task name from file
get_workflow_or_task_name() {
    local file="$1"
    # Look for workflow or task definition at the top level
    local name=""
    name=$(grep -E "^(workflow|task) [0-9A-Za-z_-]+" "$file" | head -1 | awk '{print $2}' | sed 's/{.*$//')
    echo "$name"
}

# Function to trace dependencies recursively
trace_dependencies() {
    local current_task="$1"
    local depth="${2:-0}"
    local indent=$(printf "%*s" $((depth * 2)) "")
    
    echo "${indent}Tracing callers of task: $current_task"
    
    # Find the task definition file
    local task_file
    task_file=$(find_task_definition "$current_task")
    
    if [[ -z "$task_file" ]]; then
        echo "${indent}  Task '$current_task' definition not found"
        return
    fi
    
    echo "${indent}  Task defined in: $task_file"
    
    # Find all files that call this task
    local callers=()
    while IFS= read -r line; do
        [[ -n "$line" ]] && callers+=("$line")
    done < <(find_callers "$current_task")
    
    if [[ ${#callers[@]} -eq 0 ]]; then
        echo "${indent}  No callers found for task '$current_task'"
        return
    fi
    
    echo "${indent}  Found ${#callers[@]} caller(s):"
    
    for caller_file in "${callers[@]}"; do
        # Skip if already visited to avoid cycles
        if is_visited "$caller_file"; then
            echo "${indent}    $caller_file (already visited, skipping to avoid cycle)"
            continue
        fi
        
        add_visited "$caller_file"
        
        local caller_name
        caller_name=$(get_workflow_or_task_name "$caller_file")
        
        echo "${indent}    $caller_file"
        echo "${indent}      Contains workflow/task: $caller_name"
        
        # Add to current chain
        CHAIN+=("$caller_name -> $caller_file")
        
        # Recursively trace callers of this workflow/task
        if [[ -n "$caller_name" && "$caller_name" != "$current_task" ]]; then
            trace_dependencies "$caller_name" $((depth + 1))
        fi
    done
}

# Function to print the full dependency chain
print_chain() {
    echo ""
    echo "=== DEPENDENCY CHAIN SUMMARY ==="
    if [[ ${#CHAIN[@]} -eq 0 ]]; then
        echo "No dependency chain found for task '$TASK_NAME'"
    else
        echo "Dependency chain for task '$TASK_NAME':"
        for item in "${CHAIN[@]}"; do
            echo "  $item"
        done
    fi
}

# Main execution
main() {
    if [[ $# -eq 0 ]]; then
        echo "Usage: $0 <task_name>"
        echo "Example: $0 GetTodayDate"
        exit 0
    fi
    
    echo "=================================="
    echo "=== WDL TASK DEPENDENCY TRACER ==="
    echo "=================================="
    echo -e "Searching for workflows that call task: \n\t$TASK_NAME"
    echo "WDL directory: $WDL_DIR"
    echo ""
    
    # Check if WDL directory exists
    if [[ ! -d "$WDL_DIR" ]]; then
        echo "Error: WDL directory '$WDL_DIR' not found"
        exit 1
    fi
    
    # Start tracing
    trace_dependencies "$TASK_NAME"
    
    # Print summary
    print_chain
}

main "$@"
