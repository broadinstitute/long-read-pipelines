#!/bin/bash

# Fast script to find orphan tasks and workflows
# Uses direct pattern matching instead of calling trace script for each item

set -euo pipefail

TSV_FILE="ALL_CALLABLES_IN_TASK_FOLDER.tsv"
WDL_DIR="./wdl"
OUTPUT_FILE="orphan_callables.tsv"
TEMP_CALLS_FILE="/tmp/all_calls.txt"

# Check if required files exist
if [[ ! -f "$TSV_FILE" ]]; then
    echo "Error: TSV file '$TSV_FILE' not found" >&2
    exit 1
fi

if [[ ! -d "$WDL_DIR" ]]; then
    echo "Error: WDL directory '$WDL_DIR' not found" >&2
    exit 1
fi

# Function to extract all call statements from WDL files
extract_all_calls() {
    echo "Extracting all call statements from WDL files..."
    
    # Find all call statements in WDL files
    # This captures patterns like:
    # - call TaskName
    # - call Namespace.TaskName
    # - call TaskName as Alias
    find "$WDL_DIR" -name "*.wdl" -exec grep -h "call [0-9A-Za-z_.-]\+" {} \; 2>/dev/null | \
        sed -E 's/.*call[[:space:]]+([0-9A-Za-z_.-]+).*$/\1/' | \
        sed -E 's/^[^.]*\.([^.]+)$/\1/' | \
        sort | uniq > "$TEMP_CALLS_FILE"
    
    echo "Found $(wc -l < "$TEMP_CALLS_FILE") unique callable names in call statements"
}

# Function to find orphans
find_orphans() {
    echo "Finding orphan tasks and workflows..."
    
    # Create output file with header
    echo -e "name\tclass\tfile_path\tline_number" > "$OUTPUT_FILE"
    
    local total_count=0
    local orphan_count=0
    local processed=0
    
    # Count total items (excluding header)
    total_count=$(tail -n +2 "$TSV_FILE" | wc -l)
    
    # Process each callable
    while IFS=$'\t' read -r name class file_path line_number; do
        ((processed++))
        
        # Skip header line
        if [[ "$name" == "name" ]]; then
            continue
        fi
        
        # Show progress
        if ((processed % 50 == 0)); then
            echo "Processed $processed/$total_count callables..."
        fi
        
        # Check if this callable name appears in any call statement
        if ! grep -q "^${name}$" "$TEMP_CALLS_FILE"; then
            echo -e "${name}\t${class}\t${file_path}\t${line_number}" >> "$OUTPUT_FILE"
            ((orphan_count++))
        fi
        
    done < "$TSV_FILE"
    
    echo ""
    echo "=== ORPHAN DETECTION SUMMARY ==="
    echo "Total callables processed: $total_count"
    echo "Orphan callables found: $orphan_count"
    echo "Results saved to: $OUTPUT_FILE"
    echo ""
    
    # Show orphan summary
    if [[ $orphan_count -gt 0 ]]; then
        echo "=== ORPHAN CALLABLES ==="
        echo "Format: Name | Class | File"
        echo "----------------------------"
        tail -n +2 "$OUTPUT_FILE" | while IFS=$'\t' read -r name class file_path line_number; do
            echo "$name | $class | $file_path:$line_number"
        done
        echo ""
        echo "Summary by class:"
        echo "Tasks: $(tail -n +2 "$OUTPUT_FILE" | grep -c "task")"
        echo "Workflows: $(tail -n +2 "$OUTPUT_FILE" | grep -c "workflow")"
    else
        echo "No orphan callables found!"
    fi
    
    # Cleanup
    rm -f "$TEMP_CALLS_FILE"
}

# Main execution
main() {
    echo "=== FAST ORPHAN CALLABLE FINDER ==="
    echo "TSV file: $TSV_FILE"
    echo "WDL directory: $WDL_DIR"
    echo "Output file: $OUTPUT_FILE"
    echo ""
    
    extract_all_calls
    find_orphans
}

main "$@"