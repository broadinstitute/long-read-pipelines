#!/bin/bash

# Script to find orphan tasks and workflows (those not called by anything)
# Uses the trace_task_dependencies.sh script and ALL_CALLABLES_IN_TASK_FOLDER.tsv

set -euo pipefail

TSV_FILE="ALL_CALLABLES_IN_TASK_FOLDER.tsv"
TRACE_SCRIPT="./trace_task_dependencies.sh"
OUTPUT_FILE="orphan_callables.tsv"

# Check if required files exist
if [[ ! -f "$TSV_FILE" ]]; then
    echo "Error: TSV file '$TSV_FILE' not found" >&2
    exit 1
fi

if [[ ! -x "$TRACE_SCRIPT" ]]; then
    echo "Error: Trace script '$TRACE_SCRIPT' not found or not executable" >&2
    exit 1
fi

# Function to check if a callable is an orphan
is_orphan() {
    local name="$1"
    local class="$2"
    local file_path="$3"
    
    # Run the trace script and capture output
    local trace_output
    trace_output=$("$TRACE_SCRIPT" "$name" 2>/dev/null || true)
    
    # Check if no callers were found
    if echo "$trace_output" | grep -q "No callers found for"; then
        return 0  # Is orphan
    elif echo "$trace_output" | grep -q "No dependency chain found for"; then
        return 0  # Is orphan
    else
        return 1  # Not orphan
    fi
}

# Function to find orphans
find_orphans() {
    echo "Finding orphan tasks and workflows..."
    echo ""
    
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
        
        echo "Processing ($processed/$total_count): $name ($class)"
        
        # Check if it's an orphan
        if is_orphan "$name" "$class" "$file_path"; then
            echo "  -> ORPHAN: $name"
            echo -e "${name}\t${class}\t${file_path}\t${line_number}" >> "$OUTPUT_FILE"
            ((orphan_count++))
        else
            echo "  -> has callers"
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
        echo "Name | Class | File"
        echo "---- | ----- | ----"
        tail -n +2 "$OUTPUT_FILE" | while IFS=$'\t' read -r name class file_path line_number; do
            echo "$name | $class | $file_path"
        done
    else
        echo "No orphan callables found!"
    fi
}

# Function to show only final summary (faster mode)
find_orphans_summary_only() {
    echo "Finding orphan tasks and workflows (summary only mode)..."
    echo ""
    
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
        
        # Show progress every 20 items
        if ((processed % 20 == 0)); then
            echo "Processed $processed/$total_count callables..."
        fi
        
        # Check if it's an orphan
        if is_orphan "$name" "$class" "$file_path"; then
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
        echo "Name | Class | File"
        echo "---- | ----- | ----"
        tail -n +2 "$OUTPUT_FILE" | while IFS=$'\t' read -r name class file_path line_number; do
            echo "$name | $class | $file_path"
        done
    else
        echo "No orphan callables found!"
    fi
}

# Main execution
main() {
    local mode="${1:-verbose}"
    
    echo "=== ORPHAN CALLABLE FINDER ==="
    echo "TSV file: $TSV_FILE"
    echo "Trace script: $TRACE_SCRIPT"
    echo "Output file: $OUTPUT_FILE"
    echo ""
    
    if [[ "$mode" == "summary" ]]; then
        find_orphans_summary_only
    else
        find_orphans
    fi
}

main "$@"