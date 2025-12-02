#!/bin/bash

# Check if directory argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <directory_path>" >&2
    exit 1
fi

SCAN_DIR="$1"

# Check if directory exists
if [ ! -d "$SCAN_DIR" ]; then
    echo "Error: Directory '$SCAN_DIR' does not exist" >&2
    exit 1
fi

# Extract unique subject names (without project ID or date)
# Pattern: 3<subject_name>_<project_id>[_<scan_date>]
# Exclude folders with PET or KIM in the name
subject_names=$(find "$SCAN_DIR" -maxdepth 1 -type d -name "3*_*" | \
    grep -v -E "(PET|KIM|_SI)" | \
    sed 's|.*/||' | \
    sed -E 's/^(3[^_]+)_.*/\1/' | \
    sort -u)

# Count total subjects
total_subjects=$(echo "$subject_names" | wc -l)

# Check if we have enough subjects
if [ "$total_subjects" -lt 50 ]; then
    echo "Warning: Only $total_subjects subjects found (less than 50 requested)" >&2
    selected_subjects="$subject_names"
else
    selected_subjects=$(echo "$subject_names")
fi

# For each selected subject, find their baseline scan (earliest session)
for subject in $selected_subjects; do
    # Find all folders for this subject, excluding PET and KIM
    folders=$(find "$SCAN_DIR" -maxdepth 1 -type d -name "${subject}_*" | \
        grep -v -E "(PET|KIM)" | \
        sed 's|.*/||')
    
    # Select the baseline scan (first alphabetically, which should be earliest)
    baseline=$(echo "$folders" | sort | head -n 1)
    
    if [ -n "$baseline" ]; then
        echo "$baseline"
    fi
done | sort