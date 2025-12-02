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

# Extract unique subject IDs from folder names
# Pattern: 3<subject_name>_<project_id>[_<scan_date>]
# Exclude folders with PET or KIM in the name
subject_ids=$(find "$SCAN_DIR" -maxdepth 1 -type d -name "3*_*" | \
    grep -v -E "(PET|KIM|_SI)" | \
    sed 's|.*/||' | \
    sed 's/_[0-9]\{8\}$//' | \
    sort -u)

# Count total subjects
total_subjects=$(echo "$subject_ids" | wc -l)

# Check if we have enough subjects
if [ "$total_subjects" -lt 50 ]; then
    echo "Warning: Only $total_subjects subjects found (less than 50 requested)" >&2
    selected_subjects="$subject_ids"
else
    # Randomly select 50 subjects
    selected_subjects=$(echo "$subject_ids" | shuf -n 50)
fi

# For each selected subject, find all their session folders
for subject in $selected_subjects; do
    find "$SCAN_DIR" -maxdepth 1 -type d -name "${subject}*" | \
        grep -v -E "(PET|KIM)" | \
        sed 's|.*/||'
done | sort