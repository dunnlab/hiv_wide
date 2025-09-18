#!/bin/bash

# Define the cutoff date (e.g., files older than January 1, 2024)
CUTOFF_DATE="2025-07-14"

# Create a temporary reference file with the specified date
echo touch -t $(date d "$CUTOFF_DATE" +%Y%m%d%H%M) /tmp/reference_timestamp

FILENAME_PATTERN="*mask00*" 

# Find files older than the reference timestamp and execute a command
# -type f: Only consider regular files
# ! -newer /tmp/reference_timestamp: Find files NOT newer than the reference (i.e., older or same age)
# -exec your_command {} \;: Execute 'your_command' for each found file.
#                         Replace 'your_command' with the actual command you want to run.
#                         For example, 'rm' to delete, 'mv' to move, or 'echo' to list.
find /Users/aguang/CORE/kantorlab/hiv_wide/results/trees/samples/10576383 -type f -name "$FILENAME_PATTERN" ! -newer /tmp/reference_timestamp -print0 | while IFS= read -r -d $'\0' FILE; do
    echo "Processing: {}"
done

# Clean up the temporary reference file
rm /tmp/reference_timestamp
