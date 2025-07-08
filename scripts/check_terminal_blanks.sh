#!/bin/bash
#
# check_terminal_blanks.sh
#
# Description:
# This script takes a directory path as input. For each file in that directory,
# it calculates the proportion of terminal blank characters for every sequence.
# It then computes the minimum, 1st quartile (Q1), median, and 3rd quartile (Q3)
# of these proportions for the entire file.
#
# The final output is a tab-separated value (TSV) file printed to standard
# output, with condensed filenames as columns and statistics as rows.
#
# Usage:
# ./check_terminal_blanks.sh /path/to/your/fasta_files > summary_stats.tsv
#

# --- Script Setup and Input Validation ---

# Check if a directory path was provided
if [ -z "$1" ]; then
    echo "Error: No directory path provided." >&2
    echo "Usage: $0 <path_to_directory>" >&2
    exit 1
fi

DIR_PATH="$1"

# Check if the provided path is actually a directory
if [ ! -d "$DIR_PATH" ]; then
    echo "Error: Path '$DIR_PATH' is not a valid directory." >&2
    exit 1
fi

# --- AWK Program Definitions ---

# AWK_PROPORTIONS: Calculates and prints the proportion for each sequence.
# It's the engine that feeds the statistical analysis.
# AWK_PROPORTIONS: Calculates and prints the proportion for each sequence.
# It's the engine that feeds the statistical analysis.
AWK_PROPORTIONS='
function process_sequence(sequence) {
    if (length(sequence) == 0) return;
    total_length = length(sequence);
    match(sequence, /[ACGTUacgtu]/);
    first_valid_base = RSTART;
    last_valid_base = 0;
    for (i = total_length; i >= 1; i--) {
    	if (substr(sequence, i, 1) ~ /[ACGTUacgtu]/) {
	   last_valid_base = i;
	    break;
	}
    }
    leading_blanks = 0;
    trailing_blanks = 0;
    if (first_valid_base == 0) {
       leading_blanks = total_length;
    } else {
      leading_blanks = first_valid_base - 1;
      trailing_blanks = total_length - last_valid_base;
    }
    total_terminal_blanks = leading_blanks + trailing_blanks;
    proportion = (total_length > 0) ? (total_terminal_blanks / total_length) : 0;
    printf "%.4f\n", proportion;

    # Increment counter if proportion is >= 0.8860
    if (proportion >= 0.8860) {
        high_blank_sequences++;
    }
    total_sequences++;
}
/^>/ { if (sequence != "") process_sequence(sequence); sequence=""; next; }
{ sequence = sequence $0; }
END { 
    if (sequence != "") process_sequence(sequence); 
    
    # Print the proportion of sequences with high blank characters
    if (total_sequences > 0) {
        printf "HIGH_BLANK_PROPORTION:%.4f\n", high_blank_sequences / total_sequences;
    } else {
        printf "HIGH_BLANK_PROPORTION:0.0000\n";
    }
}
'

# AWK_STATS: Reads a sorted list of numbers and calculates Min, Q1, Median, and Q3.
# The function is now defined at the top level for maximum compatibility.
AWK_STATS='
# This function calculates the median of a subarray.
# The extra variables (len, mid_idx, mid_val) are declared as parameters,
# which makes them local to the function call in awk.
function get_median(a, start, end,    len, mid_idx, mid_val) {
    len = end - start + 1
    if (len <= 0) {
        return 0
    }
    if (len % 2 == 1) {
        mid_idx = start + int((len - 1) / 2)
        mid_val = a[mid_idx]
    } else {
        mid_idx = start + int(len / 2) - 1
        mid_val = (a[mid_idx] + a[mid_idx + 1]) / 2
    }
    return mid_val
}

# Main processing block to read all numbers into an array.
{ a[NR] = $1 }

# END block to perform calculations after reading all input.
END {
    n = NR
    if (n == 0) {
        # min, q1, median, q3
        print "0.0000\t0.0000\t0.0000\t0.0000"
        exit
    }

    # Since the input from `sort -n` is sorted, the first element is the minimum.
    min = a[1]

    median = get_median(a, 1, n)

    # Determine the indices for the lower and upper halves
    if (n % 2 == 1) {
        lower_half_end = int((n - 1) / 2)
        upper_half_start = int((n + 3) / 2)
    } else {
        lower_half_end = int(n / 2)
        upper_half_start = int(n / 2) + 1
    }

    q1 = get_median(a, 1, lower_half_end)
    q3 = get_median(a, upper_half_start, n)

    # Print in order: min, q1, median, q3
    printf "%.4f\t%.4f\t%.4f\t%.4f\n", min, q1, median, q3
}
'

# --- Main Processing Loop ---

# Create a temporary file to hold the intermediate results.
# This is safer than relying on complex bash arrays.
TMP_RESULTS_FILE=$(mktemp)
TMP_HIGH_BLANK_PROPORTIONS_FILE=$(mktemp) # New temporary file for high blank proportions

# Ensure the temporary files are removed when the script exits, even on error.
trap 'rm -f "$TMP_RESULTS_FILE" "$TMP_HIGH_BLANK_PROPORTIONS_FILE"' EXIT

# Loop over all files in the provided directory
for FILE in "$DIR_PATH"/*; do
    # Check if the item is a regular file
    if [ -f "$FILE" ]; then
	# Extract the condensed filename (e.g., "mask001" from "mask001.fasta")
	# This regex captures "mask" followed by digits.
        FILENAME=$(basename "$FILE" | sed -E 's/(mask[0-9]+).*/\1/')

	# Calculate stats for the file
	# 1. Get all proportions using AWK_PROPORTIONS. This now also prints the HIGH_BLANK_PROPORTION.
        #    We need to separate these outputs.
        RAW_OUTPUT=$(awk "$AWK_PROPORTIONS" "$FILE")
        
        # Extract the statistical proportions (numerical lines)
        PROPORTIONS=$(echo "$RAW_OUTPUT" | grep -v "HIGH_BLANK_PROPORTION")
        
        # Extract the high blank proportion line
        HIGH_BLANK_LINE=$(echo "$RAW_OUTPUT" | grep "HIGH_BLANK_PROPORTION")

        # Now, sort the proportions and pipe to AWK_STATS
        STATS=$(echo "$PROPORTIONS" | sort -n | awk "$AWK_STATS")

        # Store the results for this file in our temp file
        echo -e "$FILENAME\t$STATS" >> "$TMP_RESULTS_FILE"

        # Store the high blank proportion separately
        if [ -n "$HIGH_BLANK_LINE" ]; then
            # Extract just the numerical value for high blank proportion
            HIGH_BLANK_VALUE=$(echo "$HIGH_BLANK_LINE" | sed -E 's/HIGH_BLANK_PROPORTION:(.*)/\1/')
	    echo -e "$FILENAME\t$HIGH_BLANK_VALUE" >> "$TMP_HIGH_BLANK_PROPORTIONS_FILE"
        else
            echo -e "$FILENAME\t0.0000" >> "$TMP_HIGH_BLANK_PROPORTIONS_FILE"
        fi
   fi
done

# --- Final Output Formatting ---

# Check if any results were generated
if [ ! -s "$TMP_RESULTS_FILE" ]; then
    echo "No valid sequence files found or processed in '$DIR_PATH'." >&2
    exit 0
fi

# Transpose the results from the temporary file to the final format.
# Use cut to extract columns and tr to turn newlines into tabs.
# Use sed to remove the trailing tab that tr adds.
HEADERS=$(cut -f1 "$TMP_RESULTS_FILE" | tr '\n' '\t' | sed 's/\t$//')
MIN_VALS=$(cut -f2 "$TMP_RESULTS_FILE" | tr '\n' '\t' | sed 's/\t$//')
Q1_VALS=$(cut -f3 "$TMP_RESULTS_FILE" | tr '\n' '\t' | sed 's/\t$//')
MEDIAN_VALS=$(cut -f4 "$TMP_RESULTS_FILE" | tr '\n' '\t' | sed 's/\t$//')
Q3_VALS=$(cut -f5 "$TMP_RESULTS_FILE" | tr '\n' '\t' | sed 's/\t$//')

# Get the high blank proportions
HIGH_BLANK_PROPORTIONS_VALS=$(cut -f2 "$TMP_HIGH_BLANK_PROPORTIONS_FILE" | tr '\n' '\t' | sed 's/\t$//')

# Print the final, formatted TSV to standard output
echo -e "\t$HEADERS"
echo -e "Min\t$MIN_VALS"
echo -e "Q1\t$Q1_VALS"
echo -e "Median\t$MEDIAN_VALS"
echo -e "Q3\t$Q3_VALS"
echo -e "Prop_Ge_0.8860\t$HIGH_BLANK_PROPORTIONS_VALS"
