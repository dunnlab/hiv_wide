#!/usr/bin/env python
#
# written with gemini and modified to get formatting code to work
# analyze_hiv_data.py
#
# Description:
# This script reads a tab-delimited text file containing HIV sequence metadata
# (like the one from the LANL HIV database). It calculates the proportion
# of sequences derived from proviral DNA versus viral RNA and identifies the
# publication ID associated with the most sequences.
# This version uses the pandas library for robust data parsing.
#
# Usage:
# python analyze_hiv_data.py results.txt
#
# Prerequisite:
# You must have pandas installed. If not, run:
# pip install pandas
#

import argparse
import sys
import pandas as pd

def analyze_sequence_metadata(filepath):
    """
    Analyzes a tab-delimited file using pandas to determine molecule type
    proportions and the most frequent publication ID.

    Args:
        filepath (str): The path to the input data file.

    Returns:
        A dictionary containing the analysis results, or None on error.
    """
    try:
        # Use pandas to read the tab-separated file.
        # skiprows=[0] is used to skip the first line (e.g., "Number of records retrieved: ...")
        # The next line is correctly inferred as the header.
        df = pd.read_csv(filepath, sep='\t', skiprows=[0,1])

    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An error occurred while reading the file with pandas: {e}", file=sys.stderr)
        return None

    if df.empty:
        print("Warning: No data records were found in the file.", file=sys.stderr)
        return None

    # --- Analysis using pandas ---

    # Get total record count from the DataFrame shape
    total_records = len(df)

    first_row = df.head(1)

    # Use .value_counts() to get the counts of each 'Molecule type'
    molecule_counts = df['Molecule type'].value_counts()
    dna_count = molecule_counts.get('DNA', 0)
    rna_count = molecule_counts.get('RNA', 0)

    # Use .value_counts() on the publication ID column.
    # .dropna() removes any missing values before counting.
    pub_id_counts = df['PUB id(SPL)'].dropna().value_counts()

    # Find the most common publication ID and its count
    if not pub_id_counts.empty:
        most_common_pub_id = pub_id_counts.index[0]
        most_common_pub_count = pub_id_counts.iloc[0]
    else:
        most_common_pub_id = "N/A"
        most_common_pub_count = 0

    # Prepare results
    results = {
        "total_records": total_records,
        "dna_count": int(dna_count),
        "rna_count": int(rna_count),
        "dna_proportion": (dna_count / total_records) * 100 if total_records > 0 else 0,
        "rna_proportion": (rna_count / total_records) * 100 if total_records > 0 else 0,
        "most_common_pub_id": most_common_pub_id,
        "most_common_pub_count": int(most_common_pub_count)
    }
    
    return results

def main():
    """
    Main function to parse arguments and print the analysis.
    """
    parser = argparse.ArgumentParser(
        description="Analyzes HIV sequence metadata to find molecule type proportions and top publication."
    )
    parser.add_argument(
        "filepath",
        help="Path to the tab-delimited results file (e.g., results.txt)."
    )
    args = parser.parse_args()

    analysis = analyze_sequence_metadata(args.filepath)

    if analysis:
        print("--- HIV Sequence Data Analysis ---")
        print(f"\nTotal sequences analyzed: {analysis['total_records']}")
        
        print("\n[Molecule Type Breakdown]")
        print(f"Proviral DNA sequences: {analysis['dna_count']} ({analysis['dna_proportion']:.2f}%)")
        print(f"Viral RNA sequences:    {analysis['rna_count']} ({analysis['rna_proportion']:.2f}%)")
        
        print("\n[Publication Analysis]")
        print(f"Publication ID with the most sequences: '{analysis['most_common_pub_id']}'")
        print(f"Number of sequences from this publication: {analysis['most_common_pub_count']}")
        print("\n----------------------------------")


if __name__ == "__main__":
    main()
