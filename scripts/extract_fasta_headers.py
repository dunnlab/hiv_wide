#!/usr/bin/env python
#
# fasta_header_extractor.py
#
# Description:
# This script reads a FASTA formatted file and extracts all header lines.
# It then cleans these headers by removing the leading '>' character and any
# surrounding whitespace, printing just the sequence ID to the console.
#
# Usage:
# python fasta_header_extractor.py your_file.fasta > sequence_ids.txt
#
# The output will be a list of sequence IDs, one per line, which you can
# redirect to a text file as shown in the usage example.
#

import argparse
import sys

def extract_fasta_headers(fasta_filepath):
    """
    Reads a FASTA file and yields clean sequence IDs.

    This function opens a FASTA file, iterates through it line by line,
    and yields each header it finds after cleaning it. Using a generator
    is memory-efficient, as it doesn't load the whole file into memory.

    Args:
        fasta_filepath (str): The full path to the input FASTA file.

    Yields:
        str: A cleaned sequence ID for each header line found.
    """
    try:
        with open(fasta_filepath, 'r') as fasta_file:
            for line in fasta_file:
                # Header lines in FASTA format start with a '>'
                if line.startswith('>'):
                    # 1. Remove the leading '>' character using lstrip.
                    # 2. Remove any leading/trailing whitespace (like newlines) using strip.
                    clean_header = line.lstrip('>').strip()
                    yield clean_header
    except FileNotFoundError:
        # Provide a user-friendly error message if the file doesn't exist
        print(f"Error: The file '{fasta_filepath}' was not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        # Handle other potential file reading errors
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """
    Main function to parse command-line arguments and run the extraction.
    """
    # Set up the argument parser to handle command-line inputs
    parser = argparse.ArgumentParser(
        description="Extracts and cleans all sequence IDs from a FASTA file."
    )
    parser.add_argument(
        "fasta_file",
        help="Path to the input FASTA file."
    )
    args = parser.parse_args()

    # Get the generator and loop through it, printing each ID
    header_generator = extract_fasta_headers(args.fasta_file)
    for header in header_generator:
        print(header)

if __name__ == "__main__":
    main()
