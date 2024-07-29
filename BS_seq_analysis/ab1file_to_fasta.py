import os
from Bio import SeqIO
import sys

def ab1_to_fasta(input_file, output_file):

    # Read the .ab1 file
    record = SeqIO.read(input_file, "abi")  # 'abi' is the format for .ab1 files
    # Write to a .fasta file
    SeqIO.write(record, output_file, "fasta")


def process_directory(directory):
    # Loop through all files in the specified directory
    for filename in os.listdir(directory):
        if filename.endswith(".ab1"):  # Check for .ab1 files
            input_path = os.path.join(directory, filename)
            output_path = os.path.join(directory, filename.replace(".ab1", ".fasta"))
            ab1_to_fasta(input_path, output_path)
            print(f"Converted {filename} to FASTA.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <path/to/ab1_files>")
        sys.exit(1)  # Exit the script if no directory path is provided

    directory_path = sys.argv[1]  # Take the directory path from the command line
    process_directory(directory_path)
