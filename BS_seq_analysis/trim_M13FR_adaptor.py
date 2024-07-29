import sys
import regex as re
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement(seq):
    """ Returns the reverse complement of a DNA sequence. """
    return str(Seq(seq).reverse_complement())

def trim_sequence(sequence, start_pattern, end_pattern):
    """ Trims sequence by removing everything before start_pattern and everything after end_pattern. """
    # Find the start pattern with regex for up to 1 mismatch
    start_match = re.search(f"(?e)({start_pattern}){{e<=1}}", sequence)
    if start_match:
        start_index = start_match.end()  # get the end index of the match to keep everything after the pattern
    else:
        start_index = 0  # if not found, start from beginning

    # Find the end pattern with regex for up to 1 mismatch
    end_match = re.search(f"(?e)({end_pattern}){{e<=1}}", sequence)
    if end_match:
        end_index = end_match.start()  # get the start index of the match to keep everything before the pattern
    else:
        end_index = len(sequence)  # if not found, end at the last character

    # Slice the sequence between the found indices
    return sequence[start_index:end_index]

def main(input_file):
    # Define the target sequences
    start_sequence = "CTGGAATTCGCCCTT"
    end_sequence = "AAGGGCGAATTCTGC"

    # Output filename based on the input filename
    output_file = f"trimmed_{input_file}"

    # Read the FASTA file, trim the sequences, and write to a new file
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            trimmed_seq = trim_sequence(str(record.seq), start_sequence, end_sequence)
            record.seq = Seq(trimmed_seq)
            SeqIO.write(record, output_handle, "fasta")

    print(f"Trimming complete. Output saved to '{output_file}'.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <inputfile.fasta>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    main(input_fasta)
