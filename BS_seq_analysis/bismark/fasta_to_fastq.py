import os
from Bio import SeqIO
import sys


def fasta_to_fastq(fasta_file, fastq_file):
    with open(fasta_file, 'r') as fasta, open(fastq_file, 'w') as fastq:
        for line in fasta:
            if line.startswith('>'):
                # Write the sequence identifier line for FASTQ
                fastq.write('@' + line[1:])
            else:
                sequence = line.strip()
                # Write the sequence line for FASTQ
                fastq.write(sequence + '\n')
                # Write the plus separator line for FASTQ
                fastq.write('+\n')
                # Write the quality score line for FASTQ, with 'I' for a dummy quality score of 40
                quality_score = 'I' * len(sequence)
                fastq.write(quality_score + '\n')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output.fastq")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_fastq = sys.argv[2]
    fasta_to_fastq(input_fasta, output_fastq)