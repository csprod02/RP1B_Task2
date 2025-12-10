import random
import argparse
import sys
import subprocess

random.seed(1)

## READING SEQ FROM FASTA
def read_sequence_from_fa(filename):
    with open(filename, 'r') as file:
        # Reads the header line (we don't need this info)
        file.readline()
    
        # Creates an empty String to store the sequence from the file
        sequence = ''
    
        # For each line in the file, strip that line of any whitespace (\n, \t) and then add it to the sequence String
        for line in file:
            sequence += line.strip()

    return sequence


## SIMULATE PAIRED READS
def reverse_complement(dna_seq):
    # Flips the entire String i.e. hello -> olleh
    reverse_seq = dna_seq[::-1]

    # Creates an empty string to hold the reverse complement sequence/String
    reverse_complement_seq = ''

    # Creates a dictionary containing the complementary base
    complement = {'A': 'T',
                 'T': 'A',
                 'G': 'C',
                 'C': 'G'}

    # For each base in the reversed sequence, add its complement to the reverse_complement String
    for i in reverse_seq:
        reverse_complement_seq += complement[i]

    return reverse_complement_seq


# Illumina reads are paired
# fragments are taken and then reads are taken from both ends (one in 5' to 3' and the other in the 3' to 5' direction)
# The fragments must be longer than 2 x read_length
# in this example, reads are 100bp, so fragments can be 300bp
def simulate_paired_reads(genome, depth, read_length):
    # E.g. 1000bp seq and read len is 100, 10 reads will cover, so 300 will cover 30x depth
    num_reads = (len(genome)/read_length)*depth
    num_pairs = int(num_reads/2)

    # [[r1, r2], [r1, r2], [r1,r2], etc.]
    paired_reads_arr = []

    for i in range(num_pairs):
        start = random.randint(0,len(genome)-1-read_length)
        end = start + (2*read_length)+100
        frag = genome[start:end]

        read_1 = frag[0:read_length]
        read_2 = reverse_complement(frag[-read_length:])
        
        read_1_info = {'header_line': f'@read_{i}/1',
                        'sequence_line': read_1,
                        'separator_line': '+',
                        'quality_line': read_length*'I'} # I indicates perfect score?
    
        read_2_info = {'header_line': f'@read_{i}/2',
                        'sequence_line': read_2,
                        'separator_line': '+',
                        'quality_line': read_length*'I'}

        paired_reads_arr.append([read_1_info, read_2_info])

    return paired_reads_arr


## Writing to Files
def generate_mock_fastq(paired_reads_array, contig_id, directory):
    with open(f'{directory}/{contig_id}_mutated_R1.fastq', 'w') as r1_file, open(f'{directory}/{contig_id}_mutated_R2.fastq', 'w') as r2_file:

        for read_1, read_2 in paired_reads_array:
            r1_header = read_1['header_line']
            r1_seq = read_1['sequence_line']
            r1_sep = read_1['separator_line']
            r1_qual = read_1['quality_line']

            r1_fastq_entry = f'{r1_header}\n{r1_seq}\n{r1_sep}\n{r1_qual}\n'

            r2_header = read_2['header_line']
            r2_seq = read_2['sequence_line']
            r2_sep = read_2['separator_line']
            r2_qual = read_2['quality_line']

            r2_fastq_entry = f'{r2_header}\n{r2_seq}\n{r2_sep}\n{r2_qual}\n'

            r1_file.write(r1_fastq_entry)
            r2_file.write(r2_fastq_entry)


## Argument Parsing
def parse_args():
    parser = argparse.ArgumentParser(description='Simulate illumina reads (paired, short reads)')
    
    parser.add_argument('-i', '--input', required=False, help='Path to input genome to generate reads')
    parser.add_argument('-d', '--depth', required=False, default=30, help='Required depth, default is 30')
    parser.add_argument('-rl', '--read_length', required=False, default=100, help='Required read length, default is 100')
    parser.add_argument('-c', '--contig_ID', required=True, help='Contig ID (must be consistent with ref fasta seq)')
    parser.add_argument('-j', '--job', required=False, help='Title for job- Used in naming resulting files etc.')
    
    return parser.parse_args()


## MAIN
def main():
    args = parse_args()
    
    folder = f'results_{args.job}'
        
    mkdir_master_cmd = f'mkdir {folder}'
    p0 = subprocess.run(mkdir_master_cmd, shell=True)

    
    # if input is not text from the terminal (i.e. file path), use stdin
    if not sys.stdin.isatty():
        genome = sys.stdin.read()
    
    else: # If not using stdin, use args.input
        genome = read_sequence_from_fa(args.input)

        
    paired_reads = simulate_paired_reads(genome, int(args.depth), int(args.read_length))

    generate_mock_fastq(paired_reads, args.contig_ID, folder)
    
    
    
if __name__ == '__main__':
    main()
