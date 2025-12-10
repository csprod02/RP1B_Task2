import random
import argparse
import sys
import subprocess

random.seed(1)


## READING SEQ FROM FASTA
def read_sequence_from_fa(filename):
    with open(filename, 'r') as file:
        # Reads the header line
        file.readline()
    
        # Creates an empty String to store the sequence from the file
        sequence = ''
    
        # For each line in the file, strip that line of any whitespace (\n, \t) and then add it to the sequence String
        for line in file:
            sequence += line.strip()

    return sequence


## DELETIONS
def section_overlapping(new_start, new_length, del_array):
    for i in del_array:
        if (new_start < i['start_index'] + i['length']) and (i['start_index'] < new_start + new_length):
            return True

    return False


def generate_deletions(genome, num_dels):
    deletions_arr = []
    
    for i in range(num_dels):
        length = random.randint(1,10)
        start_index = random.randint(1,len(genome)-1) # starts at 1 so first base is never deleted

        while section_overlapping(start_index, length, deletions_arr):
            length = random.randint(1,10)
            start_index = random.randint(1,len(genome)-1)

        deletion_info = {'type': 'deletion',
                        'start_index': start_index,
                        'length': length,
                         # vcf pos should be the base to the left of the deletion
                         'pos': start_index, # do not need to change for vcf 1-based indexing
                         # vcf format says you should include the base to the left of the deletion
                         'ref': genome[start_index-1:start_index+length],
                         'alt': genome[start_index-1]}

        deletions_arr.append(deletion_info)

    return deletions_arr


## SNPS
def check_for_deletion(base_pos, del_array):
    for i in del_array:
        if (base_pos >= i['start_index']) and (base_pos < i['start_index']+i['length']):
            return True
    return False


def generate_snps(genome, num_snps, del_array):
    n = ['A', 'T', 'G', 'C']

    snp_arr = []

    for i in range(num_snps): # how many SNPs do you want to make?
        base_index = random.randint(0, len(genome)-1) # len -1 becuase the final index is = len-1 because index starts at 0

        # keep generating a new base index until you find one that hasn't been used before (just in case- prevents the same base mutating twice)
        while any(snp['base_index'] == base_index for snp in snp_arr) or (check_for_deletion(base_index, del_array)):
            base_index = random.randint(0, len(genome)-1)

        new_base = n[random.randint(0,3)]
        # ensures that the mutation does not replace the base with the same one that was originally there
        while new_base == genome[base_index]:
            new_base = n[random.randint(0,3)]

        snp_info = {'type': 'snp',
                    'base_index': base_index,
                    'pos': base_index+1, # +1 for vcf 1-based indexing
                    'ref': genome[base_index],
                    'alt': new_base}
        
        snp_arr.append(snp_info)

    return snp_arr
    

## INSERTIONS
def generate_insertions(genome, num_insertions, del_array, snp_array):
    n = ['A', 'T', 'G', 'C']

    insertions_arr = []

    for i in range(num_insertions): # how many insertions do you want?
        insert_length = random.randint(1,10)
        
        insertion_string = ''
        
        for i in range(insert_length):
            insertion_string += n[random.randint(0,3)]
    
        insert_loc = random.randint(0,len(genome)-1)

        # need to check that the insertion is not in the middle of a deleted zone, or at the same loc as a snp
        while any(snp['base_index'] == insert_loc for snp in snp_array) or (check_for_deletion(insert_loc, del_array)):
            insert_loc = random.randint(0,len(genome)-1)
            
        insertion_info = {'type': 'insertion',
                          'length': insert_length,
                          'insertion_string': insertion_string,
                          'insert_loc': insert_loc, # left-anchor point
                          'pos': insert_loc+1, # +1 for vcf 1-based indexing
                          'ref': genome[insert_loc],
                          'alt': genome[insert_loc]+insertion_string}

        insertions_arr.append(insertion_info)

    return insertions_arr


## APPLYING MUTATIONS
def generate_mutations(original_genome, num_indels, num_snps):
    
    num_deletions = random.randint(1,num_indels) # starting at 1 ensures there is always at least..
    num_insertions = num_indels - num_deletions # .. one deletion and one insertion (never 0 of either)
    
    deletions = generate_deletions(original_genome, num_deletions)
    snps = generate_snps(original_genome, num_snps, deletions)
    insertions = generate_insertions(original_genome, num_insertions, deletions, snps)

    # need to collate all alterations and order them by index
    # then work right to left when applying them so that indexes of following alerations
    # are not affected by shifting around of bases etc.
    all_alterations = deletions + snps + insertions
    # each x.get() is nested so that if it can't find 'start_index', it will try 'base_index'
    # and then if it can't find 'base_index', it will find 'insert_loc'
    sorted_alterations = sorted(all_alterations, key= lambda x : x.get('start_index', x.get('base_index', x.get('insert_loc'))), reverse=True)

    return sorted_alterations


def apply_mutations(original_genome, sorted_alterations):
    
    mutated_genome = original_genome
    
    for i in sorted_alterations:
        if i['type'] == 'deletion':
            start = i['start_index']
            length = i['length']
            mutated_genome = mutated_genome[:start] + mutated_genome[start+length:]
            
        elif i['type'] == 'snp':
            loc = i['base_index']
            mutated_genome = mutated_genome[:loc] + i['alt'] + mutated_genome[loc+1:]

        elif i['type'] == 'insertion':
            insert_str = i['insertion_string']
            loc = i['insert_loc']
            mutated_genome = mutated_genome[:loc+1] + insert_str + mutated_genome[loc+1:]

    return mutated_genome


## WRITING TO FILES (FASTA AND VCF)
def generate_mock_fasta(mutated_genome, contig_id, directory):
    with open(f'{directory}/mutated_{contig_id}.fasta', 'w') as file:

        fasta_header = f'{contig_id}_MUTATED'

        file.write(f'{fasta_header}\n')

        # 
        for i in range(0, len(mutated_genome), 80):
            file.write(f'{mutated_genome[i:i+80]}\n')


# Want to create a mock vcf that has all of the TRUE variants (i.e. all of the mutations I have made)
# vcf files have the following columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER,and INFO
# I have tried to follow the format correctly, using . for unnkowns
def generate_mock_vcf(sorted_alterations, contig_id, genome_len, directory):
    with open(f'{directory}/mock_mutated_{contig_id}.vcf', 'w') as file:
        
        vcf_header_1 = '##fileformat=VCFv4.2\n'
        vcf_header_2 = f'##contig=<ID={contig_id},length={genome_len}>\n'
        vcf_header_3 = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'        

        file.write(vcf_header_1)
        file.write(vcf_header_2)
        file.write(vcf_header_3)

        asc_sorted_alterations = sorted(sorted_alterations, key= lambda x : x.get('pos'))
    
        for i in asc_sorted_alterations:
            pos = i['pos']
            ref = i['ref']
            alt = i['alt']
            vcf_line = f'{contig_id}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.'

            file.write(vcf_line + '\n')


## ARGUMENT PARSING
def parse_args():
    parser = argparse.ArgumentParser(description='Generate a mock-mutated genome.')
    
    parser.add_argument('-i', '--input', required=True, help='Path to input ref genome')
    parser.add_argument('-s', '--snps', required=False, default=300, help='Number of SNPs to make in mutated genome')
    parser.add_argument('-id', '--indels', required=False, default=20, help='Number of indels to make in mutated genome')
    parser.add_argument('-c', '--contig_ID', required=True, help='Contig ID (must be consistent with ref fasta seq)')
    parser.add_argument('-j', '--job', required=False, default='Task2_Results', help='Title for job- Used in naming resulting files etc.')
        
    return parser.parse_args()


## MAIN
def main():
    args = parse_args()  # Parse command-line arguments
    
    folder = f'results_{args.job}'
        
    mkdir_master_cmd = f'mkdir {folder}'
    p0 = subprocess.run(mkdir_master_cmd, shell=True)
    
    original_genome = read_sequence_from_fa(args.input)
    
    mutations_list = generate_mutations(original_genome, int(args.indels), int(args.snps))
    mutated_genome = apply_mutations(original_genome, mutations_list)

    generate_mock_fasta(mutated_genome, args.contig_ID, folder)
    generate_mock_vcf(mutations_list, args.contig_ID, len(mutated_genome), folder)
    
    return mutated_genome


    
if __name__ == '__main__':
    mutated_genome = main()
    
    # only prints if stdout is not the terminal (i.e. will print if being piped)
    if not sys.stdout.isatty():
        print(mutated_genome)
