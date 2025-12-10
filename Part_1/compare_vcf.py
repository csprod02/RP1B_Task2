import subprocess
import argparse


def calc_prec_recall(job):
    # these are false negatives (unique to my mock vcf- not detected by bcftools)
    with open(f'results_{job}/isec_output/0000.vcf') as file_0:
        num_var_0 = 0
        
        for line in file_0:
            if '#' not in line:
                num_var_0 += 1
        
    # these are false positives (unique to the bcftools vc vcf- variants detected that I didn't put in there)    
    with open(f'results_{job}/isec_output/0001.vcf') as file_1:
        num_var_1 = 0
        
        for line in file_1:
            if '#' not in line:
                num_var_1 += 1
                
    # these are true positives (present in my mock vcf and also detected by the bcftools vc vcf)            
    with open(f'results_{job}/isec_output/0002.vcf') as file_2:
        num_var_2 = 0
        
        for line in file_2:
            if '#' not in line:
                num_var_2 += 1
    
    # precision = how many of the variants detected by bcftools vc are actually correct??  (i.e. tp / (tp + fp))
    precision = num_var_2 / (num_var_2 + num_var_1)
    
    # recall = how many of the 320 variants does bcftools detect??  (i.e. tp / (tp + fn) OR num_var)
    recall = num_var_2 / (num_var_2 + num_var_0)
    
    with open(f'results_{job}/{job}_prec_recall.txt', 'w') as results_file:
        results_file.write(f'Precision: {precision}; Recall: {recall}\n')
        results_file.write(f'False Negatives: {num_var_0}; False Positives: {num_var_1}; True Positives: {num_var_2}')
        
    
    return precision, recall


## Argument Parsing
def parse_args():
    parser = argparse.ArgumentParser(description='Generate a mock-mutated genome and then simulate reads.')
    
    # arguments for mutate_genome.py and simulate_reads.py
    parser.add_argument('-i', '--input', required=True, help='Path to input ref genome')
    parser.add_argument('-s', '--snps', required=False, default=300, help='Number of SNPs to make in mutated genome, default is 300')
    parser.add_argument('-id', '--indels', required=False, default=20, help='Number of indels to make in mutated genome, default is 20')
    parser.add_argument('-c', '--contig_ID', required=True, help='Contig ID (must be consistent with ref fasta seq)')
    parser.add_argument('-d', '--depth', required=False, default=30, help='Required depth, default is 30')
    parser.add_argument('-rl', '--read_length', required=False, default=100, help='Required read length, default is 100')
    parser.add_argument('-j', '--job', required=True, help='Title for job- Used in naming resulting files etc.')
        
    return parser.parse_args()


# using .run() instead of .Popen() because I want it to finish this before moving to the next part
def main():
    args = parse_args()
    
    p1 = subprocess.run(['python', 'mutate_genome.py',
                           '-i', args.input,
                           '-s', str(args.snps),
                           '-id', str(args.indels),
                           '-c', str(args.contig_ID),
                           '-j', str(args.job)], capture_output=True, text=True)
    
    genome = p1.stdout.strip()

    
    p2 = subprocess.run(['python', 'simulate_reads.py',
                         '-d', str(args.depth),
                         '-rl', str(args.read_length),
                         '-c', str(args.contig_ID),
                         '-j', str(args.job)], input=genome, text=True)
    
    
    # make master folder for everything to go into
    mkdir_master_cmd = f'mkdir results_{args.job}'
    
    p0 = subprocess.run(mkdir_master_cmd, shell=True)
    
    # run the bcftools pipeline from class
    vc_1_cmd = f'minimap2 -a -x sr {args.input} results_{args.job}/{args.contig_ID}_mutated_R1.fastq results_{args.job}/{args.contig_ID}_mutated_R2.fastq | samtools view -h -F 0x900 - | samtools sort -O bam > results_{args.job}/{args.contig_ID}_mapped_reads.bam'
    
    vc_2_cmd = f'bcftools mpileup -Ou -f {args.input} results_{args.job}/{args.contig_ID}_mapped_reads.bam | bcftools call -vc -Ov > results_{args.job}/{args.contig_ID}_variants.vcf'
    
    p = subprocess.run(vc_1_cmd, shell=True)
    pp  = subprocess.run(vc_2_cmd, shell=True)
    
    # remove other cols
    filter_cmd = f'bcftools annotate -x INFO,FORMAT -Oz -o results_{args.job}/stripped_mock_mutated_{args.contig_ID}.vcf results_{args.job}/mock_mutated_{args.contig_ID}.vcf'
    pf = subprocess.run(filter_cmd, shell=True)
    
    filter_cmd2 = f'bcftools annotate -x INFO,FORMAT -Oz -o results_{args.job}/stripped_{args.contig_ID}_variants.vcf results_{args.job}/{args.contig_ID}_variants.vcf'
    pf2 = subprocess.run(filter_cmd2, shell=True)
    

    # need to normalise and compress vcf before vcfeval
    norm_1_cmd = f'bcftools norm -f {args.input} -m -both results_{args.job}/stripped_mock_mutated_{args.contig_ID}.vcf -Oz -o results_{args.job}/mock_mutated_{args.contig_ID}_norm.vcf.gz'
    norm_2_cmd = f'bcftools norm -f {args.input} -m -both results_{args.job}/stripped_{args.contig_ID}_variants.vcf -Oz -o results_{args.job}/{args.contig_ID}_variants_norm.vcf.gz'
    
    p3_1 = subprocess.run(norm_1_cmd, shell=True)
    p3_2 = subprocess.run(norm_2_cmd, shell=True)
    
    
    # need to index them too
    index_1_cmd = f'bcftools index --tbi results_{args.job}/mock_mutated_{args.contig_ID}_norm.vcf.gz'
    index_2_cmd = f'bcftools index --tbi results_{args.job}/{args.contig_ID}_variants_norm.vcf.gz'
    
    p4_1 = subprocess.Popen(index_1_cmd, shell=True)
    p4_2 = subprocess.Popen(index_2_cmd, shell=True)
    
    p4_1.wait()
    p4_2.wait()
    
    # make a directory for results of bcftools isec to go into
    mkdir_cmd = f'mkdir results_{args.job}/isec_output'
    
    pa = subprocess.run(mkdir_cmd, shell=True)
    
    # Using bcftools isec (finds intersecting sections of vcf files)
    # Also used this forum post to help https://www.biostars.org/p/136937/
    # 0000.vcf = variants unique to my mock vcf (true vcf)
    # 0001.vcf = variants unique to the bcftools-produced vcf (false positives)
    # 0002.vcf = variants present in both vcfs
    intersect_cmd = f'bcftools isec results_{args.job}/mock_mutated_{args.contig_ID}_norm.vcf.gz results_{args.job}/{args.contig_ID}_variants_norm.vcf.gz -p results_{args.job}/isec_output'
    
    pb = subprocess.run(intersect_cmd, shell=True)
    
    precision, recall = calc_prec_recall(args.job)
    
    print(f'Precision: {precision}; Recall: {recall}')


if __name__ == '__main__':
    main()
