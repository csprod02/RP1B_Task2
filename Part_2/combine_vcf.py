import subprocess
import argparse


def calc_prec_recall(job):
    # these are false negatives (unique to my mock vcf- not detected by bcftools)
    with open(f'results_{job}/isec_output_sim/0000.vcf') as file_0:
        num_var_0 = 0
        
        for line in file_0:
            if '#' not in line:
                num_var_0 += 1
        
    # these are false positives (unique to the bcftools vc vcf- variants detected that I didn't put in there)    
    with open(f'results_{job}/isec_output_sim/0001.vcf') as file_1:
        num_var_1 = 0
        
        for line in file_1:
            if '#' not in line:
                num_var_1 += 1
                
    # these are true positives (present in my mock vcf and also detected by the bcftools vc vcf)            
    with open(f'results_{job}/isec_output_sim/0002.vcf') as file_2:
        num_var_2 = 0
        
        for line in file_2:
            if '#' not in line:
                num_var_2 += 1
    
    # precision = how many of the variants detected by bcftools vc are actually correct??  (i.e. tp / (tp + fp))
    precision = num_var_2 / (num_var_2 + num_var_1)
    
    # recall = how many of the 320 variants does bcftools detect??  (i.e. tp / (tp + fn) OR num_var)
    recall = num_var_2 / (num_var_2 + num_var_0)
    
    print(f'Precision: {precision}; Recall: {recall}')
    
    with open(f'results_{job}/{job}_combined_prec_recall.txt', 'w') as results_file:
        results_file.write(f'Precision: {precision}; Recall: {recall}\n')
        results_file.write(f'False Negatives: {num_var_0}; False Positives: {num_var_1}; True Positives: {num_var_2}')
        
    
    return precision, recall

def calc_num_vars(job):
    # these are detected only by bcftools
    with open(f'results_{job}/isec_output_vc/0000.vcf') as file_0:
        num_var_0 = 0
        
        for line in file_0:
            if '#' not in line:
                num_var_0 += 1
          
    # these are detected only by snippy    
    with open(f'results_{job}/isec_output_vc/0001.vcf') as file_1:
        num_var_1 = 0
        
        for line in file_1:
            if '#' not in line:
                num_var_1 += 1
    
    # these are detected by both
    with open(f'results_{job}/isec_output_vc/0002.vcf') as file_2:
        num_var_2 = 0
        
        for line in file_2:
            if '#' not in line:
                num_var_2 += 1
                
    print(f'bcftools: {num_var_0}; snippy: {num_var_1}; both: {num_var_2}')
    
    with open(f'results_{job}/{job}_vc_compare.txt', 'w') as results_file:
        results_file.write(f'variants unique to bcftools: {num_var_0}\nvariants unique to snippy: {num_var_1}\nvariants found by both: {num_var_2}')
                
    return
                


def parse_args():
    parser = argparse.ArgumentParser(description='Run two variant callers (bcftools and snippy) and combine the resulting vcf files.')
    
    # arguments for mutate_genome.py and simulate_reads.py
    parser.add_argument('-i', '--input', required=True, help='Path to input ref genome')
    parser.add_argument('-r1', '--r1_fastq', required=True, help='Path for illumina reads fastq files R1')
    parser.add_argument('-r2', '--r2_fastq', required=True, help='Path for illumina reads fastq files R1')
    parser.add_argument('-c', '--contig_ID', required=True, help='Contig ID (must be consistent with ref fasta seq)')
    parser.add_argument('-j', '--job', required=True, help='Title for job- Used in naming resulting files etc.')

    return parser.parse_args()
    

def main():
    args = parse_args()
    
    # make master folder for everything to go into
    mkdir_master_cmd = f'mkdir results_{args.job}'
    
    p0_1 = subprocess.run(mkdir_master_cmd, shell=True)
    
    # make a folder for the bcftools results to go into (makes it clearer which vcf has come from which variant caller)
    mkdir_cmd = f'mkdir results_{args.job}/bcftools_vcf_results'
    
    p0_2 = subprocess.run(mkdir_cmd, shell=True) # .run() because I want this to be finished before I run anything else
    
    # variant caller 1 (bcftools)
    vc_1_cmd = f'minimap2 -a -x sr {args.input} {args.r1_fastq} {args.r2_fastq} | samtools view -h -F 0x900 - | samtools sort -O bam > results_{args.job}/bcftools_vcf_results/mapped_reads.bam'
    vc_2_cmd = f'bcftools mpileup -Ou -f {args.input} results_{args.job}/bcftools_vcf_results/mapped_reads.bam | bcftools call -vc -Ov > results_{args.job}/bcftools_vcf_results/variants.vcf'
    
    p1_1 = subprocess.Popen(vc_1_cmd, shell=True)
    
    p1_1.wait() # wait for p1_1 to complete before starting p1_2
    
    p1_2 = subprocess.Popen(vc_2_cmd, shell=True)
    
    # variant caller 2 (snippy)
    # https://github.com/tseemann/snippy
    # added to my conda environment using the instructions on GitHub
    # adapted the quick start cli script (assumed it was okay to remove the part about cpu usage) and added --cleanup because it produces a ton of files otherwise!!
    vc2_cmd = f'snippy --cleanup --outdir results_{args.job}/snippy_vcf_results --ref {args.input} --R1 {args.r1_fastq} --R2 {args.r2_fastq}'
    
    p2 = subprocess.Popen(vc2_cmd, shell=True)
    
    # waiting for the variant callers to both finish running
    p1_2.wait()
    p2.wait()
    
    ## Combining vcf files (bcftools merge)
    
    # need to normalise, compress, and index the vcf files first
    # https://samtools.github.io/bcftools/bcftools.html#norm
    # -f is the ref fasta path
    # -m -both means that both snps and indels will be split (-) into multiple records (makes it easier for combining vcfs?)
    # -Oz means that the output type will be compressed vcf
    # -o is the output file name
    norm_1_cmd = f'bcftools norm -f {args.input} -m -both results_{args.job}/bcftools_vcf_results/variants.vcf -Oz -o results_{args.job}/bcftools_vcf_results/bcftools_norm.vcf.gz'
    norm_2_cmd = f'bcftools norm -f {args.input} -m -both results_{args.job}/snippy_vcf_results/snps.vcf -Oz -o results_{args.job}/snippy_vcf_results/snippy_norm.vcf.gz'
    
    p3_1 = subprocess.Popen(norm_1_cmd, shell=True)
    p3_2 = subprocess.Popen(norm_2_cmd, shell=True)
    
    # waiting for these to finish before passing to the next stage
    p3_1.wait()
    p3_2.wait()
    
    # index both .vcf.gz files
    # https://samtools.github.io/bcftools/bcftools.html#index
    index_1_cmd = f'bcftools index results_{args.job}/bcftools_vcf_results/bcftools_norm.vcf.gz'
    index_2_cmd = f'bcftools index results_{args.job}/snippy_vcf_results/snippy_norm.vcf.gz'
    
    p4_1 = subprocess.Popen(index_1_cmd, shell=True)
    p4_2 = subprocess.Popen(index_2_cmd, shell=True)
    
    # wait for both files to be indexed before passing to bcftools merge
    p4_1.wait()
    p4_2.wait()
    
    # then use bcftools merge
    # https://samtools.github.io/bcftools/bcftools.html#merge
    # using merge instead of concat because these are from the same region of the genome (concat would just append one vcf at the end of the other?)
    merge_cmd = f'bcftools merge results_{args.job}/bcftools_vcf_results/bcftools_norm.vcf.gz results_{args.job}/snippy_vcf_results/snippy_norm.vcf.gz > results_{args.job}/combined.{args.job}.vcf'
    
    p5 = subprocess.run(merge_cmd, shell=True)
    
    # remove other cols (my mock vcf only has CHROM, POS, ID, REF, ALT, FILTER, INFO)
    # I am only interested in the POS, REF, and ALT columns so the extra info is not needed and it makes it easier to combine since my mock vcf is missing lots of info/cols
    filter_cmd = f'bcftools annotate -x INFO,FORMAT -Oz -o results_{args.job}/stripped_mock_mutated_{args.contig_ID}.vcf mock_mutated_{args.contig_ID}.vcf'
    pf = subprocess.run(filter_cmd, shell=True)
    
    filter_cmd2 = f'bcftools annotate -x INFO,FORMAT -Oz -o results_{args.job}/stripped_combined.{args.job}.vcf results_{args.job}/combined.{args.job}.vcf'
    pf2 = subprocess.run(filter_cmd2, shell=True)
    

    # need to normalise and compress vcf before vcfeval
    norm_1_cmd = f'bcftools norm -f {args.input} -m -both results_{args.job}/stripped_mock_mutated_{args.contig_ID}.vcf -Oz -o results_{args.job}/mock_mutated_{args.contig_ID}_norm.vcf.gz'
    norm_2_cmd = f'bcftools norm -f {args.input} -m -both results_{args.job}/stripped_combined.{args.job}.vcf -Oz -o results_{args.job}/stripped_combined.{args.job}_norm.vcf.gz'
    
    p3_1 = subprocess.run(norm_1_cmd, shell=True)
    p3_2 = subprocess.run(norm_2_cmd, shell=True)
    
    
    # need to index them too
    index_1_cmd = f'bcftools index --tbi results_{args.job}/mock_mutated_{args.contig_ID}_norm.vcf.gz'
    index_2_cmd = f'bcftools index --tbi results_{args.job}/stripped_combined.{args.job}_norm.vcf.gz'
    
    p4_1 = subprocess.Popen(index_1_cmd, shell=True)
    p4_2 = subprocess.Popen(index_2_cmd, shell=True)
    
    p4_1.wait()
    p4_2.wait()
    
    
    # Using bcftools isec (finds intersecting sections of vcf files)
    # Also used this forum post to help https://www.biostars.org/p/136937/
    # 0000.vcf = variants unique to my mock vcf (true vcf)
    # 0001.vcf = variants unique to the bcftools-produced vcf (false positives)
    # 0002.vcf = variants present in both vcfs
    intersect_cmd = f'bcftools isec results_{args.job}/mock_mutated_{args.contig_ID}_norm.vcf.gz results_{args.job}/stripped_combined.{args.job}_norm.vcf.gz -p results_{args.job}/isec_output_sim'
    
    pb = subprocess.run(intersect_cmd, shell=True)
    
    calc_prec_recall(args.job)
    
    # Now find how many variants were found by bcftools alone vs snippy alone, vs found by both
    # Can use bcftools isec for this too
    # in this case...
    # 0000.vcf is everything unique to the bcftools vcf
    # 0001.vcf is everything unique to the snippy vcf
    # 0002.vcf is variants present in both bcftools vcf and snippy vcf
    intersect_cmd_2 = f'bcftools isec results_{args.job}/bcftools_vcf_results/bcftools_norm.vcf.gz results_{args.job}/snippy_vcf_results/snippy_norm.vcf.gz -p results_{args.job}/isec_output_vc'
    
    pc = subprocess.run(intersect_cmd_2, shell=True)
    
    calc_num_vars(args.job)
    



if __name__ == '__main__' :
    main()
    
 