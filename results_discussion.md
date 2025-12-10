# Discussion of Results
## E-coli
### Part 1:
***CLI Command:*** python compare_vcf.py -i EcoliK12-MG1655.fasta -c NC_000913.3 -j Ecoli  
**Precision:** 1.0
**Recall:** 0.978  
The recall of bcftools alone was ~0.98 and there was a total of 7 false negatives (mutations missed by the variant caller).  When looking at the 0000.vcf file for these results, it appears that all of these missed variants were SNPs.  When using samtoolsl tview to view these locations, it seems that these missing SNPs appear to be only present in some of the reads, but not all mapped to that location.  Considering these are simulated reads, this doesn’t appear to make much sense, however I assume that due to repetitive sections within the genome, some reads could have matches in multiple places with minimal mismatch and so these variants have maybe been missed as a result.  

### Part 2:
#### Simulated Data
***CLI Command:*** python combine_vcf.py -i EcoliK12-MG1655.fasta -r1 NC_000913.3_mutated_R1.fastq -r2 NC_000913.3_mutated_R2.fastq -c NC_000913.3 -j Ecoli_sim  
**Combined Precision:** 1.0
**Combined Recall:** 0.978  
***Variants unique to bcftools:*** 4  
***Variants unique to snippy:*** 0  
***Variants found by both:*** 309  
The precision and recall of the combined pipeline were equal to bcftools alone, with a precision of 1.0 and a recall of ~0.978.  A total of 313 out of 320 variants were detected across both snippy and bcftools combined.  This was the same as bcftools alone.  
When looking at the results in the Ecoli_sim_vc_compare.txt file, you can see that bcftools actually detected 4 variants that snippy missed.  These results also show that snippy failed to detect any additional variants that bcftools missed.  
#### Real Data
***CLI Command:*** python combine_vcf.py -i EcoliK12-MG1655.fasta -r1 SRR25083113_1.fastq.gz -r2 SRR25083113_2.fastq.gz -c NC_000913.3 -j Ecoli_real  
The Ecoli_real_combined_prec_recall.txt can be ignored as it is comparing/combining with the simulated data as part of the pipeline.  The precision and recall here shows as 0.0 for both because there isn’t any overlap with the truth vcf file because the reads are from a different real dataset.  
However, when looking at the Ecoli_real_vc_compare.txt file you can see that both bcftools and snippy detected 47,354 variants.  There were an additional 22,609 variants detected by bcftoools alone, and a further 8,557 variants detected by snippy alone.  This would bring the total number of variants to 78,520.  

## Plasmodium
### Part 1:
***CLI Command:*** python compare_vcf.py -i NC_037282.1.fasta -c NC_037282.1 -j plasmodium  
**Precision:** 1.0
**Recall:** 1.0  
When running my pipeline on the plasmodium genome, both the precision and recall of the bcftools variant caller were 1.0.  This means that the variant caller detected all of the simulated variants, and that there were no additional variants detected (i.e. no false positives).  Of all the variants detected, all were true positives.  

### Part 2:
***CLI Command:*** python combine_vcf.py -i NC_037282.1.fasta -r1 NC_037282.1_mutated_R1.fastq -r2 NC_037282.1_mutated_R2.fastq -c NC_037282.1 -j plasmodium  
**Combined Precision:** 0.997
**Combined Recall:** 1.0  
***Variants unique to bcftools:*** 6  
***Variants unique to snippy:*** 1  
***Variants found by both:*** 314  
When comparing the combined vcf with the truth vcf, the precision of the ‘combined variant caller’ (bcftools and snippy) was actually slightly lower than bcftools alone.  
This combined method found a total of 321 variants and given that there was a total of 320 mutations applied to the plasmodium genome (300 SNPs and 20 indels), it can be assumed that at least one of these detected variants must be a false positive.  
Considering that the bcftools only pipeline had both a precision and recall of 1.0 each, we can assume that the introduction of snippy could be resulting in this error.  There were 6 variants detected by bcftools which were not found by snippy, and we know that these must be genuine mistakes because they were correctly detected previously using the same bcftools pipeline during Part 1.  Because the same random seed was used each time either pipeline was run, we can eliminate the possibility that these differences could be resulting from differences in mutations applied to the reference genome when simulating reads.  There was one false positive reported in the plasmodium_combined_prec_recall.txt file and it makes sense that this would be the additional variant detected by snippy.  

## How much do I trust each method?
When looking at the plasmodium data, bcftools performed very well, finding all variants with no false positives, or negatives.  The combined pipeline had a marginally lower precision due to an extra variant detected by snippy which was a probable false positive.  
When looking at the simulated E-coli data, bcftools performed well, and the combined pipeline did not appear to add much benefit since precision and recall were the same as bcftools alone.  There were no additional variants detected when compared to bcftools alone.  
Therefore, when looking at the real E-coli data, I would be hesitant to rely on the combined pipeline due to the reduction in precision seen in the plasmodium data, and the absence of any additional benefit seen with the simulated E-coli data.  However, this is only a small sample of results.  To get a better understanding of the precision and recall, I could run these pipelines again on more sets of mutated genomes and simulated reads and take an average across these repeats.  This would give me more confidence in using the combined pipeline with snippy.  
The real reads also have many more variants detected (~78k vs ~320 in simulated dataset) by the combined pipeline.  Bcftools detected a total of 69,963 variants, whilst snippy detected a total of 55,911.  This is an additional 14,052 variants detected by bcftools.  There were 47,354 variants detected by both variant callers.  Testing the pipelines on more simulated data but with a significantly higher number of variants could also be interesting to see how this could potentially affect precision and recall scores.
