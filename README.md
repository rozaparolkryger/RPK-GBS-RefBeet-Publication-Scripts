# RPK-GBS-RefBeet-Publication-Scripts

Summary:

Following scripts have been created within GBS pipeline for genome assembly improvement.


Create_ABH.pl
Script for generation of SNP table, contiaining SNP information for all individuals with available VCF files. 
This table genarate basis file was basis file used for anchoring approach.

INPUT: txt files converted from VCF files (for each individuum)

OUTPUT: ABH_table.txt 
EXEC:


ABH_cover_Pos.pl

Script for assignment of genotypes, coverage verification and noise reduction.
This script needs pileup file containing a coverage information for each SNP position in the reference sequence as well as the primarily constructed ABH table. 

Reference allele has been assign to 'A', non-reference allele to ‘B’ and both alleles to ‘H’. The positions with no SNP evidence were labeled as missing allele. 

To avoid noise while calling the genotype along chromosomes of each individual, the created ABH table was scanned for the expected accuracy. The SNP positions with over 33% missing information or those showing vast deviation from expected segregation pattern (25%A:50%H:25%B), like for example position containing >42% of alternative allele were searched using bash commands. All those matched SNP positions were removed from the resulting ABH table file. 


INPUT: pileup file and ABH_table.txt
OUTPUT:
EXEC:

evaluate_ABH.py

The ABH table was evaluated using a sliding window approach to call genotypes along the scaffolds aligned to chromosomes for each individual. Genotypes were called based on the allele ratio in one window of certain length and shifted by an overlap. A homozygous genotype for the reference was identified, if the proportion was greater than 90%, a homozygous genotype for alternative allele (K1P2) was called, when the proportion of SNPs was over 80%. If the proportion of missing information was over 60%, no genotype was called and in all other cases a heterozygous genotype was automatically assigned. 

INPUT: 
OUTPUT:
EXEC:

pattern_search.py

This script has simultaneously compared all pattern from query list with all patterns from core pattern list. The anchoring algorithm were used to search for the similarities of the pattern due to different thresholds. The results were divided in the three matching groups according to the similarity threshold. Perfect matches presented 100 % similarity and were split to unique match (query sequence hits only one core pattern) or multi-match (query sequence hits more than one core pattern). Optimal matches were considered with over 95 % identity of the pattern and also divided in unique and multi-matches. The last group of suboptimal matches were achieved by similarity of > 80 % and < 95 % and was useful only to define anchoring position to at least linkage group for query patterns.


INPUT: 
OUTPUT:
EXEC:

