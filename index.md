# prediction target genes of miRNA using Lasso

## Data preparation
### Data resources

|Resource | URL|
| ------------- | ------------- |
|Annotated genes | https://www.gencodegenes.org/human/release_32lift37.html|
|RNA-seq by tissues | https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1733/samples/?s_page=1&s_pagesize=500|
|FANTOM by tissues | https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz|

### Converting Data Format

| Convert | Command |
| ------------- | ------------- |
| fastq → sam (hi-sat) | hisat2 -p 4 -x genome_tran -1 {fastq1} – 2 {fastq2} -S {sam} |
| sam → bam (samtools) | samtools sorted -@ 8 -o {bam} {sam} |
| bam → gtf (stringtie)	| stringtie -p 4 -G genes.gtf -o {gtf} -i {bam} |

### File entries by Tissues
| resource	| URL |
| ------------- | ------------- |
| FANTOM	| https://drive.google.com/open?id=1y3w-rOvgbxtlCKWss1wAzcHD884Na__tijNWL7Yjb2I |
| RNA-seq	| https://drive.google.com/open?id=1lKZN5hn2e6Zq3eKzn5r3H5NqdEfjgUmy |

### 22 common tissues between FANTOM and RNA-seq
The common tissues are listed in the below table. 
| RNA-seq |	FANTOM |
| ------------- | ------------- |
| appendix | appendix |
| urinarybladder | bladder |
| brain	| brain|
| colon	| colon|
| esophagus	| esophagus|
| fat	| adipose|
| gallbladder	| gall_bladder|
| heart	| heart|
| kidney | kidney|
| liver	| liver|
| lung	| lung|
| lymphnode	| lymph_node|
| ovary	| ovary|
| pancreas	| pancreas|
| placenta	| placenta|
| prostate	| prostate|
| salivarygland	| salivary_gland|
| smallintestine	| small_intestine|
| spleen	| spleen|
| testis	| testis|
| thyroid	| thyroid|
| endometrium	| uterus|

### Pre-processing
#### Averaged data by tissue
One tissue can have multiple replicates. Therefore, RNA-seq and FANTOM were averaged on same tissues. The average were calculated by the below equation. For this, we  only used chr1, chr2, …chr 22, chrX and chrY. The average is calculated by each chromosome and strand.

#### Expression data by 22 tissues
RNA-seq data already contains FPKM while FANTOM does not have expression data. To get the expression data from FANTOM, we calculated expression level of (+/-)100 bp from annotated gene TSS. The expression level is summed by the score sum of cage tags within the region. When the summation, every tags, which are overlapped to the region;  [TSS-100, TSS+100] should be used.


#### High consistent genes
RNA-seq and FANTOM data have expression level by a gene on the 22 tissues. Each gene has two 22 length vectors (RNA-seq and FANTOM), which an element show expression level on a tissue. By the two vectors, the correlation coefficient is calculated to observe how consistent between RNA-seq and FANTOM on a gene. Each gene has the coefficient and high coefficient shows two data source are consistent on the 22 tissues. For the correlation, spearman method was used. We round the coefficient at the decimal point with two digit. To extract only high correlated one, we set threshold as 0.75. Finally, 4,781 transcripts are received as high consistent ones. the transcripts are grouped by gene name and if there are multiple transcripts, we picked one with the maximum score. After this, we got 3,519 genes. The list of 3,519 transcripts is in the below link. only double type is used for every calculation. Transcript region are only used that gene and transcript type are protein coding gene from genecode transcripts.  
https://drive.google.com/open?id=1FAYRAa746bWeeN6G3tsKlTJ6CDXRNkp3


#### High consistent miRNA
For the consistent miRNAs, we compared eleven papers and retrieved the TSS locations by the miRNA name. The consistent miRNA TSSs are predicted by all papers and the TSS regions are close then we chose the meddle point as the TSS. The number of miRNA TSSs is 330. The name and TSS of 330 miRNAs are in the below link.  
https://drive.google.com/open?id=1qZHBqubcYeJfVk7uvE8lFyzVXZeSknmIOf-qfS5zBmM


### Comparison
#### CAGE raw data vs. CAGE tag processed data
#### Before  vs. After Filter

##### There are three filters

maximum value of FANTOM or RNA-seq on 22 tissues is zero
median of FANTOM or RNA-seq on 22 tissues is zero
maximum value is greater than sum of others

##### before filter

The total number of transcripts for 51,686 out of 83,866.

83,866 is the number of transcripts from GENCODE. The reason why only 51,686 transcripts exist if RNA-seq data is missed in specific tissue, it is skipped. So I only took when every 22 tissue data available in RNA-seq.

##### after filter

The total number of transcripts for 31,344 out of 83,866. the procedure is same as before filter but the three filters.


#### RNA-seq data comparison
First, processed cage data are same as Amlan’s and mine.
Therefore, we checked RNA-seq data.
I aligned sequence using same reference gene and compared but it is still different.
So, we think only reason is stringtie version difference.
I updated stringtie as what Amlan use and sent the updated gtf file to Amlan.

## Methodology
### FANTOM cell line data
we received cell line specific data. The below table shows the url of the data source and the table including file id, cell line and so on.

| label | URL|
| ------------- | ------------- |
| FANTOM url | http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.tissue.hCAGE/|
| cell line data list | https://drive.google.com/open?id=1frNwIbWzFwWdvKfThP8DT1mdbPfGRERy|

### Expression data by 240 cell lines
sum of cell line data have multiple replications. Thus, we averaged if replications exist as we did before. From the data preparation step, we got high consistent transcripts on gene and miRNA. We extracted the expression level from 240 cell lines. The score sum of the region where (+/-)100bp from annotated gene TSS is the expression level. We calculated this levels about the high consistent genes and miRNAs.

### Lasso linear regression
To predict target genes on a miRNA, we used lasso regression. the consistent genes and miRNAs have 240 vector, which are expression level on each cell line. Assume the number of genes and miRNA are m and n respectively. gene matrix **X**: (m x 240) and miRNA matrix **Y**: (n x 240). Lasso calculated the relationship between **X** and **Y**. The coefficient matrix by the Lasso result is **B** (n x m); **Y = BX**. Before, processing Lasso, **X** and **Y** were centered by subtracting mean value by each row. 

From the coefficient matrix, each row shows coefficient of one miRNA and multiple genes. Therefore, the non-zero coefficient values are assumed target genes, which are related to the corresponding miRNA.

Figure 1 shows the input and output data of Lasso regression and target genes.



