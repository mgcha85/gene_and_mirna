# prediction target genes of miRNA using Lasso


Table of Contents
=================

  * [Data preparation](#Data-preparation)
  * [Pre-processing](#Pre-processing)
  * [Comparison](#Comparison)
  * [Methodology](#Methodology)
  * [Validation](#Validation)
  * [Summary](#Summary)
  
  
  

Data-preparation
============
## Data-resources


|Resource | URL|
| ------------- | ------------- |
|Annotated genes | https://www.gencodegenes.org/human/release_32lift37.html|
|RNA-seq by tissues | https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1733/samples/?s_page=1&s_pagesize=500|
|FANTOM by tissues | https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz|

FANTOM CAGE peaks expression tables for the robust (expression>10TPM) set of peaks.  
It based expression table (RLE normalized) for human with annotation.  
RLE normalization is explained on [this paper](https://www.nature.com/articles/sdata2017112).  

## Converting Data Format

| Convert | Command |
| ------------- | ------------- |
| fastq → sam (hi-sat) | hisat2 -p 4 -x genome_tran -1 {fastq1} – 2 {fastq2} -S {sam} |
| sam → bam (samtools) | samtools sorted -@ 8 -o {bam} {sam} |
| bam → gtf (stringtie)	| stringtie -p 4 -G genes.gtf -o {gtf} -i {bam} |

## File entries by Tissues

|Resource | URL|
| ------------- | ------------- |
|FANTOM | https://drive.google.com/open?id=1y3w-rOvgbxtlCKWss1wAzcHD884Na__tijNWL7Yjb2I|
|RNA-seq | https://drive.google.com/open?id=1lKZN5hn2e6Zq3eKzn5r3H5NqdEfjgUmy|

## 22 common tissues between FANTOM and RNA-seq
The common tissues are listed in the below table. 


| RNA-seq | FANTOM |
| ------------- | ------------- |
| appendix | appendix |
| urinarybladder | bladder |
| brain	| brain|
| colon	| colon|
| esophagus | esophagus|
| fat | adipose|
| gallbladder | gall_bladder|
| heart | heart|
| kidney | kidney|
| liver | liver|
| lung | lung|
| lymphnode | lymph_node|
| ovary | ovary|
| pancreas | pancreas|
| placenta | placenta|
| prostate | prostate|
| salivarygland | salivary_gland|
| smallintestine | small_intestine|
| spleen | spleen|
| testis | testis|
| thyroid | thyroid|
| endometrium | uterus|


Pre-processing
============
## Averaged data by tissue
One tissue can have multiple replicates. Therefore, RNA-seq and FANTOM were averaged on same tissues. For this, we  only used chr1, chr2, …chr 22, chrX and chrY. The average is calculated by each chromosome and strand.

## Expression data by 22 tissues
RNA-seq data already contains FPKM while FANTOM does not have expression data. To get the expression data from FANTOM, we calculated expression level of (+/-)100 bp from annotated gene TSS. The expression level is summed by the score sum of cage tags within the region. When the summation, every tags, which are overlapped to the region;  [TSS-100, TSS+100] should be used.

## CAGE raw data vs. CAGE tag processed data
### Before  vs. After Filter

There are three filters.
1. maximum value of FANTOM or RNA-seq on 22 tissues is zero.  
2. median of FANTOM or RNA-seq on 22 tissues is zero.  
3. maximum value is greater than sum of others.  

### after filter
After the filtering, the total number of transcripts for 27,493 out of 84,088, which is the number of TSS in gencode.


## High consistent genes
RNA-seq and FANTOM data have expression level by a gene on the 22 tissues. Each gene has two 22 length vectors (RNA-seq and FANTOM), which an element show expression level on a tissue. By the two vectors, the correlation coefficient is calculated to observe how consistent between RNA-seq and FANTOM on a gene. Each gene has the coefficient and high coefficient shows two data source are consistent on the 22 tissues.  
For the correlation, **spearman method** was used. We round the coefficient at the decimal point with two digit.  
To extract only high correlated one, we set **threshold as 0.75**. Finally, **5,116 transcripts** are received as high consistent ones from **27,493** transcripts. the transcripts are grouped by gene name and if there are multiple transcripts, we picked one with the maximum score. After this, we got **3,738 genes**.  
The list of 3,738 transcripts is in the below link. only double type is used for every calculation. Transcript region are only used that gene and transcript type are protein coding gene from genecode transcripts.  
[link to see](https://drive.google.com/file/d/1Q9PvJdm1jVIW7zFH2ZE00rXkXMqLzZm4/view?usp=sharing)


## High consistent miRNA
For the consistent miRNAs, we compared eleven papers and retrieved the TSS locations by the miRNA name. The consistent miRNA TSSs are predicted by all papers and the TSS regions are close then we chose the meddle point as the TSS. The number of miRNA TSSs is 330. The name and TSS of 330 miRNAs are in the below link.  
[link to see](https://drive.google.com/open?id=1qZHBqubcYeJfVk7uvE8lFyzVXZeSknmIOf-qfS5zBmM)


Methodology
============
### FANTOM cell line data
we received cell line specific data. The below table shows the url of the data source and the table including file id, cell line and so on.

| label | URL|
| ------------- | ------------- |
| FANTOM url | http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.tissue.hCAGE/|
| cell line data list | https://drive.google.com/open?id=1frNwIbWzFwWdvKfThP8DT1mdbPfGRERy|

### Expression data by 240 cell lines
From the data preparation step, we got high consistent transcripts on gene and miRNA. We extracted the expression level from 240 cell lines. The score sum of the region where (+/-)100bp from annotated gene TSS is the expression level. We calculated this levels about the high consistent genes and miRNAs. 
some of 240 data does not have chromosomeY so, we did not consider the chromosome. Therefore, the total number of genes and miRNA is **2,313** and **330** as we mentioned. 

| label | URL|
| ------------- | ------------- |
| miRNA expression data url | https://drive.google.com/file/d/13pG6rZnK-Wn7Tn1cPrQe67-piAowJ5Jl/view?usp=sharing|
| gene expression data url | https://drive.google.com/file/d/12CPzPnH-d8bRTFvAtL-mUlsN6mQDnIJV/view?usp=sharing|


### Lasso linear regression
To predict target genes on a miRNA, we used lasso regression. the consistent genes and miRNAs have 240 vector, which are expression level on each cell line. Assume the number of genes and miRNA are m and n respectively. gene matrix **X**: (m x 240) and miRNA matrix **Y**: (n x 240).  
Lasso calculated the relationship between **X** and **Y**. The coefficient matrix by the Lasso result is **B** (n x m); **Y = BX**.  
Before, processing Lasso, **X** and **Y** were centered by subtracting mean value by each row.  

From the coefficient matrix (**B**), each row shows coefficient of one miRNA and multiple genes. Therefore, the non-zero coefficient values are assumed target genes, which are related to the corresponding miRNA. The coefficient of Lasso is [here](https://drive.google.com/file/d/1Y2sA1EE5KC4ZGyXqjyQnHA877p7NfC6O/view?usp=sharing)  

Figure 1 shows the mean of **B** as a sample.
![Image](/images/Figure_1.png)
Figure1. The mean of **B**    

When the coefficient is not zero between gene and miRNA From the result, we assume those are target genes by a miRNA.   
The target genes are [here](https://drive.google.com/file/d/1gQlfWrGh_GZEGHGTyeiVxajql9d7dSKL/view?usp=sharing).  
257 miRNA has target genes out of 330. This means 73 miRNAs has zero coefficients becasue they failed to Lasso convergence.  

### Gene ontology analysis
From the lasso result, we got the target genes by a miRNA. To investigate the reliability of the result, we processed gene ontology analysis using the result. For the GO analysis, the high consistent genes are used as background genes. Target genes is the lasso result.  For the analysis, we use the below website.  
[link to see](http://cbl-gorilla.cs.technion.ac.il/)
From this analysis, we recevied another gene set based on the Lasso result.  
the result is [here](https://drive.google.com/file/d/1JXhQbQLv6scJV88k3499GNvHKP8zWBtk/view?usp=sharing).  

Figure 2 shows the input of GO by the website.
![Image](/images/Figure_2.png)  
Figure2. GO input set

![Image](/images/Figure_3.png)
Figure3. GO output  

### Gene Set Enrichment Analysis
For GSEA analysis, three inputs are required;  

1. Expression dataset  

2. phenotype labels  

3. gene sets  

For the expression dataset, we used the table which consist of gene rows and cell line columns. The gene names
of the dataset are required only Hugo symbol although GECODE contains
multi-symbols. Thus, the symbols were converted to Hugo symbols by
HGNC.  [https://www.genenames.org/tools/multi-symbol-checker/](https://www.genenames.org/tools/multi-symbol-checker/) 

For phenotype label, we selected
continuous file format as peak profile. The peak values are obtained by
averaging expression levels by cell lines and converted to integer value.  

Finally, molecular signature database (MSigDb) is used for gene set. We
had target genes by using Lasso. Through GSEA,
we searched statistically important gene set.  

In this experiment, we used some specific parameters. Since the expres-
sion dataset is already gene level data, we did not collapse dataset to gene
symbols. GSEA ranks the genes and analyzes ranked list of genes. For the
rank, various options exist. In this experiment, we used continuous phe-
notypes and Pearson, Cosine, Manhattan and Euclidean are available for
the phenotypes. For this model, Manhattan is proper ranking method since
it calculates similarity through distance with only one direction. The
GSEA software is v3.0 from broadinstitute.  

The GSEA analysis result are in the below link.  

[link to see](https://drive.google.com/open?id=1SfNwJtYHWc1oRLbuy4lcu8Fd9cTsgpJj)  

Figure 4. shows input data of GSEA and Figure 5 shows the parameters of GSEA.  
![Image](/images/Figure_4.png)
Figure4. GSEA input profile  

![Image](/images/Figure_5.png)
Figure5. GSEA parameter set  

![Image](/images/Figure_6.png)
Figure6. GSEA heat map  


Validation
============
#### Lasso 
10 fold cross-validation  

For the cross-validation, we split the 240 vector to 216 (90%) as a train and 24 (10%) as a test. by the two set, we calculated the distance between actual Y and predicted Y. The below equation shows the distance.  

![Image](/images/formula1.gif)  
  
![Image](/images/formula2.gif)  

where Yactual is miRNA expression level, X is gene expression level and N is the size of Y matrix (n x 216 or n x 24). The below table shows 10-fold cross validation result.
(n: 330, N: 71,280 or 79,20)  


| # | test | train |
| ------------- | ------------- | ------------- |
| 0 | 0.0043 | 0.0553 |
| 1 | 0.0026 | 0.0565 |
| 2 | 0.0022 | 0.0580 |
| 3 | 0.0021 | 0.0578 |
| 4 | 0.0026 | 0.0572 |
| 5 | 0.0025 | 0.0570 |
| 6 | 0.0026 | 0.0561 |
| 7 | 0.0026 | 0.0565 |
| 8 | 0.0021 | 0.0572 |
| 9 | 0.0027 | 0.0562 |


### FANTOM cell lines vs. RNA-seq tissues
To validate 240 cell lines data, we checked out RNA-seq data as well. From EMBL-EBI, we have RNA-seq data on 27 tissues. The number of high consistent transcripts of gene and miRNA are m and n respectively. Then, X: (m x 27) and Y: (n x 27). The RNA-seq data contains only 264 miRNAs out of 330 miRNAs. By using this **X and Y**, we calculated coefficient matrix **B** by Lasso.  

To compare FANTOM and RNA-seq data, we calculated three distances as the below formula.

distance1: ![Image](/images/formula1-1.gif) and ![Image](/images/formula1-2.gif)  
distance2: ![Image](/images/formula2-1.gif) and ![Image](/images/formula2-2.gif)  

where B_{CELL} is calculated by using the 264 miRNAs’ expression data on FANTOM 240 cell lines, X_{RNA} is RNA-seq gene expression matrix with size (m x 27), Y_{RNA} is RNA-seq miRNA expression matrix with size (n x 27).  

 

### Lasso result vs. other researches
we use three researches; MIRANDA, RNA22 and TargerScan for the comparison. From the three researches, we extracted target genes by a miRNA. One version is the intersection of three researches and another version is the union of the researches. For example, miRNA1 has gene1, gene2 and gene3 by research1, gene2, gene4 and gene5 by research2 and gene2, gene6 and gene7 by research3. Then the intersection of the researches is gene2 and the union of the researches is gene1, gene2, gene3, gene4, gene5, gene6 and gene7. We applied GSEA and GO on this result.

We checked the correlation of expression level between RNA-seq and FANTOM tag using target genes of other researches. The union of target genes are used for this. The expression level are obtained as we have done in *"Expression data by 240 cell lines"*. 

We checked the Lasso regression if the target genes have many non-zeros values from the coefficient. Based on the result, only 4.06% target genes has non-zeros coefficient from intersection result and 3.55% target genes has non-zeros coefficient from union result

	
| | % non-zeros | mean of correlation |
| ------------- | ------------- | ------------- |
| Intersection | 4.06% | 0.02 |
| Union | 3.55% | 0.02 |


### Hypertest
we calculated hypergeometic test by comparing this result to other research or other result such as RNA-seq or GO. Since each miRNA has target genes, we can compare the target genes based on each miRNA. For this, we use survival function to calculate p-value. The function has four inputs, which q, m, n and k; q is # of the intersection of predicted genes by Lasso and predicted genes by another research. m is # of genes by another research, n is the # of high correlated genes – m and k is # of target genes by Lasso result. For the survival function, q-1 is used instead of q because the function is the inverse of the cumulative distribution function (1-cdf). Therefore the survival function excludes till q-1 instead of q.  

Not only we calculated p-value between this result and another result but we also consider the common genes by all three researches by a miRNA and common genes, which are overlapped by at least two methods and the union of the target gene are calculated. Also, GO has q-value by a miRNA so it is used for this comparison.  

The below link contains the all comparison for the hypergeometric test.  
[link to see](https://drive.google.com/file/d/1CA7FK5FSHs_Gd3LjONZLau5QnFNSILx_/view?usp=sharing)  


q significant is the multiplication of q-value (GO) and row number.  

The # significant column shows statistically importance of miRNA. The number shows how many research is important. The importance is considered if the p-value is smaller than 1 / # miRNA; threshold. Therefore, 2 means p-values of two researches are smaller than threshold and important.  
The result shows only 13 miRNAs are statistically important out of 327.  

Summary
=========
For data preparation,  

I used raw data from FANTOM and processed data from FANTOM. I calculated correlation  coefficient by a transcript by each data set. The correlation of the correlation coefficient between them is 0.6866. And I also check random generated value, which has same min and max value as raw data. The random generated values for miRNA and gene are computed for correlation and I compared this to correlation from raw data. The result is 0.00229, which means they are almost orthogonal.  
For the validation, 
1. 10 cross validation by using lasso result and calculate distance
2. lasso result by 240 cell lines and lasso result by 22 tissues comparison: X_cell Y_tissue = Y_pred is compared to Y_cell 
3. Lasso result vs. other researches: calculated correlation coefficient of expression level between miRNA and target genes and lasso regression to check how many non-zero (target genes for this research) values exist.  
4. GSEA and GO are calculated to evaluate this result.
