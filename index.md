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




## Welcome to GitHub Pages

You can use the [editor on GitHub](https://github.com/mgcha85/gene_and_mirna/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/mgcha85/gene_and_mirna/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
