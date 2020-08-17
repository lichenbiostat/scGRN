# scGRN

A computational pipeline for predicting cell-type gene regulatory networks via multi-omics data integration and analyses. 

## Introduction

Identifying GRN (gene regulatory network) can help us better understand gene activities and biological processes. Nowadays, single-cell gene expression data have been an important resource for GRN computing, on account of reasons mentioned in [1]. So far a bunch of algorithms has been developed to decipher gene regulatory networks such as GENIE3, PIDC and GRN. Here we introduce **scGRN** - a computational pipeline to integrate single cell multi-omics data including chromatin interactions, epigenomics and transcriptomics, for predicting cell-type gene regulatory networks linking transcription factors, distal regulatory elements and target genes. For our project, we used data for the human brain (e.g., excitatory and inhibitory neurons, microglia, oligodendrocyte). The pipeline is depicted in the below flow chart.

![Pipeline](https://github.com/daifengwanglab/scGRN/blob/master/pipeline.png)


## Prerequisite

Our package needs to import the following R packages to work.

> GenomicInteractions, GenomicRanges, GenomicFeatures, GenomeInfoDb, IRanges, S4Vectors,
biomaRt, TFBSTools, glmnet, motifmatchr, data.table, dplyr, parallel, doParallel, foreach

Besides, data from packages *BSgenome.Hsapiens.UCSC.hg19*, *TxDb.Hsapiens.UCSC.hg19.knownGene* and *JASPAR2018* are also used for our project.

All of these packages can be installed either in Bioconductor or CRAN.

```R
# For example, I want to download R package - GenomicInteractions
# The following code would help
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicInteractions")
```

## Install scGRN R package
The devtools package provides install_github() that enables people to install packages from GitHub. You can use the following code to download the scGRN package.

```R
install.packages("devtools")
library(devtools)
install_github("daifengwanglab/scGRN")
```

After the installation, you can directly use **library(scGRN)** to load the package.

## Functions in the package

There are three functions inside the package.

1. **scGRN_interaction** :
    * main input:
        * hic_interaction: data frame containing the variables chr1,start1,end1,chr2,start2,end2
          or chr,start1,end1,start2,end2 
          chr1,ch2 or chr should have the following format 'chrD'(D represents digit 1-22 or X Y)                          
          start1,end1,start2,end2 should be integer type
          
        * enhancer: a data frame containing chr,start,end for enhancers
        
    * output: a data frame containing gene, gene_chr, promoter_start, promoter_end,               
      enh_chr,enh_start,enh_end
      
The function, scGRN_interaction inputs the chromatin interaction data (e.g., Hi-C) and predicts all possible interactions between enhancers and promoters in the data or the user-provided list; i.e., ones from Topologically Associating Domains (TADs) in Hi-C data. In addition, the function uses an R package, GenomicInteractions [4] to annotate interacting regions and link them to genes.
 
2. **scGRN_getTF** :
    * main input: 
        * df: output data frame from **scGRN_interaction**
    
    * output: a data.table which contains transcription factors for each region
    
Function scGRN_getTF infers the transcription factor binding sites (TFBS) based on consensus binding site sequences in the enhancers and promoters that potentially interact from the previous step, scGRN_interaction. It outputs a reference gene regulatory network linking these TF, enhancers and/or promoters of genes. In particular, this function uses TFBSTools [5] to obtain the position weight matrices of the TFBS motifs from the JASPAR database [6] and predict the TFBS locations on the enhancers and promoters via mapping TF motifs. It further links TFs with binding sites on all possible interacting enhancers and promoters, and outputs the reference regulatory network. Furthermore, this function can run on a parallel computing version via an R package, motifmatchr [7] for computational speed-up.

3. **scGRN_getNt** :
    * main input: 
        * df: output data.table from **scGRN_getTF**
        * gexpr: gene expression data in which each row represents a gene and each column represents an observation
        
    * output:  a data frame containing TG, TF, promoter, enhancer and coef
    
The function, scGRN_getNt predicts the final gene regulatory network based on the TF-target gene expression relationships in the reference network. The reference gene regulatory network from the previous step provides all possible regulatory relationships (wires) between TF, enhancers, and target genes. However, changes in gene expression may trigger different regulatory wires. To refine our maps and determine the activity status of regulatory wires, this function applies elastic net regression, a machine learning method that has successfully modelled gene regulatory networks in our previous work [3].  In particular, given a gene expression dataset and a reference network from scGRN_getTF, the function uses the TF expression to predict each target gene expression, and finds the TF with high regression coefficients, indicating an active regulatory influence on the target gene’s expression in the gene expression data. The final gene regulatory network consists of the TF with high elastic net coefficients, target genes and the linked enhancers from their reference network links if any. 
    
## Example

To present how our pipeline works, we will use microglia interactome data [2], microglia enhancers [2], gene expression data [3] to illustrate the steps.

To begin with, we prepare and format the data in R.

```{r}
library(readxl)
interactome_data = read_xlsx("PLAC-seq promoter interactome map.xlsx",
                       sheet = 'Microglia interactome',skip = 2)[,1:6]

enhancers = read_xlsx("PLAC-seq promoter interactome map.xlsx",
                       sheet = 'Microglia enhancers',skip = 2,
                       col_names = c('chr','start','end'))
    
gexpr = read.table('DER-22_Single_cell_expression_raw_UMI.tsv',
                       header = T,row.names = 1,sep = '\t')
# Here we do some simple data processing without normalizing the data
# Select microglia cells 
# Remove genes that are not expressed (based on the gene expression data)
# For more details about normalizing single cell gene expression data, users can check Seurat 3.0 [8].
gexpr = gexpr[,grep('Micro',colnames(gexpr))]
gexpr = gexpr[rowSums(gexpr != 0) > 5,]
```

Now we will start using the functions in our package.

**Step1: find interaction**

```{r}
df1 <- scGRN_interaction(interactome_data,enhancers)
# For convenience, I use a subset of df1 and show relevant outputs.
df1 <- df1[sample(nrow(df1),20),]
head(df1)
```

    ##         gene gene_chr promoter_start promoter_end enh_chr enh_start   enh_end
    ## 70681  PICK1    chr22       38451253     38456253   chr22  38731974  38732426
    ## 59812  CBLL1     chr7      107381779    107386779    chr7 107515988 107516634
    ## 25923  GSK3A    chr19       42743546     42748546   chr19  42629613  42633317
    ## 4621    CHD3    chr17        7810319      7815319   chr17   7942326   7946725
    ## 62217 PICALM    chr11       85689772     85694772   chr11  85859249  85859502
    ## 28046  RAB43     chr3      128838148    128843148    chr3 128952814 128954805



**Step2: get TFs for each promoter and enhancer**

```{r}
df2 <- scGRN_getTF(df1)
head(df2,1)
```
   
    ##     gene                promoter                enhancer
    ## 1: PICK1 chr22:38451253-38456253 chr22:38731974-38732426
    ##                             promoter_TF                          enhancer_TF
    ## 1: IRF2,PPARG,RXRA,RREB1,NR1H2,REST,... CREB1,MZF1,PPARG,RXRA,RORA,NR1H2,...
 
**Step3: predict TF -> TG via Elastic net**

Last, what we need to do is to input the gene expression data and the output of function *scGRN_getTF*.

```{r}
df3 <- scGRN_getNt(df = df2, gexpr = gexpr)
head(df3)
```
    ##      TG    TF                enhancer                promoter     TFbs          coef
    ## 1 PICK1 CREB1 chr22:38731974-38732426 chr22:38451253-38456253 enhancer   0.004038054
    ## 2 PICK1  MZF1 chr22:38731974-38732426 chr22:38451253-38456253 enhancer   0.016943537
    ## 3 PICK1 PPARG chr22:38731974-38732426 chr22:38451253-38456253     both  -0.008975126
    ## 4 PICK1  RXRA chr22:38731974-38732426 chr22:38451253-38456253     both  -0.002222325
    ## 6 PICK1 NR1H2 chr22:38731974-38732426 chr22:38451253-38456253     both   0.002630548
    ## 7 PICK1   SP1 chr22:38731974-38732426 chr22:38451253-38456253 enhancer  -0.002453652


## Reference and source
1. Pratapa A, Jalihal AP, Law JN, Bharadwaj A, Murali TM. Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data. Nat Methods. 2020;17:147–54. 

2. Nott A, Holtman IR, Coufal NG, Schlachetzki JCM, Yu M, Hu R, et al. Brain cell type-specific enhancer-promoter interactome maps and disease-risk association. Science. 2019;366:1134–9. 

3. Wang D, Liu S, Warrell J, Won H, Shi X, Navarro FCP, et al. Comprehensive functional genomic resource and integrative model for the human brain. Science. 2018;362. 

4. Harmston, N., Ing-Simmons, E., Perry, M., Baresic, A., Lenhard, B. GenomicInteractions: R package for handling genomic interaction data [Internet]. 2020. Available from: https://github.com/ComputationalRegulatoryGenomicsICL/GenomicInteractions/

5. Tan G, Lenhard B. TFBSTools: an R/bioconductor package for transcription factor binding site analysis. Bioinforma Oxf Engl. 2016;32:1555–6. 

6. Fornes O, Castro-Mondragon JA, Khan A, van der Lee R, Zhang X, Richmond PA, et al. JASPAR 2020: update of the open-access database of transcription factor binding profiles. Nucleic Acids Res. 2020;48:D87–92. 

7. Schep, Alicia. motifmatchr: Fast Motif Matching in R [Internet]. 2019. Available from: https://www.bioconductor.org/packages/release/bioc/html/motifmatchr.html

8. Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM, et al. Comprehensive Integration of Single-Cell Data. Cell. 2019;177:1888-1902.e21. 
