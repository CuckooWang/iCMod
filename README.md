# iCMod
An efficient pipeline for computationally integrating time-course multi-omics data to identify normalized circadian phosphorylation sites (NCPs) that are oscillating in a circadian manner potentially due to rhythmic phosphorylation/dephosphorylation events. Potentially circadian kinases that phosphorylate and regulate these NCPs could also be inferred. 

## The description of each source code
### FPKM1.java
The code was used to remove the protein sequences whose corresponding FPKM values < 1.0 in all time points of each batch, and generate new reference protein sequence databases for the MaxQuant-based MS/MS searching of the proteomic and phosphoproteomic data.
### CorrectByWT3.java
WT_3 samples quantified in both Batch 1 and Batch 4 were used as the normalization control for the proteomic and phosphoproteomic data.
### PhosNormByPro.java
The intensities of p-sites were normalized by the intensities of their corresponding proteins.
### GetFasta.java
To generate the FASTA files for the prediction of site-specific kinase-substrate relations (ssKSRs) by GPS 2.1. 
### MetaCycleRun.R
The code calls MetaCycle, an integrative R package that incorporated three computational programs including ARSER, JTK_CYCLE and Lomb-Scargle, for the detection of circadian oscillations at different levels.
### CalCricadianSiteKa.java
Based on the predicted ssKSRs by GPS 2.1, the number of potential substrates of each protein kinase were calculated.
### enrichment.java
The hypergeometric test was adopted for the identification of potential circadian kinases.
### tools.java
Contains some routine methods that the above programs need to call for multiple times.
### test_data
This folder contains example data and files for testing above programs. It should be noted that some files are partially present due to the limitation of uploading size, such as "gene_exp.diff".

## Contact
Yu Xue: xueyu@hust.edu.cn  
Luoying Zhang: zhangluoying@hust.edu.cn  
Chenwei Wang: wangchenwei@hust.edu.cn  
Ke Shui: shuike@hust.edu.cn  
