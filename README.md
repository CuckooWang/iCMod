# iCMod
An efficient pipeline for computationally integrating circadian multi-omics data to acquire normalized cricadian p-sites that are oscillating in a circadian manner truly due to rhythmic phosphorylation/dephosphorylation events.

## Information of each code
### FPKM1.java
These codes were used to remove the protein sequence whose corresponding FPKM value < 1.0 in all time points of each batch, and generate the new reference protein sequence databases for the Maxquant search of Protomic and Phosphoproteomic data.
### CorrectByWT3.java
WT_3 sample quantified in both Batch 1 and Batch 4 were used as the normalization control for proteomic and phosphoproteomic data.
### PhosNormByPro.java
The intensities of p-sites were normalized by the intensities of their corresponding proteins quantified by proteome.
### GetFasta.java
To generate the fasta files for the predictions of GPS.
### MetaCycleRun.R
Circadian oscillations at different levels were identified by MetaCycle, an integrative R package that incorporated three computational programs including ARSER, JTK_CYCLE and Lomb-Scargle.
### CalCricadianSiteKa.java
Based on the predictions of GPS, the number of potential substrates of each protein kinase were calculated firstly.
### enrichment.java
The hypergeometric test was adopted for identification of potential circadian kinases.
### tools.java
Contains some methods that the above programs need to call.
### test_data
This folder contains the files that needed for the test of above programs, it is should noting that some files are incomplete due to the limitation of upload size, such as "gene_exp.diff".

## Contact
Yu Xue: xueyu@hust.edu.cn  
Luoying Zhang: zhangluoying@hust.edu.cn  
Chenwei Wang: wangchenwei@hust.edu.cn  
Ke Shui: shuike@hust.edu.cn  
