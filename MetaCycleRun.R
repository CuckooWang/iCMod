library(MetaCycle)

meta2d(infile="KAA_wt_16.txt", filestyle="txt", 
       outdir="WT_out_16", timepoints=rep(seq(0, 45, by=3), each=1),
       minper = 20,cycMethod=c("ARS","JTK","LS"),analysisStrategy = "selfUSE",
       outIntegration="noIntegration")
