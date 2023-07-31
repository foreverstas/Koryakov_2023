#SalmonTE output EXPR.csv DataFrame was used for differential expression analysis using DESeq2 algorithm

#import R libraries
library(DESeq2)
library(tidyverse)
library(IHW)
library(ashr)
library(ggplot2)

#importing SalmonTE pseudocounts
data<-read.csv("SalmonTE_output/EXPR.csv",sep = ",",header=T)%>%column_to_rownames("TE")%>%mutate_if(is.numeric, round)
colnames(data)<-c("SG_dm3","SG_dm2", "SG_dm1","SG_OR3","SG_OR2", "ID_dm3","ID_dm2","ID_dm1", "ID_OR3","ID_OR2","ID_OR1", "SG_OR1")

meta<-data.frame(Sample=colnames(data),Tissue=c(rep("SG",5),rep("ID",6),"SG"),Type=c(rep("dm",3),rep("OR",2),rep("dm",3),rep("OR",4)))

meta$Tissue<-factor(meta$Tissue)
meta$Type<-factor(meta$Type)

#Convert data dataframe to DESeqDataSet object
dds<-DESeqDataSetFromMatrix(data,meta,design=~Tissue+Type+Tissue*Type)
dds<-DESeq(dds)

#Modeling contrasts
modMat<-model.matrix(design(dds), colData(dds))
SG_wt<-colMeans(modMat[dds$Tissue=="SG" & dds$Type=="OR",])
SG_DKO<-colMeans(modMat[dds$Tissue=="SG" & dds$Type=="dm",])
ID_wt<-colMeans(modMat[dds$Tissue=="ID" & dds$Type=="OR",])
ID_DKO<-colMeans(modMat[dds$Tissue=="ID" & dds$Type=="dm",])

#Calculate differental expression using contrasts
lfcThreshold=1
alpha=0.05
altHypothesis="greaterAbs"
SG_DE_TE<-results(dds,contrast=(SG_DKO-SG_wt),alpha=alpha,lfcThreshold=lfcThreshold,altHypothesis=altHypothesis,filterFun = ihw)
ID_DE_TE<-results(dds,contrast=(ID_DKO-ID_wt),alpha=alpha,lfcThreshold=lfcThreshold,altHypothesis=altHypothesis,filterFun = ihw)

#Shrink LFC estimates
SG_DE_TE_shr <- lfcShrink(dds, contrast = SG_DKO-SG_wt, res=SG_DE_TE, type="ashr",lfcThreshold = 1)
ID_DE_TE_shr <- lfcShrink(dds, contrast = ID_DKO-ID_wt, res=ID_DE_TE, type="ashr",lfcThreshold = 1)

#Convert DESeq2 results to single data.frame

LFC_table_DKO<-SG_DE_TE_shr%>%data.frame(.)%>%dplyr::select(SG=log2FoldChange,SGp=padj)%>%rownames_to_column("row")
LFC_table_DKO<-ID_DE_TE_shr%>%data.frame(.)%>%dplyr::select(ID=log2FoldChange,IDp=padj)%>%rownames_to_column("row")%>%left_join(LFC_table_DKO)

LFC_table_DKO$color="gray"
LFC_table_DKO[which(LFC_table_DKO$IDp<0.05&LFC_table_DKO$SGp<0.05),"color"]<-"red"
LFC_table_DKO[which(!(LFC_table_DKO$IDp<0.05)&LFC_table_DKO$SGp<0.05),"color"]<-"green"
LFC_table_DKO[which(LFC_table_DKO$IDp<0.05&!(LFC_table_DKO$SGp<0.05)),"color"]<-"blue"
LFC_table_DKO$color<-factor(LFC_table_DKO$color,levels = c("gray","blue","green","red"))
LFC_table_DKO<-LFC_table_DKO%>%arrange(color)

#Plot row image for Supplementary Figure 2B
theme_set(theme_bw())

p<-ggplot(LFC_table_DKO%>%arrange(color), aes(x=SG, y=ID)) + 
  geom_point(aes(color=color),size=0.5) +
  scale_color_manual(values=c("gray","blue","green","red"),
                     name=element_blank(), 
                     labels=c("Not Affected TE","ID DE TE","SG DE TE","Common DE TE")) + 
  coord_fixed()+
  ylim(-7,7)+
  xlim(-7,7)+ 
  xlab("LFC SG") + 
  ylab("LFC ID")
ggsave("TE_LFC_comparison_ID_vs_SG.pdf",device = cairo_pdf,plot = p)

#sessionInfo()
#R version 4.1.2 (2021-11-01)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Debian GNU/Linux bookworm/sid

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.19.so

#locale:
 #[1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 #[4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 #[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
#[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
# [1] ashr_2.2-54                 pheatmap_1.0.12             ggrastr_1.0.1              
# [4] ggpubr_0.4.0                genomation_1.24.0           GenomicFeatures_1.44.2     
# [7] AnnotationDbi_1.56.2        lubridate_1.9.2             forcats_1.0.0              
# [10] stringr_1.5.0               dplyr_1.1.2                 purrr_1.0.1                
# [13] readr_2.1.4                 tidyr_1.3.0                 tibble_3.2.1               
# [16] ggplot2_3.4.2               tidyverse_2.0.0             DESeq2_1.34.0              
# [19] SummarizedExperiment_1.24.0 Biobase_2.54.0              MatrixGenerics_1.6.0       
# [22] matrixStats_0.63.0          GenomicRanges_1.46.1        GenomeInfoDb_1.30.1        
# [25] IRanges_2.28.0              S4Vectors_0.32.3            BiocGenerics_0.40.0        
# [28] IHW_1.20.0                 

#loaded via a namespace (and not attached):
#  [1] backports_1.4.1          circlize_0.4.15          BiocFileCache_2.0.0      plyr_1.8.8              
# [5] splines_4.1.2            BiocParallel_1.28.3      gridBase_0.4-7           lpsymphony_1.20.0       
# [9] digest_0.6.31            invgamma_1.1             foreach_1.5.2            htmltools_0.5.5         
# [13] SQUAREM_2021.1           fansi_1.0.4              magrittr_2.0.3           memoise_2.0.1           
# [17] BSgenome_1.60.0          cluster_2.1.2            doParallel_1.0.17        tzdb_0.4.0              
# [21] ComplexHeatmap_2.10.0    Biostrings_2.62.0        annotate_1.72.0          timechange_0.2.0        
# [25] prettyunits_1.1.1        colorspace_2.1-0         blob_1.2.4               rappdirs_0.3.3          
# [29] countToFPKM_1.0          xfun_0.39                crayon_1.5.2             RCurl_1.98-1.12         
# [33] genefilter_1.76.0        impute_1.66.0            survival_3.2-13          iterators_1.0.14        
# [37] glue_1.6.2               gtable_0.3.3             zlibbioc_1.40.0          XVector_0.34.0          
# [41] GetoptLong_1.0.5         DelayedArray_0.20.0      car_3.0-12               shape_1.4.6             
# [45] abind_1.4-5              scales_1.2.1             DBI_1.1.3                rstatix_0.7.0           
# [49] Rcpp_1.0.10              plotrix_3.8-2            xtable_1.8-4             progress_1.2.2          
# [53] clue_0.3-64              bit_4.0.4                DT_0.20                  truncnorm_1.0-9         
# [57] htmlwidgets_1.6.2        httr_1.4.6               RColorBrewer_1.1-3       pkgconfig_2.0.3         
# [61] XML_3.99-0.14            dbplyr_2.3.2             locfit_1.5-9.4           utf8_1.2.2              
# [65] tidyselect_1.2.0         rlang_1.1.1              reshape2_1.4.4           munsell_0.5.0           
# [69] tools_4.1.2              cachem_1.0.6             cli_3.6.1                generics_0.1.3          
# [73] RSQLite_2.2.9            broom_1.0.4              fdrtool_1.2.17           evaluate_0.21           
# [77] fastmap_1.1.0            yaml_2.3.7               knitr_1.42               bit64_4.0.5             
# [81] KEGGREST_1.34.0          slam_0.1-50              xml2_1.3.3               biomaRt_2.48.3          
# [85] compiler_4.1.2           rstudioapi_0.14          beeswarm_0.4.0           filelock_1.0.2          
# [89] curl_5.0.0               png_0.1-8                ggsignif_0.6.3           geneplotter_1.72.0      
# [93] stringi_1.7.12           lattice_0.20-45          Matrix_1.5-4.1           vctrs_0.6.2             
# [97] pillar_1.9.0             lifecycle_1.0.3          GlobalOptions_0.1.2      irlba_2.3.5.1           
#[101] data.table_1.14.8        bitops_1.0-7             rtracklayer_1.52.1       R6_2.5.1                
#[105] BiocIO_1.2.0             KernSmooth_2.23-20       vipor_0.4.5              codetools_0.2-18        
#[109] rjson_0.2.21             withr_2.5.0              GenomicAlignments_1.28.0 Rsamtools_2.8.0         
#[113] GenomeInfoDbData_1.2.7   parallel_4.1.2           hms_1.1.3                rmarkdown_2.21          
#[117] carData_3.0-4            mixsqp_0.3-48            seqPattern_1.24.0        ggbeeswarm_0.6.0        
#[121] restfulr_0.0.13
