# Example PLSR Analysis for AHBA gene expression vs. cortical thickness (CT) deviance for 22q11.2 Deletion Syndrome

#Reading and writing tables
library(openxlsx)

#Data wrangling
library(plyr)
library(tidyverse)
library(purrr)

#Statistical Analyses
library(pls)
library(caret)

setwd("/home/jovyan/curriculum/forsyth-genetics")
data_folder = "/home/jovyan/shared/forsyth-genetics/"

###########################################################################################################################################################################
################################################## PLSR for 22q1Del CT Deviance with AHBA Gene Expression Data ############################################################

#Load AHBA data in left hemisphere DK atlas parcellation for protein-coding, brain expressed genes
#Adapted from French & Paus, 2015 (https://doi.org/10.3389%2Ffnins.2015.00323) 
geneexp.lh.protcoding.brainexprs.zscore.mat <- readRDS(paste0(data_folder, "geneexp.lh.protcoding.brainexprs.zscore.mat.rda"))

#Look at structure of AHBA gene expression data - how many genes + brain structures do you see?
dim(geneexp.lh.protcoding.brainexprs.zscore.mat)
str(geneexp.lh.protcoding.brainexprs.zscore.mat)
rownames(geneexp.lh.protcoding.brainexprs.zscore.mat)

#Load 22q11.2 deletion subjects versus healthy control map of z-score group effect sizes across left hemisphere DK atlas regions for cortical thickness
#From Supplementary Table 3 of Forsyth et al., 2021 (https://doi.org/10.1093/cercor/bhab008) 
ENIGMA22q.CT.deviance <- readRDS(paste0(data_folder, "ENIGMA22q_DelvsHC_ZScoreDeviance_CT.rda"))

#Look at structure of 22qDel deviance scores - how many brain structures do you see?
str(ENIGMA22q.CT.deviance)

#Make data frame for plsr with variable for 22qDel vs HC CT effect sizes across regions and matrix with AHBA gene expression across same regions as a variable
ENIGMA22q.CT.deviance.geneexp.zscore.reformat <- data.frame(ENIGMA.Dev.DelvsHC.ZScore.CT = ENIGMA22q.CT.deviance, AHBA.geneexp = I(as.matrix(geneexp.lh.protcoding.brainexprs.zscore.mat)))

#Alternatively, can load in already made dataframe 
#ENIGMA22q.CT.deviance.geneexp.zscore.reformat <- readRDS(paste0(data_folder, "ENIGMA22q.CT.deviance.geneexp.zscore.reformat.rda"))

#Run partial least square regression for 22qDel vs HC effect size map and AHBA gene expression data
pls.AHBA.geneexp.DelvsHC.ZScore.CT.model = plsr(ENIGMA.Dev.DelvsHC.ZScore.CT ~ AHBA.geneexp, data = ENIGMA22q.CT.deviance.geneexp.zscore.reformat, validation="LOO")

#Check LOO CV-based RMSEP output for models with different numbers of principal components
summary(pls.AHBA.geneexp.DelvsHC.ZScore.CT.model)

#RMSEP plot for CT Deviance
pdf("plsr.22qDelvsHC.ZScore.CT.numcomponents.byRMSEP.pdf")
validationplot(pls.AHBA.geneexp.DelvsHC.ZScore.CT.model, val.type="RMSEP", legendpos = "bottomright") #create plot
dev.off() #Close the pdf file

#R2 plot for CT Deviance
pdf("plsr.22qDelvsHC.ZScore.CT.numcomponents.byR2.pdf")
plot(pls.AHBA.geneexp.DelvsHC.ZScore.CT.model, "validation", val.type="R2") #create plot
dev.off() #Close the pdf file 

# Find the number of dimensions/components with lowest cross validation error
CV.CTDev = RMSEP(pls.AHBA.geneexp.DelvsHC.ZScore.CT.model)
CTDev.best.dims = which.min(CV.CTDev$val[estimate = "adjCV", , ]) - 1

# Re-run plsr specifying 3 components
pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel = plsr(ENIGMA.Dev.DelvsHC.ZScore.CT ~ AHBA.geneexp, data = ENIGMA22q.CT.deviance.geneexp.zscore.reformat, ncomp = 3)

#Compute Observed RMSEP & R2 from Final model 
pls.pred.Del_HC.CT.zscore = predict(pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel, ENIGMA22q.CT.deviance.geneexp.zscore.reformat, ncomp=1)
pls.pred.Del_HC.CT.zscore.eval = data.frame(obs=ENIGMA22q.CT.deviance.geneexp.zscore.reformat$ENIGMA.Dev.DelvsHC.ZScore.CT, pred=pls.pred.Del_HC.CT.zscore[,1,1])
pls.pred.Del_HC.CT.zscore.eval.summary <- t(defaultSummary(pls.pred.Del_HC.CT.zscore.eval)) %>% as.data.frame()

#Extract gene loadings for all components output from plsr
pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel.Yscores <- Yscores(pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel)

#Extract gene loadings for 1st component output from plsr and format into percentiles
pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights <- data.frame(Loading.Weight = loading.weights(pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel)[,1])
pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights$Gene <- colnames(ENIGMA22q.CT.deviance.geneexp.zscore.reformat$AHBA.geneexp)
pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights$Loading.Weight.round <- round(pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights$Loading.Weight, 5)
pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights$Loading.Weight.perrank <- round(percent_rank(pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights$Loading.Weight.round),5)

####Note: Common downstream functional annotation of top PC1 loading genes includes gene ontology (GO) term enrichment, such as through g:profiler (https://biit.cs.ut.ee/gprofiler/gost) or similar tools
#Can also test enrichment of top PC1 load genes for cell-type specific expression, such as through pSI tool (http://doughertytools.wustl.edu/CSEAtool.html, https://sites.wustl.edu/doughertylab/psi_package-page/)

#Save output for all brain-expressed, protein-coding AHBA gene loading weights onto first PC for 22qDel deviance
write.xlsx(pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights, "plsr.AHBA.DelvsHC.ZScore.CT.loadingweights.xlsx")


#Read in list of Gene symbols of 22q genes
GeneList22q <- read.xlsx(paste0(data_folder, "genelist22q.MarshallSandersUnion.xlsx"), colNames = TRUE)

#Extract percentile loadings onto 1st PC of 22q genes
pls.22q.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights <- pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights[pls.AHBA.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights$Gene %in% GeneList22q$HGNCUpdatedSymbol, ]

#Save output for brain-expressed, protein-coding 22q gene loading weights onto first PC for 22qDel deviance
write.xlsx(pls.22q.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights, "plsr.22q.DelvsHC.ZScore.CT.loadingweights.xlsx")







###########################################################################################################################################################################
#################################################### PLSR for CT with AHBA data re-processed with Abagen ##################################################################

#Prep Abagen Version of gene expression data
#Read in DK region labels used by abagen
abagen.DK.labels <- read_csv(paste0(data_folder, "atlas-desikankilliany_copy.csv"))

#subset to left cortical regions
abagen.DK.labels.LHcort <- abagen.DK.labels$label[abagen.DK.labels$structure == "cortex" & abagen.DK.labels$hemisphere == "L"]

#Read in abagen version of AHBA parcellated to surface-based DK atlas
abagen.geneexp.matrix <- read_csv(paste0(data_folder, "Expression_AHBA_abagen_maxprobe.csv")) %>% as.data.frame()

#Look at structure of expression object
dim(abagen.geneexp.matrix)
str(abagen.geneexp.matrix)

#Get list of genes in abagen expression object
abagen.genes <- colnames(abagen.geneexp.matrix)

#Subset expression data to only first 34 rows for left hemisphere regions, based on full 6 subjects (vs only 2 subjects with right hemisphere data in AHBA)
abagen.geneexp.matrix.LHcort <- abagen.geneexp.matrix[1:34,] %>% as.data.frame()

#Set rownames to DK atlas cortical regions
rownames(abagen.geneexp.matrix.LHcort) <- abagen.DK.labels.LHcort

#Load 22q11.2 deletion subjects versus healthy control map of z-score group effect sizes across left hemisphere DK atlas regions for cortical thickness
#From Supplementary Table 3 of Forsyth et al., 2021 (https://doi.org/10.1093/cercor/bhab008) 
ENIGMA22q.CT.deviance <- readRDS(paste0(data_folder, "ENIGMA22q_DelvsHC_ZScoreDeviance_CT.rda"))

#Look at structure of 22qDel deviance scores - how many brain structures do you see?
str(ENIGMA22q.CT.deviance)

#Make data frame for plsr with variable for 22qDel vs HC CT effect sizes across regions and matrix with AHBA gene expression across same regions as a variable
ENIGMA22q.CT.deviance.abagen.geneexp.zscore.reformat <- data.frame(ENIGMA.Dev.DelvsHC.ZScore.CT = ENIGMA22q.CT.deviance, AHBA.abagen.geneexp = I(as.matrix(abagen.geneexp.matrix.LHcort)))

#Run partial least square regression for 22qDel vs HC effect size map and AHBA gene expression data
pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.model = plsr(ENIGMA.Dev.DelvsHC.ZScore.CT ~ AHBA.abagen.geneexp, data = ENIGMA22q.CT.deviance.abagen.geneexp.zscore.reformat, validation="LOO")

#Check LOO CV-based RMSEP output for models with different numbers of principal components
summary(pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.model)

#RMSEP plot for CT Deviance
pdf("plsr.abagen.22qDelvsHC.ZScore.CT.numcomponents.byRMSEP.pdf")
validationplot(pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.model, val.type="RMSEP", legendpos = "bottomright") #create plot
dev.off() #Close the pdf file

#R2 plot for CT Deviance
pdf("plsr.abagen.22qDelvsHC.ZScore.CT.numcomponents.byR2.pdf")
plot(pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.model, "validation", val.type="R2") #create plot
dev.off() #Close the pdf file 

# Find the number of dimensions/components with lowest cross validation error
CV.CTDev = RMSEP(pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.model)
CTDev.best.dims = which.min(CV.CTDev$val[estimate = "adjCV", , ]) - 1

# Re-run plsr specifying 1 component
pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel = plsr(ENIGMA.Dev.DelvsHC.ZScore.CT ~ AHBA.abagen.geneexp, data = ENIGMA22q.CT.deviance.abagen.geneexp.zscore.reformat, ncomp = 1)

#Compute Observed RMSEP & R2 from Final model 
pls.pred.Del_HC.CT.zscore = predict(pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel, ENIGMA22q.CT.deviance.abagen.geneexp.zscore.reformat, ncomp=1)
pls.pred.Del_HC.CT.zscore.eval = data.frame(obs=ENIGMA22q.CT.deviance.abagen.geneexp.zscore.reformat$ENIGMA.Dev.DelvsHC.ZScore.CT, pred=pls.pred.Del_HC.CT.zscore[,1,1])
pls.pred.Del_HC.CT.zscore.eval.summary <- t(defaultSummary(pls.pred.Del_HC.CT.zscore.eval)) %>% as.data.frame()

#Extract gene loadings for all components output from plsr
pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel.Yscores <- Yscores(pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel)

#Extract gene loadings for 1st component output from plsr and format into percentiles
pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights <- data.frame(Loading.Weight = loading.weights(pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel)[,1])
pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights$Gene <- colnames(ENIGMA22q.CT.deviance.abagen.geneexp.zscore.reformat$AHBA.abagen.geneexp)
pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights$Loading.Weight.round <- round(pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights$Loading.Weight, 5)
pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights$Loading.Weight.perrank <- round(percent_rank(pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights$Loading.Weight.round),5)

####Note: Common downstream functional annotation of top PC1 loading genes includes gene ontology (GO) term enrichment, such as through g:profiler (https://biit.cs.ut.ee/gprofiler/gost) or similar tools
#Can also test enrichment of top PC1 load genes for cell-type specific expression, such as through pSI tool (http://doughertytools.wustl.edu/CSEAtool.html, https://sites.wustl.edu/doughertylab/psi_package-page/)

#Save output for all brain-expressed, protein-coding AHBA gene loading weights onto first PC for 22qDel deviance
write.xlsx(pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights, "plsr.abagen.AHBA.DelvsHC.ZScore.CT.loadingweights.xlsx")

#Extract percentile loadings onto 1st PC of 22q genes
pls.22q.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights <- pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights[pls.AHBA.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights$Gene %in% GeneList22q$HGNCUpdatedSymbol, ]

#Save output for brain-expressed, protein-coding 22q gene loading weights onto first PC for 22qDel deviance
write.xlsx(pls.22q.abagen.geneexp.DelvsHC.ZScore.CT.finalmodel.loadingweights, "plsr.abagen.22q.DelvsHC.ZScore.CT.loadingweights.xlsx")








