# ===============================================================================
# ACE.r
# Purpose: Perform ancestral state reconstruction for nematode life history traits
#   (parasitic vs free-living) using HiSSE models across multiple phylogenetic trees.
#   Includes model fitting, plotting, and rate comparison analysis.
# ===============================================================================

# Set the working directory for output files and data loading.
setwd('/To/Your/Directory/nematoda/lifehistory/res20260405')
# Load the saved R workspace containing previous analysis results.
load('./NemaLifeHist.RData')
# Source custom functions for HiSSE modeling and state rate testing.
source('./HiSSEModels.r')
source('./SSEStateRateTest.r')

# Load required R packages for data manipulation, phylogenetics, and visualization.
library(readxl)
library(ape)
library(vegan);library(ips);library(phytools)
library(ggtree);library(treeio)
library(future);library(doParallel)
library(hisse) # install_github(repo = "thej022214/hisse", ref = "master")
library(geiger)
library(ggplot2)

## Input data----
# Read the classification data for taxa.
Classification <- read.csv('/To/Your/Directory/nematoda/timing/Classification.csv', header = T, row.names = 1)

# Load phylogenetic trees from BEAST output files (184-taxon trees).
cLGPMSF <- treeio::read.beast('/To/Your/Directory/nematoda/timing/Res_posttrees/PRD-DBSCAN-AverageTree_Grp_2.tre')
C60PMSF <- treeio::read.beast('/To/Your/Directory/nematoda/timing/CalibrSensitivityAnalysis/AllPosteriorTimeTree/NMtimetree_C60PMSF.CST1.5P.IR.rg20_combined_posterior_tree_good.tre')
LG <- treeio::read.beast('/To/Your/Directory/nematoda/timing/CalibrSensitivityAnalysis/AllPosteriorTimeTree/NMtimetree_LG.CST1.5P.IR.rg20_combined_posterior_tree_good.tre')
Trad <- treeio::read.beast('/To/Your/Directory/nematoda/timing/CalibrSensitivityAnalysis/AllPosteriorTimeTree/NMtimetree_Trad.CST1.5P.IR.rg20_combined_posterior_tree_good.tre')

# Prune trees to focus on Nematoda and related groups, removing outgroups.
cLGPMSF <- drop.tip(cLGPMSF@phylo, tip = cLGPMSF@phylo$tip.label[1:25]) # Nematoda + Tardigrada
C60PMSF <- drop.tip(C60PMSF@phylo, tip = C60PMSF@phylo$tip.label[1:22]) # Nematoda + Tardigrada + Nematomorpha
LG <- drop.tip(LG@phylo, tip = LG@phylo$tip.label[1:22]) # Nematoda + Tardigrada + Nematomorpha
Trad <- drop.tip(Trad@phylo, tip = Trad@phylo$tip.label[1:25]) # Nematoda + Nematomorpha

# Load life history data (parasitic/free-living traits).
LifeHist<-read.csv('./NematodaTreeTaxonLifeStyle.csv',header = T)
rownames(LifeHist)<-LifeHist[,1]

## Data cleaning----
# Check which taxa in life history data match the tree tip labels.
# Some taxa may be excluded due to phylogenetic scope.
rownames(LifeHist)%in%cLGPMSF$tip.label # three FALSEs, Nematomorpha excluded
rownames(LifeHist)%in%C60PMSF$tip.label # all TRUE
rownames(LifeHist)%in%LG$tip.label # all TRUE
rownames(LifeHist)%in%Trad$tip.label # three FALSEs, Tardigrada excluded

# Define color schemes for plotting ancestral states.
ParasiticColorS<-c('#0000aa','#ffaadd','#4dda2a')
# HabitatColorS<-c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#f781bf', '#a65628')

## HiSSE----
# data of species richness: 
# Tardigrada: 1464 spp., according to (Degma and Guidetti 2023) online dataset, 2023-Jan-09
# Nematomorpha: 356 spp., according to Catalogue of Life online dataset, 2024-May-20
# Nematoda: 28537 described spp. according to (Hodda 2022, Phylum Nematoda: a classification). In a book published in 2013, Gardner (2013) estimated approximately 14,000â€?6,000 species of parasitic nematodes are parasitic (Encyclopedia of Biodiversity (Second Edition), when there were about 25043 species described (Zhang 2013, Zootaxa). From 2011 to 2011 to 2019 there are 3590 species were found as new, in which 1595 are free-living (44.42%) and 1995 are parasitic (55.57%) (Hodda 2022, Phylum Nematoda: trends in species descriptions). So, in a binary life style diviation, we estimate the number of parasitic species as (28537 - 25043) * 0.5557 + 16000 = 17942, that is the number of free-living species is 28537 - 17942 = 10595.
# Rhabdiasidae: 113 spp.; Strongyloididae: 82 spp.; Diplogasteromorpha: 483 spp. according to (Hodda 2022, Phylum Nematoda: a classification), totally 678 spp.

# For cLGPMSF, f=c(sum(LifeHist[cLGPMSF$tip.label, 'Parasitic']==0)/(10595+1464), sum(LifeHist[cLGPMSF$tip.label, 'Parasitic'])/(17942)) = c(0.005, 0.006)
# For C60PMSF, f=c(sum(LifeHist[C60PMSF$tip.label, 'Parasitic']==0)/(10595+1464), sum(LifeHist[C60PMSF$tip.label, 'Parasitic'])/(17942+356)) = c(0.005, 0.006)
# For LG, f=c(sum(LifeHist[LG$tip.label, 'Parasitic']==0)/(10595+1464), sum(LifeHist[LG$tip.label, 'Parasitic'])/(17942+356)) = c(0.005, 0.006)
# For Trad, f=c(sum(LifeHist[Trad$tip.label, 'Parasitic']==0)/(10595), sum(LifeHist[Trad$tip.label, 'Parasitic'])/(17942+356)) = c(0.005, 0.006)
# So, we use f=c(0.005, 0.006) for all the four HiSSE estimations.
ncrs <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
ncrs <- ifelse(ncrs > 81, 81, ncrs-1)

# Fit HiSSE models for parasitic lifestyle across different phylogenetic trees.
Parasitic_cLGPMSF_hisse <- HiSSEModels(phy=ladderize(cLGPMSF), data=LifeHist[cLGPMSF$tip.label, c('Tip_label', 'Parasitic')], f=c(0.005, 0.006), multiprocessors = ncrs)

Parasitic_C60PMSF_hisse <- HiSSEModels(phy=ladderize(C60PMSF), data=LifeHist[C60PMSF$tip.label, c('Tip_label', 'Parasitic')], f=c(0.005, 0.006), multiprocessors = ncrs)

Parasitic_LG_hisse <- HiSSEModels(phy=ladderize(LG), data=LifeHist[LG$tip.label, c('Tip_label', 'Parasitic')], f=c(0.005, 0.006), multiprocessors = ncrs)

Parasitic_Trad_hisse <- HiSSEModels(phy=ladderize(Trad), data=LifeHist[Trad$tip.label, c('Tip_label', 'Parasitic')], f=c(0.005, 0.006), multiprocessors = ncrs)

# Save workspace after parasitic models.
save.image('./NemaLifeHist.RData')

# For cLGPMSF, f=c(sum(LifeHist[cLGPMSF$tip.label, 'FreeLiving'])/(10595+1464+678), sum(LifeHist[cLGPMSF$tip.label, 'FreeLiving']==0)/(17942-678)) = c(0.005, 0.005)
# For C60PMSF, f=c(sum(LifeHist[C60PMSF$tip.label, 'FreeLiving'])/(10595+1464+678), sum(LifeHist[C60PMSF$tip.label, 'FreeLiving']==0)/(17942+356-678)) = c(0.005, 0.005)
# For LG, f=c(sum(LifeHist[LG$tip.label, 'FreeLiving'])/(10595+1464+678), sum(LifeHist[LG$tip.label, 'FreeLiving']==0)/(17942+356-678)) = c(0.005, 0.005)
# For Trad, f=c(sum(LifeHist[Trad$tip.label, 'FreeLiving'])/(10595+678), sum(LifeHist[Trad$tip.label, 'FreeLiving']==0)/(17942+356-678)) = c(0.006, 0.005)
# So, we use f=c(0.005, 0.005) for all the four HiSSE estimations but f=c(0.006, 0.005) for the Trad-tree-based reconstruction.

LHX <- LifeHist[,c('Tip_label', 'FreeLiving')]
LHX[, 'FreeLiving'] <- 1-LHX[, 'FreeLiving']

FreeLiving_cLGPMSF_hisse <- HiSSEModels(phy=ladderize(cLGPMSF), data=LHX[cLGPMSF$tip.label, ], f=c(0.005, 0.005), multiprocessors = ncrs)

FreeLiving_C60PMSF_hisse <- HiSSEModels(phy=ladderize(C60PMSF), data=LHX[C60PMSF$tip.label, ], f=c(0.005, 0.005), multiprocessors = ncrs)

FreeLiving_LG_hisse <- HiSSEModels(phy=ladderize(LG), data=LHX[LG$tip.label, ], f=c(0.005, 0.005), multiprocessors = ncrs)

FreeLiving_Trad_hisse <- HiSSEModels(phy=ladderize(Trad), data=LHX[Trad$tip.label, ], f=c(0.006, 0.005), multiprocessors = ncrs)

# Save workspace after free-living models.
save.image('./NemaLifeHist.RData')

## Plotting HiSSE models----
# Generate SVG plots for each lifestyle and phylogenetic model combination.
# Includes phylogram and fan tree views with state and rate information.
for (fp in c('Parasitic', 'FreeLiving')){
    for (PhyModel in c('cLGPMSF', 'C60PMSF', 'LG', 'Trad')){
        # Phylogram plot with tip labels.
        svg(paste('./', fp, '_', PhyModel, '_hisse_p.svg', sep = ''), onefile = T, width = 10, height = 20, pointsize = 10, bg = "transparent")
            assign(paste(fp, '_', PhyModel, '_hisse_p_treeplots', sep = ''), plot.hisse.states(x = get(paste(fp, '_', PhyModel, '_hisse', sep = ''))$Pred.RefinedModels[[1]], rate.param = "net.div", rate.colors=c('#559900', '#ff5500'), type = "phylogram", show.tip.label = TRUE, legend ='internal'))
        dev.off()
        
        # Side-by-side state and rate trees (phylogram).
        svg(paste('./', fp, '_', PhyModel, '_hisse_p_rs.svg', sep = ''), onefile = T, width = 20, height = 20, pointsize = 10, bg = "transparent")
            layout(matrix(1:2,1,2))
            get(paste(fp, '_', PhyModel, '_hisse_p_treeplots', sep = ''))$state.tree %>% plot(., direction='rightwards')
            get(paste(fp, '_', PhyModel, '_hisse_p_treeplots', sep = ''))$rate.tree %>% plot(., direction='leftwards')
        dev.off()
        
        # Fan plot without tip labels.
        svg(paste('./', fp, '_', PhyModel, '_hisse_c.svg', sep = ''), onefile = T, width = 10, height = 15, pointsize = 10, bg = "transparent")
            assign(paste(fp, '_', PhyModel, '_hisse_c_treeplots', sep = ''), plot.hisse.states(x = get(paste(fp, '_', PhyModel, '_hisse', sep = ''))$Pred.RefinedModels[[1]], rate.param = "net.div", rate.colors=c('#559900', '#ff5500'), type = "fan", show.tip.label = FALSE, legend = 'internal'))
        dev.off()
        
        # Side-by-side state and rate trees (fan view).
        svg(paste('./', fp, '_', PhyModel, '_hisse_c_rs.svg', sep = ''), onefile = T, width = 20, height = 15, pointsize = 10, bg = "transparent")
            layout(matrix(1:2,1,2))
            get(paste(fp, '_', PhyModel, '_hisse_c_treeplots', sep = ''))$state.tree %>% plot(., type='fan', show.tiplabels=FALSE)
            get(paste(fp, '_', PhyModel, '_hisse_c_treeplots', sep = ''))$rate.tree %>% plot(., type='fan', show.tip.label=FALSE)
        dev.off()
    }
}

## Phylogenetic comparison----

RateData_Comb <- NULL
RateData_KN_Comb <- NULL

for (AmbiLifeStyle in c('Parasitic', 'FreeLiving')){
    PhyModel_Comb <- NULL
    PhyModel_KN_Comb <- NULL
    for (PhyModel in c('cLGPMSF', 'C60PMSF', 'LG', 'Trad')){
        
        SSEStateRateTest(get(paste(AmbiLifeStyle,'_',PhyModel,'_hisse', sep = '')), bin.wd=10, cut.point=90, drop.outgroups=grep("^(NM_|TG_)", get(paste(AmbiLifeStyle,'_',PhyModel,'_hisse', sep = ''))$BestModel$phy$tip.label, value=TRUE), tag=paste(AmbiLifeStyle,'_',PhyModel, sep = ''), Classification=Classification)

        RateDataBest$Model <- 'Best'
        RateDataBest_KN$Model <- 'Best'
        RateDataModAvg$Model <- 'Averaged'
        RateDataModAvg_KN$Model <- 'Averaged'
        
        PhyModel_ALL <- rbind(RateDataModAvg, RateDataBest)
        PhyModel_KN <- rbind(RateDataModAvg_KN, RateDataBest_KN)
        
        PhyModel_ALL$AmbiLifeStyle <- AmbiLifeStyle
        PhyModel_KN$AmbiLifeStyle <- AmbiLifeStyle
        PhyModel_ALL$PhyModel <- PhyModel
        PhyModel_KN$PhyModel <- PhyModel
        
        PhyModel_Comb <- rbind(PhyModel_Comb, PhyModel_ALL)
        PhyModel_KN_Comb <- rbind(PhyModel_KN_Comb, PhyModel_KN)
        
        # in case that the function itself failed to print the figures
        # pdf(paste('./Rate Distributions Histogram for Nodes and Tips - ', AmbiLifeStyle, '_', PhyModel,  '.pdf', sep = ''), width = 60, height = 64, pointsize = 12)
        # ggsave(file = paste('./Rate Distributions Histogram for Nodes and Tips - ', AmbiLifeStyle, '_', PhyModel,  '.pdf', sep = ''), plot = mplt, width = 60, height = 64, limitsize = FALSE)        
        # dev.off()
    }
    RateData_Comb <- rbind(RateData_Comb, PhyModel_Comb)
    RateData_KN_Comb <- rbind(RateData_KN_Comb, PhyModel_KN_Comb)
}

RateData_Comb <- RateData_Comb[,c('PhyModel', 'AmbiLifeStyle', 'Model', 'nodeid', 'time', 'state', 'net.div', 	'turnover', 'speciation', 'extinction', 'state2', 'state4')]
RateData_KN_Comb <- RateData_KN_Comb[,c('Node', 'PhyModel', 'AmbiLifeStyle', 'Model', 'time', 'state', 'net.div', 	'turnover', 'speciation', 'extinction', 'state2', 'state4')]

RateData_KN_Comb <- RateData_KN_Comb %>% arrange(Node, PhyModel, AmbiLifeStyle, Model)

write.csv(RateData_Comb, file = './HiSSE_RateDataForBothNodesAndTips.csv', row.names = F)
write.csv(RateData_KN_Comb, file = './HiSSE_RateDataForKeyNodes.csv', row.names=F)

rm(AmbiLifeStyle, fp, PhyModel, RateDataBest, RateDataBest_KN, RateDataModAvg, RateDataModAvg_KN, PhyModel_Comb, PhyModel_KN_Comb, PhyModel_ALL, PhyModel_KN)

system('bash ./svg2pdf.sh r')

save.image('./NemaLifeHist.RData')
