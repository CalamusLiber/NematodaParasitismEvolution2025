# ==============================================================================
# SensAna.r - Sensitivity Analysis for MCMCtree Calibration Strategies
# ==============================================================================
#
# PURPOSE:
#   This R script performs comprehensive sensitivity analysis visualization for
#   different MCMCtree calibration strategies in molecular clock dating. It
#   compares divergence time estimates across various calibration approaches
#   and generates publication-quality figures showing the impact of different
#   fossil calibration strategies on phylogenetic dating results.
#
# WORKFLOW OVERVIEW:
#   1. Load multiple MCMCtree output files from different calibration strategies
#   2. Process tree files and extract divergence times for key clades
#   3. Generate comparative visualizations (individual trees and streamer plots)
#   4. Create geological time scale overlays and taxonomic color schemes
#   5. Export results as SVG files for publication
#
# INPUT FILES:
#   - Multiple .mcmc files from MCMCtree runs with different calibrations
#   - Classification.csv: Taxonomic classification data
#   - RDS files containing grouped tree data (GrpCntTrees_*.rds, GrpAvgTrees_*.rds)
#
# OUTPUT FILES:
#   - Individual tree plots: SA_TimeTree_*.svg
#   - Streamer plot: Key_Clade_Ages.svg
#   - Divergence time data: KeyCladeTime_SATrees.csv
#   - Combined analysis: GrpAvg_SA_Timetrees_simple.svg
#
# DEPENDENCIES:
#   - ape: Phylogenetic analysis
#   - ggtree: Tree visualization
#   - treeio: Tree data import/export
#   - grid: Grid graphics layout for composite plotting
#   - ggalt: Extended ggplot geoms for stream plots
#
# AUTHOR: [Molecular Clock Analysis Pipeline]
# DATE: [Current Version]
# ==============================================================================

## Summarize Sensitivity Analysis----
# Loading libraries----
library(ape)        # Core phylogenetic functions
library(ggtree)     # Advanced tree plotting
library(treeio)     # Tree data I/O and manipulation
library(grid)       # Grid graphics system
library(ggalt)      # Additional ggplot extensions
library(ggplot2)    # Plotting framework used for streamer plots
library(magrittr)   # Pipe operator (%>%) used in data processing

## Setting directory and loading image----
# setwd('/mnt/sda/nematoda/timing/CalibrSensitivityAnalysis/AllPosteriorTimeTree')
setwd('/To/Your/Directory/nematoda/timing/CalibrSensitivityAnalysis/AllPosteriorTimeTree')
load('/To/Your/Directory/nematoda/timing/SensAna.RData')
# save.image('/To/Your/Directory/nematoda/timing/SensAna.RData')

# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

# Extended color palette for taxonomic groups
# Provides additional colors when the main ClassColors palette is insufficient
ColorX <- c('#cc8cbc','#95a7a4','#948f6c','#609c70','#f05d4a','#aedbe1','#90447c','#002468','#194f1d','#5d1d07','#edc720','#dccaa8','#b59af2','#b56710','#dc6626','#42e476','#d0aca2','#6fc8d8','#a68a8d','#32a04d','#d7c1c0','#9d0000','#7c664e','#1d4d66','#8b4c39','#2982a2','#ba8547','#846671','#5f6b6b','#584f31','#4525b1','#cdcd00','#c4bda5','#ebd261','#ac9f77','#627f89','#9d5d8c','#aea88f','#146d1e','#71542b','#bfbfb4','#8b7a50','#1cbcd5','#bab590','#553e2d','#92310f','#00868b','#591618','#978466','#a4915c','#e29873','#36648b','#693e5a','#327570','#8c571b','#ba6452','#a787a1','#b5a49e','#ba3000','#646d4f','#d21111','#7c9d9e','#3c3520','#be9e1e','#c380d9','#579261','#90988e','#7f7f66','#ffdcca','#917a63','#b6cd2b','#bda683')

# EdgeColors function: Assigns colors to tree edges based on taxonomic classification
# This function creates color-coded branches for different taxonomic groups
# Parameters:
#   tree: phylogenetic tree object
#   cat: vector of taxonomic categories for each tip
# Returns: list containing edge colors and legend information
EdgeColors <- function(tree, cat){
    tips <- tree$tip.label
    Colors <- rep('black', nrow(tree$edge))
    CAT0 <- cbind(ColorX[1:length(unique(cat))], unique(cat))
    rownames(CAT0) <- unique(cat)
    cat <- cat[tips]
    CAT0 <- CAT0[unique(cat), ]
    CAT1 <- NULL
    for (i in 1:nrow(CAT0)){
        xxx <- which(cat==CAT0[i,2])
        icat <- c(1, which(diff(xxx) > 1) + 1)
        if (length(icat) > 1){
            newcat <- paste(CAT0[i,2], '(', as.character(1:length(icat)), ')', sep = '')
            CAT1 <- rbind(CAT1, cbind(CAT0[i,1], newcat))
            for (j in 1:length(icat)){
                cat[xxx[icat[j]:length(xxx)]] <- newcat[j]
            }
        }else{
            CAT1 <- rbind(CAT1, cbind(CAT0[i,1], CAT0[i,2]))
        }
    }
    SppYPos <- node.height(tree)[1:length(tips)]
    names(SppYPos) <- tips
    ncat1 <- nrow(CAT1)
    BarTxtPos <- data.frame(by0 = rep(NA, nrow = ncat1),
                            by1 = rep(NA, nrow = ncat1),
                            ty = rep(NA, nrow = ncat1),
                            btc = rep('black', nrow = ncat1))
    m <- mrca(tree)  # Most Recent Common Ancestor matrix
    for (i in 1:ncat1){
        # Find edges belonging to each taxonomic group
        ed <- unique(as.vector(m[tips[cat==CAT1[i,2]],tips[cat==CAT1[i,2]]]))
        Colors[tree$edge[,2] %in% ed] <- CAT1[i,1]
        BarTxtPos[i,4] <- CAT1[i,1]
        BarTxtPos[i,3] <- max(SppYPos[tips[cat==CAT1[i,2]]]) 
        BarTxtPos[i,2] <- BarTxtPos[i,3] + 0.1
        BarTxtPos[i,1] <- min(SppYPos[tips[cat==CAT1[i,2]]]) - 0.1
    }
    return(list(Colors, BarTxtPos, CAT0, CAT1[,2]))
}

# ==============================================================================
# DATA LOADING AND PREPARATION
# ==============================================================================

## Reading data----
# Load taxonomic classification data with full path
Classification <- read.csv('/To/Your/Directory/nematoda/timing/Classification.csv', header = T, row.names = 1)

# Create long and short classification names for different display purposes
ClassLong <- c(Classification[1:28,'Phylum'],Classification[-(1:28),'Order'])
ClassShort <- c(Classification[1:28,'Phylum'],Classification[-(1:28),'Clade'])
names(ClassLong) <- names(ClassShort) <- row.names(Classification)

# Identify sensitivity analysis tree files in the working directory
# Pattern matches MCMCtree output files with various calibration strategies
safiles <- dir(getwd())
satreefiles <- safiles[grep("^NMtimetree[A-Za-z0-9_\\.]*\\.tre", safiles)]

# Rearrange tree files for systematic plotting order
# Groups files by calibration strategy (CST1, CST2, CST3, CST4) for comparison
satreefiles <- satreefiles[c(grep('NMtimetree_cLGPMSF.CST[0-9]w.5P', satreefiles), 
                             grep('NMtimetree_cLGPMSF.CST[0-9]x.5P', satreefiles), 
                             grep('NMtimetree_cLGPMSF.CST[0-9]y.5P', satreefiles), 
                             grep('NMtimetree_C60PMSF', satreefiles), 
                             grep('NMtimetree_LG', satreefiles), 
                             grep('NMtimetree_Trad', satreefiles)
                             )]

# Initialize data structures for sensitivity analysis results
SATrees <- NULL  # Will store loaded tree objects
mcmc.sa <- data.frame(matrix(NA,nrow = length(satreefiles),ncol = 7))  # MCMC convergence data
colnames(mcmc.sa) <- c('cvg','mod','cst','fp','part','corr','rg')
mcmc.sa$cvg <- rep(NA,length(satreefiles))
BrTime <- matrix(NA,nrow = length(satreefiles), ncol = 184)  # Branching times matrix

# Load and process each sensitivity analysis tree file
for (i in 1:length(satreefiles)){
    # Read MCMCtree output file using treeio
    tre <- read.beast(satreefiles[i])
    
    # Parse filename to extract analysis parameters
    nm <- strsplit(substr(satreefiles[i], nchar('NMtimetree_')+1, nchar(satreefiles[i])-4), '_combined_posterior_tree_')[[1]]
    xxx <- strsplit(nm[1], '[.]')[[1]]
    
    # Record convergence status and analysis parameters
    mcmc.sa[i,] <- c(ifelse(nm[2]=='good','G','B'),  # G=good convergence, B=bad
                     # ifelse(nm[2]!='bad','G','B'), # Alternative convergence criterion
                     xxx[1],  # Model
                     substr(xxx[2],4,4),
                     toupper(gsub('CST[1-4]','',xxx[2])),
                     toupper(xxx[3]),
                     toupper(xxx[4]),
                     gsub('rg','',xxx[5]))
    BrTime[i,] <- branching.times(tre@phylo)
    SATrees[[i]] <- tre
}
mcmc.sa[mcmc.sa[,4]=='',4] <- '0'

MaxBrTimeCompan_SA <- unlist(lapply(SATrees, function(x){max(branching.times(x@phylo))}))
MaxBrTimeCompan_SA <- max(MaxBrTimeCompan_SA)-MaxBrTimeCompan_SA

# ==============================================================================
# VISUALIZATION SECTION
# ==============================================================================

# writing results----

# Change to results directory for output files
setwd('/To/Your/Directory/nematoda/timing/Res_posttrees')

# Create comprehensive SVG plot comparing all sensitivity analysis trees
# Each tree is plotted on its own panel to compare calibration strategies in series
svg(filename = './SA_TimeTrees.svg', width = 40, height = 240, pointsize = 12)
layout(mat = matrix(c(1:length(SATrees)),nrow = length(SATrees)))

# Plot each sensitivity analysis tree
for (i in 1:length(SATrees)){
    # Ladderize tree for consistent orientation and readability
    lad.tree <- ladderize(SATrees[[i]]@phylo, right = T)
    # Use taxonomic short labels for colored branches and legend
    EdgeBarTxtColors <- EdgeColors(lad.tree, ClassShort)
    
    # Compute the 95% HPD interval positions for each internal node
    nodebar_pos <- max(branching.times(lad.tree)) - Reduce(rbind, SATrees[[i]]@data$`0.95HPD`[order(as.numeric(SATrees[[i]]@data$node))])
    nodebar_pos <- cbind(nodebar_pos, node.height(lad.tree)[-c(1:(Ntip(lad.tree)+1))])
    
    # Plot the time tree without tip labels, using consistent axis limits
    plot(lad.tree, edge.color = EdgeBarTxtColors[[1]], show.tip.label = FALSE, x.lim = c(0,750)-MaxBrTimeCompan_SA[i], main = paste('Sensitivity Analysis Time Tree ', as.character(i), ': ', SATrees[[i]]@file, sep = ''), cex.main = 3, edge.width = 4)
    abline(v=max(branching.times(lad.tree))-seq(0,600,100), lty = 'solid', lend = 'butt', col = "gray30", lwd = 3)
    abline(v=max(branching.times(lad.tree))-seq(50,650,100), lty = 'dotdash', lend = 'butt', col = "gray60", lwd = 3)
    # Draw 95% HPD intervals for each internal node as horizontal blue segments
    segments(x0=nodebar_pos[,1], y0=nodebar_pos[,3], x1=nodebar_pos[,2], y1=nodebar_pos[,3], lwd = 4, col = 'royalblue', lend = 'butt')
    # Annotate calibrated nodes with a symbol if the node label contains '@'
    nodelabels(node = as.numeric(unlist(strsplit(lad.tree$node.label[grep('@',lad.tree$node.label)], '@'))) + ifelse(is.rooted(lad.tree), 1, 0), frame = 'none', pch = 8, adj = c(1.5,1), cex = 3)
    segments(x0 = max(branching.times(lad.tree)) + 7, y0 = EdgeBarTxtColors[[2]]$by0, x1 = max(branching.times(lad.tree)) + 7, y1 = EdgeBarTxtColors[[2]]$by1, lwd = 8, col = EdgeBarTxtColors[[2]]$btc, lend = 'butt')
    text(x = max(branching.times(lad.tree)) + 14, y = EdgeBarTxtColors[[2]]$ty, labels = EdgeBarTxtColors[[4]], cex = 3, col = EdgeBarTxtColors[[2]]$btc, adj = c(0,0.5))
    rect(max(branching.times(lad.tree))-c(66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8, 650), 185, max(branching.times(lad.tree))-c(0, 66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8), 191, border = NA, col = c('#8ecfc9','#ffbe7a','#fa7f6f','#82b0d2cc','#beb8dc','#e7dad2','#999900','#c497b2','#1450908a','#a9b8c6','#8e8bfe'))
    text(x=max(branching.times(lad.tree))-c(33, 105.5, 173.2, 226.651, 275.401, 328.9, 389.05, 431.5, 464.6, 512.1, 570), y=188, labels = c('Cz', 'K', 'J', 'T', 'P', 'Cn', 'D', 'S', 'O', 'Cm', 'PreCm'), cex = 2)
    if (i == length(SATrees)){
        axisPhylo(cex.axis=3, lwd = 3, padj = 1)
        mtext('Ma', side=1, at = 640, cex = 2, padj = 1.7)
        legend(x = 'bottomleft', legend = c('Calibrated nodes', EdgeBarTxtColors[[3]][,2]), pch = c(8, rep(15, nrow(EdgeBarTxtColors[[3]]))), col = c('black', EdgeBarTxtColors[[3]][,1]), text.col = c('black', EdgeBarTxtColors[[3]][,1]), cex = 2.5, bg = 'white', box.lwd =3)
    }
}
dev.off()

svg(filename = './SA_TimeTrees_simple.svg', width = 40, height = 8*length(SATrees), pointsize = 12)
layout(mat = matrix(c(1:length(SATrees)),nrow = length(SATrees)))
for (i in 1:length(SATrees)){
    lad.tree <- ladderize(SATrees[[i]]@phylo, right = T)
    EdgeBarTxtColors <- EdgeColors(lad.tree, ClassShort)
    plot(lad.tree, edge.color = EdgeBarTxtColors[[1]], show.tip.label = FALSE, x.lim = c(0,750)-MaxBrTimeCompan_SA[i], main = paste('Sensitivity Analysis Time Tree ', as.character(i), ': ', SATrees[[i]]@file, sep = ''), cex.main = 3, edge.width = 4)
    abline(v=max(branching.times(lad.tree))-seq(0,600,100), lty = 'solid', lend = 'butt', col = "gray30", lwd = 3)
    abline(v=max(branching.times(lad.tree))-seq(50,650,100), lty = 'dotdash', lend = 'butt', col = "gray60", lwd = 3)
    segments(x0 = max(branching.times(lad.tree)) + 7, y0 = EdgeBarTxtColors[[2]]$by0, x1 = max(branching.times(lad.tree)) + 7, y1 = EdgeBarTxtColors[[2]]$by1, lwd = 8, col = EdgeBarTxtColors[[2]]$btc, lend = 'butt')
    text(x = max(branching.times(lad.tree)) + 14, y = EdgeBarTxtColors[[2]]$ty, labels = EdgeBarTxtColors[[4]], cex = 3, col = EdgeBarTxtColors[[2]]$btc, adj = c(0,0.5))
    rect(max(branching.times(lad.tree))-c(66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8, 650), 185, max(branching.times(lad.tree))-c(0, 66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8), 191, border = NA, col = c('#8ecfc9','#ffbe7a','#fa7f6f','#82b0d2cc','#beb8dc','#e7dad2','#999900','#c497b2','#1450908a','#a9b8c6','#8e8bfe'))
    text(x=max(branching.times(lad.tree))-c(33, 105.5, 173.2, 226.651, 275.401, 328.9, 389.05, 431.5, 464.6, 512.1, 570), y=188, labels = c('Cz', 'K', 'J', 'T', 'P', 'Cn', 'D', 'S', 'O', 'Cm', 'PreCm'), cex = 2)
    if (i == length(SATrees)){
        axisPhylo(cex.axis=3, lwd = 3, padj = 1)
        mtext('Ma', side=1, at = 640, cex = 2, padj = 1.7)
        legend(x = 'bottomleft', legend = EdgeBarTxtColors[[3]][,2], pch = rep(15, nrow(EdgeBarTxtColors[[3]])), col = EdgeBarTxtColors[[3]][,1], text.col = EdgeBarTxtColors[[3]][,1], cex = 2.5, bg = 'white', box.lwd =3)
    }
}
dev.off()

# Create a 2-column layout SVG for the same sensitivity analysis trees
# This alternative layout facilitates side-by-side comparison of trees
svg(filename = './SA_TimeTrees_2C.svg', width = 80, height = 240, pointsize = 12)
layout(mat = matrix(c(1:length(SATrees)), ncol = 2, nrow = length(SATrees)/2, byrow = T))
for (i in 1:length(SATrees)){
    lad.tree <- ladderize(SATrees[[i]]@phylo, right = T)
    EdgeBarTxtColors <- EdgeColors(lad.tree, ClassShort)
    nodebar_pos <- max(branching.times(lad.tree)) - Reduce(rbind, SATrees[[i]]@data$`0.95HPD`[order(as.numeric(SATrees[[i]]@data$node))])
    nodebar_pos <- cbind(nodebar_pos, node.height(lad.tree)[-c(1:(Ntip(lad.tree)+1))])
    
    # Plot in a compact two-column arrangement without tip labels
    plot(lad.tree, edge.color = EdgeBarTxtColors[[1]], show.tip.label = FALSE, x.lim = c(0,750)-MaxBrTimeCompan_SA[i], main = paste('Sensitivity Analysis Time Tree ', as.character(i), ': ', SATrees[[i]]@file, sep = ''), cex.main = 3, edge.width = 4)
    abline(v=max(branching.times(lad.tree))-seq(0,600,100), lty = 'solid', lend = 'butt', col = "gray30", lwd = 3)
    abline(v=max(branching.times(lad.tree))-seq(50,650,100), lty = 'dotdash', lend = 'butt', col = "gray60", lwd = 3)
    segments(x0=nodebar_pos[,1], y0=nodebar_pos[,3], x1=nodebar_pos[,2], y1=nodebar_pos[,3], lwd = 4, col = 'royalblue', lend = 'butt')
    nodelabels(node = as.numeric(unlist(strsplit(lad.tree$node.label[grep('@',lad.tree$node.label)], '@'))) + ifelse(is.rooted(lad.tree), 1, 0), frame = 'none', pch = 8, adj = c(1.5,1), cex = 3)
    segments(x0 = max(branching.times(lad.tree)) + 7, y0 = EdgeBarTxtColors[[2]]$by0, x1 = max(branching.times(lad.tree)) + 7, y1 = EdgeBarTxtColors[[2]]$by1, lwd = 8, col = EdgeBarTxtColors[[2]]$btc, lend = 'butt')
    text(x = max(branching.times(lad.tree)) + 14, y = EdgeBarTxtColors[[2]]$ty, labels = EdgeBarTxtColors[[4]], cex = 3, col = EdgeBarTxtColors[[2]]$btc, adj = c(0,0.5))
    rect(max(branching.times(lad.tree))-c(66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8, 650), 185, max(branching.times(lad.tree))-c(0, 66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8), 191, border = NA, col = c('#8ecfc9','#ffbe7a','#fa7f6f','#82b0d2cc','#beb8dc','#e7dad2','#999900','#c497b2','#1450908a','#a9b8c6','#8e8bfe'))
    text(x=max(branching.times(lad.tree))-c(33, 105.5, 173.2, 226.651, 275.401, 328.9, 389.05, 431.5, 464.6, 512.1, 570), y=188, labels = c('Cz', 'K', 'J', 'T', 'P', 'Cn', 'D', 'S', 'O', 'Cm', 'PreCm'), cex = 2)
    if (i == length(SATrees)){
        # Add axis labels and legend only on the last panel to avoid redundancy
        axisPhylo(cex.axis=3, lwd = 3, padj = 1)
        mtext('Ma', side=1, at = 640, cex = 2, padj = 1.7)
        legend(x = 'bottomleft', legend = c('Calibrated nodes', EdgeBarTxtColors[[3]][,2]), pch = c(8, rep(15, nrow(EdgeBarTxtColors[[3]]))), col = c('black', EdgeBarTxtColors[[3]][,1]), text.col = c('black', EdgeBarTxtColors[[3]][,1]), cex = 2.5, bg = 'white', box.lwd =3)
    }
}
dev.off()

# Save each sensitivity analysis tree as an individual SVG file
# These individual plots include tip labels for detailed inspection
for (i in 1:length(SATrees)){
    svg(filename = paste('./SA_TimeTree-',as.character(i),'.svg', sep = ''), width = 50, height = 120, pointsize = 20)
    lad.tree <- ladderize(SATrees[[i]]@phylo, right = T)
    EdgeBarTxtColors <- EdgeColors(lad.tree, ClassLong)
    nodebar_pos <- max(branching.times(lad.tree)) - Reduce(rbind, SATrees[[i]]@data$`0.95HPD`[order(as.numeric(SATrees[[i]]@data$node))])
    nodebar_pos <- cbind(nodebar_pos, node.height(lad.tree)[-c(1:(Ntip(lad.tree)+1))])
    hpd95 <- unlist(SATrees[[i]]@data[order(as.numeric(SATrees[[i]]@data$node)),'0.95HPD']) %>% matrix(ncol = 2, byrow = T)
    plot(lad.tree, edge.color = EdgeBarTxtColors[[1]], show.tip.label = TRUE, x.lim = c(0,800)-MaxBrTimeCompan_SA[i], main = paste('Sensitivity Analysis Tree ', as.character(i), ': ', SATrees[[i]]@file, sep = ''), cex.main = 2, edge.width = 4)
    abline(v=max(branching.times(lad.tree))-seq(0,600,100), lty = 'solid', lend = 'butt', col = "gray30", lwd = 3)
    abline(v=max(branching.times(lad.tree))-seq(50,650,100), lty = 'dotdash', lend = 'butt', col = "gray60", lwd = 3)
    segments(x0=nodebar_pos[,1], y0=nodebar_pos[,3], x1=nodebar_pos[,2], y1=nodebar_pos[,3], lwd = 14, col = 'royalblue', lend = 'butt')
    nodelabels(node = as.numeric(unlist(strsplit(lad.tree$node.label[grep('@',lad.tree$node.label)], '@'))) + ifelse(is.rooted(lad.tree), 1, 0), frame = 'none', pch = 8, adj = c(1.5,1), cex = 2)
    segments(x0 = max(branching.times(lad.tree)) + 75, y0 = EdgeBarTxtColors[[2]]$by0, x1 = max(branching.times(lad.tree)) + 75, y1 = EdgeBarTxtColors[[2]]$by1, lwd = 16, col = EdgeBarTxtColors[[2]]$btc, lend = 'butt')
    text(x = max(branching.times(lad.tree)) + 80, y = EdgeBarTxtColors[[2]]$ty, labels = EdgeBarTxtColors[[4]], cex = 2, col = EdgeBarTxtColors[[2]]$btc, adj = c(0,0.5))
    rect(max(branching.times(lad.tree))-c(66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8, 650), 186, max(branching.times(lad.tree))-c(0, 66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8), 190, border = NA, col = c('#8ecfc9','#ffbe7a','#fa7f6f','#82b0d2cc','#beb8dc','#e7dad2','#999900','#c497b2','#1450908a','#a9b8c6','#8e8bfe'))
    text(x=max(branching.times(lad.tree))-c(33, 105.5, 173.2, 226.651, 275.401, 328.9, 389.05, 431.5, 464.6, 512.1, 570), y=188, labels = c('Cz', 'K', 'J', 'T', 'P', 'Cn', 'D', 'S', 'O', 'Cm', 'PreCm'), cex = 2)
    axisPhylo(cex.axis=3, lwd = 3, padj = 1)
    mtext('Ma', side=1, at = 640, cex = 3, padj = 1.7)
    legend(x = 'bottomleft', legend = c('Calibrated nodes', EdgeBarTxtColors[[3]][,2]), pch = c(8, rep(15, nrow(EdgeBarTxtColors[[3]]))), col = c('black', EdgeBarTxtColors[[3]][,1]), text.col = c('black', EdgeBarTxtColors[[3]][,1]), cex = 3, bg = 'white', box.lwd =3)
    dev.off()
}

rm(tre, nm, xxx, lad.tree, EdgeBarTxtColors, nodebar_pos, hpd95)

# ==============================================================================
# NODE AGE EXTRACTION FOR KEY CLADE COMPARISON
# ==============================================================================
# Define the clade hierarchy used for sensitivity comparison
# This list covers major metazoan lineages and focal nematode clades
keyclades <- c('root', 'EuMetazoa', 'Ecdysozoa', 'Arthropoda', 'Tardigrada', 'Nematomorpha', 'Nematoda', 'Enoplia (Clade II)', 'Dorylaimia (Clade I)', 'Chromadoria', 'Rhabditida','Spirurina (Clade III)', 'Tylenchina (Clade IV)', 'Rhabditina (Clade V)', 'Trichinellida', 'Tylenchomorpha', 'Diplogastromorpha', 'Strongyloididae', 'Ascarididae', 'Onchocercidae', 'Strongyloidea')

# GetTime: Extracts divergence times and 95% HPD intervals for key clades
# Arguments:
#   phy: a treeio phylo object containing tree and node age data
#   tre: 'cnt' for centroid tree HPDs, 'avg' for average tree HPDs
#   out: 'char' to return formatted strings, 'num' to return numeric matrix
GetTime <- function(phy, tre='cnt', out = 'char'){
    keynodes <- c(getMRCA(phy@phylo, tip = row.names(Classification)), 
                  getMRCA(phy@phylo, tip = row.names(Classification)[-1]),
                  getMRCA(phy@phylo, tip = row.names(Classification)[-1:-9]),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Clade=='Arthropoda']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Clade=='Tardigrada']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Clade=='Nematomorpha']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Phylum %in% c('Chromadoria','Dorylaimia','Enoplia')]),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Clade=='Enoplia (Clade II)']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Clade=='Dorylaimia (Clade I)']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Phylum=='Chromadoria']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Clade %in% c('Spirurina (Clade III)','Tylenchina (Clade IV)','Rhabditina (Clade V)')]),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Clade=='Spirurina (Clade III)']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Clade=='Tylenchina (Clade IV)']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Clade=='Rhabditina (Clade V)']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Order=='Trichinellida']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Order=='Tylenchomorpha']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Order=='Diplogastromorpha']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Family=='Strongyloididae']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Family=='Ascarididae']),
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Family=='Onchocercidae']), 
                  getMRCA(phy@phylo, tip = row.names(Classification)[Classification$Order=='Strongyloidea'])
    )
    if (out == 'char'){
        # Return formatted text strings with node age and HPD interval for each target clade
        if (tre=='cnt'){
            hdi <- unlist(Map(function(x){paste(round(unlist(phy@data[phy@data$node==x, '0.95HPD']), digit=2), collapse=',')}, keynodes))
        }else if (tre=='avg'){
            hdi <- unlist(Map(function(x){paste(round(unlist(phy@data[phy@data$node==x, '0.95HPD_RLX']), digit=2), collapse=',')}, keynodes))
        }
        bt <- round(branching.times(phy@phylo)[keynodes - Ntip(phy@phylo)],2)
        res <- paste(bt, '[', hdi, ']', sep='')
        names(res) <- keyclades
    }else if (out == 'num'){
        # Return numeric matrix with branch age and HPD interval bounds
        if (tre=='cnt'){
            res <- cbind(bt = branching.times(phy@phylo)[keynodes - Ntip(phy@phylo)], t(sapply(Map(function(x){unlist(phy@data[phy@data$node==x, '0.95HPD'])}, keynodes), c)))
        }else if (tre=='avg'){
            res <- cbind(bt = branching.times(phy@phylo)[keynodes - Ntip(phy@phylo)], t(sapply(Map(function(x){unlist(phy@data[phy@data$node==x, '0.95HPD_RLX'])}, keynodes), c)))
        }
        rownames(res) <- keyclades
    }
    return(res)
}

# Export key divergence times for each sensitivity analysis tree
# The output file includes clade age and HPD intervals for major clades
write.csv(cbind(Cluster=c(1:length(SATrees)), Reduce(rbind, lapply(SATrees, GetTime, tre='cnt')), Treefile=unlist(Map(function(x){SATrees[[x]]@file}, 1:length(SATrees)))), file = './KeyCladeTime_SATrees.csv', row.names = F)

# saveRDS(SATrees, file = './SATrees.rds')

## streamer plot----
# Load precomputed grouped tree results for major cluster summaries
GrpCntTrees_tSNE <- readRDS('./GrpCntTrees_tSNE.rds')
GrpAvgTrees_tSNE <- readRDS('./GrpAvgTrees_tSNE.rds')
GrpCntTrees_PRD <- readRDS('./GrpCntTrees_PRD.rds')
GrpAvgTrees_PRD <- readRDS('./GrpAvgTrees_PRD.rds')
GrpCntTrees_AEC <- readRDS('./GrpCntTrees_AEC.rds')
GrpAvgTrees_AEC <- readRDS('./GrpAvgTrees_AEC.rds')

# Convert grouped tree objects to key clade age matrices for plotting
GrpCntTrees_tSNE_CladeAges <- lapply(GrpCntTrees_tSNE, GetTime, tre='cnt', out='num')
GrpAvgTrees_tSNE_CladeAges <- lapply(GrpAvgTrees_tSNE, GetTime, tre='avg', out='num')
GrpCntTrees_PRD_CladeAges <- lapply(GrpCntTrees_PRD, GetTime, tre='cnt', out='num')
GrpAvgTrees_PRD_CladeAges <- lapply(GrpAvgTrees_PRD, GetTime, tre='avg', out='num')
GrpCntTrees_AEC_CladeAges <- lapply(GrpCntTrees_AEC, GetTime, tre='cnt', out='num')
GrpAvgTrees_AEC_CladeAges <- lapply(GrpAvgTrees_AEC, GetTime, tre='avg', out='num')

# Extract corresponding clade ages for the original sensitivity analysis trees
SATrees_CladeAges <- lapply(SATrees, GetTime, tre='cnt', out='num')

# Build the combined list of clade age matrices for streamer plotting
StrmPlotCladeAges <- NULL
for (i in 1:length(GrpCntTrees_tSNE_CladeAges)) {
    StrmPlotCladeAges <- append(StrmPlotCladeAges, GrpCntTrees_tSNE_CladeAges[i])
    StrmPlotCladeAges <- append(StrmPlotCladeAges, GrpAvgTrees_tSNE_CladeAges[i])
}
for (i in 1:length(GrpCntTrees_PRD_CladeAges)) {
    StrmPlotCladeAges <- append(StrmPlotCladeAges, GrpCntTrees_PRD_CladeAges[i])
    StrmPlotCladeAges <- append(StrmPlotCladeAges, GrpAvgTrees_PRD_CladeAges[i])
}

for (i in 1:length(GrpCntTrees_AEC_CladeAges)) {
    StrmPlotCladeAges <- append(StrmPlotCladeAges, GrpCntTrees_AEC_CladeAges[i])
    StrmPlotCladeAges <- append(StrmPlotCladeAges, GrpAvgTrees_AEC_CladeAges[i])
}

# Append the original sensitivity analysis trees to the streamer dataset
StrmPlotCladeAges <- append(StrmPlotCladeAges, SATrees_CladeAges)

# Build convergence status vector used to color points by convergence quality
# 'G' indicates good convergence; bad runs use red points
StrmPlotCladeCvgS <- c(rep('G', length(GrpCntTrees_tSNE_CladeAges) + length(GrpCntTrees_PRD_CladeAges) + length(GrpCntTrees_AEC_CladeAges), each = 2), mcmc.sa$cvg)

# Generate the streamer plot showing clade ages across cluster summaries and sensitivity analysis trees
svg('./Key_Clade_Ages.svg', width = 18, height = 18, pointsize = 12)
pushViewport(viewport(layout = grid.layout(nrow = nrow(StrmPlotCladeAges[[1]]), ncol = 1)))
for (i in 1:nrow(StrmPlotCladeAges[[1]])){
    df <- as.data.frame(Reduce(rbind, lapply(StrmPlotCladeAges, function(x, r){x[r,]}, r = i)))
    if (i == nrow(StrmPlotCladeAges[[1]])){
        cldsp0 <- ggplot()+
            stat_xspline(aes(x=1:length(StrmPlotCladeAges), y= df[,2]), spline_shape = -0.5, col = 'gray80')+
            stat_xspline(aes(x=1:length(StrmPlotCladeAges), y= df[,3]), spline_shape = -0.5, col = 'gray80')+
            theme_bw()+
            theme(panel.grid = element_blank(),
                  axis.text = element_text(size=6),
                  axis.title.y = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=6),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = 'black'))+
            ylab(keyclades[i])
    } else {
        cldsp0 <- ggplot()+
            stat_xspline(aes(x=1:length(StrmPlotCladeAges), y= df[,2]), spline_shape = -0.5, col = 'gray80')+
            stat_xspline(aes(x=1:length(StrmPlotCladeAges), y= df[,3]), spline_shape = -0.5, col = 'gray80')+
            theme_bw()+
            theme(panel.grid = element_blank(),
                  axis.text = element_text(size=6),
                  axis.title.y = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=6),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.line.y = element_line(color = 'black'),
                  axis.line.x = element_blank())+
            ylab(keyclades[i])
    }
    gg <- ggplot_build(cldsp0)
    cldsp1 <- cldsp0+
        geom_polygon(mapping=aes(x=c(gg$data[[1]]$x, rev(gg$data[[2]]$x)),
                                 y=c(gg$data[[1]]$y, rev(gg$data[[2]]$y))),
                     fill='gray80', alpha=1)+
        geom_xspline(mapping=aes(x=1:length(StrmPlotCladeAges), y= df[,1]),
                     spline_shape = -0.5, col = 'gray40', linetype = 1)+
        geom_point(mapping=aes(x=1:length(StrmPlotCladeAges), y= df[,1]),
                   size = 0.8, col = ifelse(StrmPlotCladeCvgS=='G','gray20','#e43d32ff'))
    if (i == nrow(StrmPlotCladeAges[[1]])){
        # Showing the X-label under the last row
        cldsp2 <- cldsp1+
            scale_x_discrete(
                name = NULL,   # no axis name
                limits = factor(1:length(StrmPlotCladeAges)),
                labels = c(rep(c('C','A'), 0.5 * (length(StrmPlotCladeAges)-length(SATrees_CladeAges))),
                           rep(c('CST1','CST2','CST3','CST4'), length(SATrees_CladeAges)/4)))+ # tick labels
            xlab(NULL) 
    } else {
        # Masking the X-label under the rows before last
        cldsp2 <- cldsp1+
            scale_x_discrete(
                name = NULL, 
                limits = factor(1:length(StrmPlotCladeAges)),
                labels = NULL) +  # tick labels (hidden)
            xlab(NULL)+ 
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank())
    }
    print(cldsp2, vp = viewport(layout.pos.row = i, layout.pos.col = 1))
}
dev.off()

# Plot average time trees for cluster summaries plus representative sensitivity analysis trees
# This small-multiples figure allows visual comparison of average cluster trees and select SA trees
GrpAvg_SA_Timetrees <- c(GrpAvgTrees_PRD, SATrees[grep('CST1', satreefiles)]) # Use CST1 representative trees for the sensitivity analysis block

# Align x-axis limits across all trees for consistent time scaling
MaxBrTimeAVG_SA <- unlist(lapply(GrpAvg_SA_Timetrees, function(x){max(branching.times(x@phylo))}))
MaxBrTimeAVG_SA <- max(MaxBrTimeAVG_SA)-MaxBrTimeAVG_SA

svg(filename = './GrpAvg_SA_Timetrees_simple.svg', width = 40, height = 100, pointsize = 12)
layout(mat = matrix(c(1:length(GrpAvg_SA_Timetrees)),nrow = length(GrpAvg_SA_Timetrees)))
j <- 1
for (i in 1:length(GrpAvg_SA_Timetrees)){
    lad.tree <- ladderize(GrpAvg_SA_Timetrees[[i]]@phylo, right = T)
    EdgeBarTxtColors <- EdgeColors(lad.tree, ClassShort)
    plot(lad.tree, edge.color = EdgeBarTxtColors[[1]], show.tip.label = FALSE, x.lim = c(0,750)-MaxBrTimeAVG_SA[i], main = '', cex.main = 3, edge.width = 4)
    abline(v=max(branching.times(lad.tree))-seq(0,600,100), lty = 'solid', lend = 'butt', col = "gray30", lwd = 3)
    abline(v=max(branching.times(lad.tree))-seq(50,650,100), lty = 'dotdash', lend = 'butt', col = "gray60", lwd = 3)
    segments(x0 = max(branching.times(lad.tree)) + 7, y0 = EdgeBarTxtColors[[2]]$by0, x1 = max(branching.times(lad.tree)) + 7, y1 = EdgeBarTxtColors[[2]]$by1, lwd = 8, col = EdgeBarTxtColors[[2]]$btc, lend = 'butt')
    if (i == 1){
        # Add legend labels and geological period bar only on first panel for clarity
        text(x = max(branching.times(lad.tree)) + 14, y = EdgeBarTxtColors[[2]]$ty, labels = EdgeBarTxtColors[[4]], cex = 3, col = EdgeBarTxtColors[[2]]$btc, adj = c(0,0.5))
        rect(max(branching.times(lad.tree))-c(66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8, 650), 185, max(branching.times(lad.tree))-c(0, 66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8), 191, border = NA, col = c('#8ecfc9','#ffbe7a','#fa7f6f','#82b0d2cc','#beb8dc','#e7dad2','#999900','#c497b2','#1450908a','#a9b8c6','#8e8bfe'))
        text(x=max(branching.times(lad.tree))-c(33, 105.5, 173.2, 226.651, 275.401, 328.9, 389.05, 431.5, 464.6, 512.1, 570), y=188, labels = c('Cz', 'K', 'J', 'T', 'P', 'Cn', 'D', 'S', 'O', 'Cm', 'PreCm'), cex = 2)
    }
    if (i == length(GrpAvg_SA_Timetrees)){
        # Add axis labels, geological time bar, and legend on the final panel
        rect(max(branching.times(lad.tree))-c(66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8, 650), -6, max(branching.times(lad.tree))-c(0, 66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8), 0, border = NA, col = c('#8ecfc9','#ffbe7a','#fa7f6f','#82b0d2cc','#beb8dc','#e7dad2','#999900','#c497b2','#1450908a','#a9b8c6','#8e8bfe'))
        text(x=max(branching.times(lad.tree))-c(33, 105.5, 173.2, 226.651, 275.401, 328.9, 389.05, 431.5, 464.6, 512.1, 570), y=-3, labels = c('Cz', 'K', 'J', 'T', 'P', 'Cn', 'D', 'S', 'O', 'Cm', 'PreCm'), cex = 2)
        axisPhylo(cex.axis=3, lwd = 3, padj = 6)
        mtext('Ma', side=1, at = 640, cex = 2, padj = 1.7)
        legend(x = 'bottomleft', legend = EdgeBarTxtColors[[3]][,2], pch = rep(15, nrow(EdgeBarTxtColors[[3]])), col = EdgeBarTxtColors[[3]][,1], text.col = EdgeBarTxtColors[[3]][,1], cex = 2.5, bg = 'white', box.lwd =3)
    }
}
dev.off()

# Clean up temporary plotting variables
rm(i, j, gg, df, cldsp0, cldsp1, cldsp2, lad.tree)

# Convert generated SVGs to PDF using an external script
system('bash ./svg2pdf.sh r')

# Save the current R workspace for future continuation of analysis
save.image('/To/Your/Directory/nematoda/timing/SensAna.RData')
