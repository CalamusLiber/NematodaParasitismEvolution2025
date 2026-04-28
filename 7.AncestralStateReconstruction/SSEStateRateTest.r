# SSEStateRateTest function: Tests and visualizes state-specific evolutionary rates from HiSSE models
#
# This function processes results from HiSSE (Hidden State Speciation and Extinction) models to:
# - Bin phylogenetic nodes by time periods
# - Categorize ancestral states (free-living vs parasitic lifestyles)
# - Perform statistical tests comparing rates between states and time periods
# - Generate plots of rate changes through time and rate distributions
#
# Parameters:
#   SSE: HiSSE model results object containing fitted models and data
#   bin.wd: Width of time bins in million years (default: 10)
#   cut.point: Optional numeric time point (Ma) to split analysis into before/after periods
#   tag: Identifier string for output files (e.g., "analysis_run1")
#   drop.outgroups: Optional vector of outgroup tip names to exclude from analysis
#   Classification: Data frame with taxonomic classification for tips
#
# Outputs:
#   - CSV files with rate data and statistical test results
#   - Text files with ANCOVA and ANOVA reports
#   - SVG plots of rate changes through time and rate distributions
#
# Author: Liang Lü
# Date: Implementation for nematode evolutionary analysis

SSEStateRateTest <- function(SSE, bin.wd=10, cut.point=NULL, tag=NULL, drop.outgroups=NULL, Classification=NULL){
    # Load required R packages for data manipulation, plotting, and statistical analysis
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(agricolae)
    
    options(warn = -1)
    
    # Extract phylogenetic tree from HiSSE results
    Phy <- SSE$BestModel$phy
    
    # Define geological time scale boundaries (in million years ago)
    GeoScale <- c(0, 66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8)
    # Colors for geological period visualization
    GeoColor <- c('gray90','#8ecfc9','#ffbe7a','#fa7f6f','#82b0d2cc','#beb8dc','#e7dad2','#999900','#c497b2','#1450908a','#a9b8c6','#8e8bfe')
    # Labels for geological periods
    GeoLabel <- c('Pst', 'Cz', 'K', 'J', 'T', 'P', 'Cn', 'D', 'S', 'O', 'Cm', 'PreCm')
    
    # Define key nematode clades for focused analysis
    keyclades <- c('Nematoda', 'Enoplia (Clade II)', 'Dorylaimia (Clade I)', 'Chromadoria', 'Rhabditida','Spirurina (Clade III)', 'Tylenchina (Clade IV)', 'Rhabditina (Clade V)', 'Trichinellida', 'Tylenchomorpha', 'Diplogastromorpha', 'Strongyloididae', 'Ascarididae', 'Onchocercidae', 'Strongyloidea')
    
    # Get node IDs for key clades using most recent common ancestors (MRCA)
    keynodes <- c(getMRCA(Phy, tip = row.names(Classification)[Classification$Phylum %in% c('Chromadoria','Dorylaimia','Enoplia')]),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Clade=='Enoplia (Clade II)']),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Clade=='Dorylaimia (Clade I)']),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Phylum=='Chromadoria']),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Clade %in% c('Spirurina (Clade III)','Tylenchina (Clade IV)','Rhabditina (Clade V)')]),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Clade=='Spirurina (Clade III)']),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Clade=='Tylenchina (Clade IV)']),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Clade=='Rhabditina (Clade V)']),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Order=='Trichinellida']),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Order=='Tylenchomorpha']),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Order=='Diplogastromorpha']),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Family=='Strongyloididae']),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Family=='Ascarididae']),
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Family=='Onchocercidae']), 
                  getMRCA(Phy, tip = row.names(Classification)[Classification$Order=='Strongyloidea'])
    ) %>% as.character
    
    # Separate nodes into time bins and handle outgroup removal ----
    # Determine which tips and nodes to keep (exclude outgroups if specified)
    if (all(c(!is.null(drop.outgroups),drop.outgroups%in%Phy$tip.label))){
        KeepTipIDs <- which(!Phy$tip.label%in%drop.outgroups)
        KeepTips <- Phy$tip.label[KeepTipIDs]
        KeepNodeIDs <- getMRCA(Phy, KeepTips):max(Phy$edge[Phy$edge[,2]%in%KeepTipIDs,1]) %>% as.character
        KeepTipIDs <- as.character(KeepTipIDs)
        KeepParentIDs <- as.character(Phy$edge[,1])
        names(KeepParentIDs) <- as.character(Phy$edge[,2])
        KeepParentIDs <- KeepParentIDs[c(KeepNodeIDs, KeepTipIDs)]
        Phy <- keep.tip(Phy, KeepTips)
    }else{
        KeepTipIDs <- as.character(1:Ntip(Phy))
        KeepTips <- Phy$tip.label
        KeepNodeIDs <- as.character(Ntip(Phy)+1:Phy$Nnode)
        KeepParentIDs <- as.character(c(Phy$edge[1,1],Phy$edge[,1]))
        names(KeepParentIDs) <- as.character(c(Phy$edge[1,1],Phy$edge[,2]))
        KeepParentIDs <- KeepParentIDs[c(KeepNodeIDs, KeepTipIDs)]
    }
    
    # Calculate node ages and create time bins ----
    # Set tip ages to present (0 Ma, but use -bin.wd for binning)
    tiptime <- rep(-bin.wd, length(KeepTips))
    names(tiptime) <- KeepTipIDs
    # Get branching times for internal nodes
    time <- c(tiptime, branching.times(SSE$BestModel$phy)[KeepNodeIDs])
    
    # Initialize bin containers
    bins <- NULL
    bin.idx <- rep(NA, length(time))
    names(bin.idx) <- names(time)
    
    # Assign nodes to time bins based on their ages
    for (i in 1:ceiling(diff(range(time))/bin.wd)){
        bins[[i]] <- names(time[(time >= bin.wd * (i - 2)) & (time < bin.wd * (i - 1))])
        bin.idx[(time >= bin.wd * (i - 2)) & (time < bin.wd * (i - 1))] <- i
    }
    
    # Create time-rate-state datasets for nodes ----
    if (is.numeric(cut.point)){
        FUN01 <- function(x){x <- x %>% mutate(state2=case_when((state < 0.5) ~ 'FreeLiving',
                                                                (state >= 0.5) ~ 'Parasitic'), 
                                               state4=case_when((state <= 0.05) ~ 'HardFreeLiving',
                                                                (state > 0.05 & state < 0.5) ~ 'SoftFreeLiving',
                                                                (state >= 0.5 & state < 0.95) ~ 'SoftParasitic',
                                                                (state >= 0.95) ~ 'HardParasitic'),
                                               BAcut=case_when((time >= cut.point) ~ 'Before',
                                                               (time < cut.point) ~ 'After', 
                                                               (time < 0) ~ 'Present'))}
    }else{
        FUN01 <- function(x){x <- x %>% mutate(state2=case_when((state < 0.5) ~ 'FreeLiving',
                                                                (state >= 0.5) ~ 'Parasitic'), 
                                               state4=case_when((state <= 0.05) ~ 'HardFreeLiving',
                                                                (state > 0.05 & state < 0.5) ~ 'SoftFreeLiving',
                                                                (state >= 0.5 & state < 0.95) ~ 'SoftParasitic',
                                                                (state >= 0.95) ~ 'HardParasitic'))}
    }
    
    # Define rate variables to analyze
    MeasureVars <- c('net.div','turnover','speciation','extinction')
    
    # Prepare scaled parent rates (baseline rates for ANCOVA)
    prt <- rbind(apply(SSE$ModelAveRates.BestModel$nodes[KeepParentIDs[KeepTipIDs], MeasureVars], 2, scale), 
                 apply(SSE$ModelAveRates.BestModel$nodes[KeepParentIDs[KeepNodeIDs], MeasureVars], 2, scale))
    colnames(prt) <- paste(MeasureVars, '0', sep = '')
    # Get descendant rates and states
    dgt <- rbind(SSE$ModelAveRates.BestModel$tips[KeepTipIDs,c('state', MeasureVars)], 
                 SSE$ModelAveRates.BestModel$nodes[KeepNodeIDs,c('state', MeasureVars)])
    # Create best model rate dataset
    RateDataBest <<- cbind(nodeid=names(time), time, bin=bin.idx, prt, dgt) %>% FUN01
    RateDataBest_KN <<- RateDataBest[keynodes,]
    RateDataBest_KN$Node <<- keyclades
    
    # Same for model-averaged rates
    prt <- rbind(apply(SSE$ModelAveRates$nodes[KeepParentIDs[KeepTipIDs], MeasureVars], 2, scale),
                 apply(SSE$ModelAveRates$nodes[KeepParentIDs[KeepNodeIDs], MeasureVars], 2, scale))
    colnames(prt) <- paste(MeasureVars, '0', sep = '')
    dgt <- rbind(SSE$ModelAveRates$tips[KeepTipIDs,c('state', MeasureVars)], 
                 SSE$ModelAveRates$nodes[KeepNodeIDs,c('state', MeasureVars)])
    RateDataModAvg <<- cbind(nodeid=names(time), time, bin=bin.idx, prt, dgt) %>% FUN01
    RateDataModAvg_KN <<- RateDataModAvg[keynodes,]
    RateDataModAvg_KN$Node <<- keyclades
    
    # Function to extract rates from individual models
    FUN02 <- function(x){
        RT <- GetModelAveRates(x, AIC.weights=NULL, type="both")
        res <- FUN01(cbind(time, bin=bin.idx, rbind(RT$tips[KeepTipIDs,c('state', MeasureVars)],RT$nodes[KeepNodeIDs,c('state', MeasureVars)])))
        return(res)
    }
    
    # Generate rate datasets for all refined models (commented out for models with weight > 1e-6)
    # RateDataModAll <- lapply(SSE$Pred.RefinedModels[SSE$ModelSelection[,'w'] > 1e-6], FUN02)
    RateDataModAll <- lapply(SSE$Pred.RefinedModels, FUN02)
    # RateDataModAll <- lapply(RateDataModAll, FUN01)
    # RateDataModAll <- array(unlist(RateDataModAll), dim=c(dim(RateDataModAll[[1]]), length(RateDataModAll)),dimnames = list(rownames(RateDataModAll[[1]]), colnames(RateDataModAll[[1]]), as.character(1:length(RateDataModAll))))
    
    # Export model selection and rate data to CSV files ----
    write.csv(SSE$ModelSelection, file = paste('./HiSSE_ModelSelection - ', tag, '.csv', sep = ''), row.names = TRUE)
    
    write.csv(RateDataBest, file = paste('./HiSSE_RateDataBestForBothNodesAndTips - ', tag, '.csv', sep = ''), row.names = FALSE)
    write.csv(RateDataBest_KN, file = paste('./HiSSE_RateDataBestForKeyNodes - ', tag, '.csv', sep = ''), row.names = FALSE)
    
    write.csv(RateDataModAvg, file = paste('./HiSSE_RateDataModAvgForBothNodesAndTips - ', tag, '.csv', sep = ''), row.names = FALSE)
    write.csv(RateDataModAvg_KN, file = paste('./HiSSE_RateDataModAvgForKeyNodes - ', tag, '.csv', sep = ''), row.names = FALSE)
    
    ## Plot rate change through time ----
    svg(filename = paste('./Rate Change Through Time - ', tag, '.svg', sep = ''), width = 60, height = 56, pointsize = 12)
    layout(mat = matrix(1:28, 7, 4, byrow = T))
    
    MeanBest <- StdBest <- MeanAvg <- StdAvg <- MaxAll0 <- MinAll0 <- data.frame(matrix(NA, length(bins), 4), row.names = as.character(1:length(bins)))
    colnames(MeanBest) <- colnames(StdBest) <- colnames(MeanAvg) <- colnames(StdAvg) <- colnames(MaxAll0) <- colnames(MinAll0) <- MeasureVars
    
    # rate-for-nodes dataset ----
    for (i in 1:length(bins)){
        if (any(bin.idx==i)){
            MeanBest[i,] <- apply(RateDataBest[bin.idx==i,MeasureVars], 2, mean, na.rm=T)
            StdBest[i,] <- apply(RateDataBest[bin.idx==i,MeasureVars], 2, sd, na.rm=T)
            MeanAvg[i,] <- apply(RateDataModAvg[bin.idx==i,MeasureVars], 2, mean, na.rm=T)
            StdAvg[i,] <- apply(RateDataModAvg[bin.idx==i,MeasureVars], 2, sd, na.rm=T)
            MaxAll0[i,] <- bind_rows(lapply(RateDataModAll, FUN=function(x){apply(x[bin.idx==i,MeasureVars], 2, max, na.rm=T)})) %>% apply(2, max, na.rm=T)
            MinAll0[i,] <- bind_rows(lapply(RateDataModAll, FUN=function(x){apply(x[bin.idx==i,MeasureVars], 2, min, na.rm=T)})) %>% apply(2, min, na.rm=T)
        }
    }
    
    for (var in MeasureVars){
        plot(x=-bin.wd*(1:length(bins)-1)+0.5*bin.wd, y=MeanAvg[,var], 
             type='n', main=ifelse(is.null(tag),'Rate Change Through Time (All)',paste('Rate Change Through Time (', tag, '-All)', sep='')), 
             xlab='Million Years Before Present', ylab=var, 
             xlim=c(-bin.wd*(length(bins)-1), bin.wd), ylim=c(min(MinAll0[,var], na.rm = T)*0.9, max(MaxAll0[,var], na.rm = T)*1.1), cex=4)
        if (! is.null(cut.point)){
            abline(v=-cut.point, col = 'gray30', lty = 2, lwd = 2)
        }
        rect(-bin.wd*(1:length(bins)-1), MinAll0[,var], -bin.wd*(1:length(bins)-2), MaxAll0[,var], col = c('gray90', rep('#ffa50070', length(bins)-1)), border = NA)
        points(x=-bin.wd*(1:length(bins)-1)+0.4*bin.wd, y=MeanAvg[,var], pch=20, type='b', lwd=3, col='#c21e56',cex=3)
        segments(x0=-bin.wd*(1:length(bins)-1)+0.4*bin.wd, y0=MeanAvg[,var]-StdAvg[,var], 
                 x1=-bin.wd*(1:length(bins)-1)+0.4*bin.wd, y1=MeanAvg[,var]+StdAvg[,var], col='#c21e56', lwd=2, lend=1)
        points(x=-bin.wd*(1:length(bins)-1)+0.6*bin.wd, y=MeanBest[,var], pch=18, type='b', lwd=2, col='#000080',cex=3)
        segments(x0=-bin.wd*(1:length(bins)-1)+0.6*bin.wd, y0=MeanBest[,var]-StdBest[,var], 
                 x1=-bin.wd*(1:length(bins)-1)+0.6*bin.wd, y1=MeanBest[,var]+StdBest[,var], col='#000080', lwd=2, lend=1)
        rect(-c(GeoScale[GeoScale < bin.wd*(length(bins)-1)], bin.wd*(length(bins)-1)), max(MaxAll0[,var], na.rm = T)*1.01, 
             -c(-bin.wd, GeoScale[GeoScale < bin.wd*(length(bins)-1)]), max(MaxAll0[,var], na.rm = T)*1.05, 
             border = NA, col = GeoColor[1:(sum(GeoScale < bin.wd*(length(bins)-1))+1)])
        text(x=-0.5*(c(GeoScale[GeoScale < bin.wd*(length(bins)-1)], bin.wd*(length(bins)-1))-c(-bin.wd, GeoScale[GeoScale < bin.wd*(length(bins)-1)]))-c(-bin.wd, GeoScale[GeoScale < bin.wd*(length(bins)-1)]), y=max(MaxAll0[,var], na.rm = T)*1.03, labels = GeoLabel[1:(sum(GeoScale < bin.wd*(length(bins)-1))+1)], cex = 2)
    }
    
    for (p in c("FreeLiving", "Parasitic", "HardFreeLiving", "SoftFreeLiving", "SoftParasitic", "HardParasitic")){
        MeanBest <- StdBest <- MeanAvg <- StdAvg <- MaxAll <- MinAll <- data.frame(matrix(NA, length(bins), 4), row.names = as.character(1:length(bins)))
        colnames(MeanBest) <- colnames(StdBest) <- colnames(MeanAvg) <- colnames(StdAvg) <- colnames(MaxAll) <- colnames(MinAll) <- MeasureVars
        
        # rate-state-for-nodes dataset ----
        for (i in 1:length(bins)){
            if (any(bin.idx==i)){
                if (any((RateDataBest$state2==p | RateDataBest$state4==p) & bin.idx==i)) {
                    MeanBest[i,] <- apply(RateDataBest[(RateDataBest$state2==p | RateDataBest$state4==p) & bin.idx==i,MeasureVars], 2, mean, na.rm=T)
                    StdBest[i,] <- apply(RateDataBest[(RateDataBest$state2==p | RateDataBest$state4==p) & bin.idx==i,MeasureVars], 2, sd, na.rm=T)
                }
                if (any((RateDataModAvg$state2==p | RateDataModAvg$state4==p) & bin.idx==i)) {
                    MeanAvg[i,] <- apply(RateDataModAvg[(RateDataModAvg$state2==p | RateDataModAvg$state4==p) & bin.idx==i,MeasureVars], 2, mean, na.rm=T)
                    StdAvg[i,] <- apply(RateDataModAvg[(RateDataModAvg$state2==p | RateDataModAvg$state4==p) & bin.idx==i,MeasureVars], 2, sd, na.rm=T)
                }
                MaxAll[i,] <- bind_rows(lapply(RateDataModAll, FUN=function(x){apply(x[(x$state2==p | x$state4==p) & bin.idx==i,MeasureVars], 2, max, na.rm=T)})) %>% apply(2, max, na.rm=T)
                MinAll[i,] <- bind_rows(lapply(RateDataModAll, FUN=function(x){apply(x[(x$state2==p | x$state4==p) & bin.idx==i,MeasureVars], 2, min, na.rm=T)})) %>% apply(2, min, na.rm=T)
            }
        }
        MaxAll[MaxAll==-Inf] <- NA
        MinAll[MinAll==-Inf] <- NA
        
        for (var in MeasureVars){
            plot(x=-bin.wd*(1:length(bins)-1)+0.5*bin.wd, y=MeanAvg[,var], 
                 type='n', main=ifelse(is.null(tag),paste('Rate Change Through Time (', p, ')', sep = ''),paste('Rate Change Through Time (', tag, '-', p, ')', sep='')), 
                 xlab='Million Years Before Present', ylab=var, 
                 xlim=c(-bin.wd*(length(bins)-1), bin.wd), ylim=c(min(MinAll0[,var], na.rm = T)*0.9, max(MaxAll0[,var], na.rm = T)*1.1), cex=4)
            if (! is.null(cut.point)){
                abline(v=-cut.point, col = 'gray30', lty = 2, lwd = 2)
            }
            rect(-bin.wd*(1:length(bins)-1), MinAll[,var], -bin.wd*(1:length(bins)-2), MaxAll[,var], col = c('gray90', rep('#ffa50070', length(bins)-1)), border = NA)
            points(x=-bin.wd*(1:length(bins)-1)+0.4*bin.wd, y=MeanAvg[,var], pch=20, type='b', lwd=3, col='#c21e56',cex=3)
            segments(x0=-bin.wd*(1:length(bins)-1)+0.4*bin.wd, y0=MeanAvg[,var]-StdAvg[,var], 
                     x1=-bin.wd*(1:length(bins)-1)+0.4*bin.wd, y1=MeanAvg[,var]+StdAvg[,var], col='#c21e56', lwd=2, lend=1)
            points(x=-bin.wd*(1:length(bins)-1)+0.6*bin.wd, y=MeanBest[,var], pch=18, type='b', lwd=2, col='#000080',cex=3)
            segments(x0=-bin.wd*(1:length(bins)-1)+0.6*bin.wd, y0=MeanBest[,var]-StdBest[,var], 
                     x1=-bin.wd*(1:length(bins)-1)+0.6*bin.wd, y1=MeanBest[,var]+StdBest[,var], col='#000080', lwd=2, lend=1)
            rect(-c(GeoScale[GeoScale < bin.wd*(length(bins)-1)], bin.wd*(length(bins)-1)), max(MaxAll0[,var], na.rm = T)*1.01, 
                 -c(-bin.wd, GeoScale[GeoScale < bin.wd*(length(bins)-1)]), max(MaxAll0[,var], na.rm = T)*1.05, 
                 border = NA, col = GeoColor[1:(sum(GeoScale < bin.wd*(length(bins)-1))+1)])
            text(x=-0.5*(c(GeoScale[GeoScale < bin.wd*(length(bins)-1)], bin.wd*(length(bins)-1))-c(-bin.wd, GeoScale[GeoScale < bin.wd*(length(bins)-1)]))-c(-bin.wd, GeoScale[GeoScale < bin.wd*(length(bins)-1)]), y=max(MaxAll0[,var], na.rm = T)*1.03, labels = GeoLabel[1:(sum(GeoScale < bin.wd*(length(bins)-1))+1)], cex = 2)
        }
    }
    dev.off()
    
    ## Perform state-rate tests using ANCOVA ----
    # Test for differences in evolutionary rates between time periods and ancestral states
    if (is.numeric(cut.point)){
        xdata <- RateDataModAvg[KeepNodeIDs,]
        RateModAvg01 <- lm(net.div ~ BAcut + net.div0, data = xdata)
            Yadj <- RateModAvg01$coefficients['net.div0'] * xdata$net.div0
            RateModAvg01 <- lm(net.div - Yadj ~ BAcut, data = xdata)
        RateModAvg02 <- lm(turnover ~ BAcut + turnover0, data = xdata)
            Yadj <- RateModAvg02$coefficients['turnover0'] * xdata$turnover0
            RateModAvg02 <- lm(turnover - Yadj ~ BAcut, data = xdata)
        RateModAvg03 <- lm(speciation ~ BAcut + speciation0, data = xdata)
            Yadj <- RateModAvg03$coefficients['speciation0'] * xdata$speciation0
            RateModAvg03 <- lm(speciation - Yadj ~ BAcut, data = xdata)
        RateModAvg04 <- lm(extinction ~ BAcut + extinction0, data = xdata)
            Yadj <- RateModAvg04$coefficients['extinction0'] * xdata$extinction0
            RateModAvg04 <- lm(extinction - Yadj ~ BAcut, data = xdata)
        
        xdata <- RateDataModAvg[KeepNodeIDs,]
        RateState2ModAvg01 <- lm(net.div ~ state2 + BAcut + net.div0, data = xdata)
            Yadj <- RateState2ModAvg01$coefficients['net.div0'] * xdata$net.div0
            RateState2ModAvg01 <- lm(net.div - Yadj ~ state2 + BAcut, data = xdata)
        RateState2ModAvg02 <- lm(turnover ~ state2 + BAcut + turnover0, data = xdata)
            Yadj <- RateState2ModAvg02$coefficients['turnover0'] * xdata$turnover0
            RateState2ModAvg02 <- lm(turnover - Yadj ~ state2 + BAcut, data = xdata)
        RateState2ModAvg03 <- lm(speciation ~ state2 + BAcut + speciation0, data = xdata)
            Yadj <- RateState2ModAvg03$coefficients['speciation0'] * xdata$speciation0
            RateState2ModAvg03 <- lm(speciation - Yadj ~ state2 + BAcut, data = xdata)
        RateState2ModAvg04 <- lm(extinction ~ state2 + BAcut + extinction0, data = xdata)
            Yadj <- RateState2ModAvg04$coefficients['extinction0'] * xdata$extinction0
            RateState2ModAvg04 <- lm(extinction - Yadj ~ state2 + BAcut, data = xdata)
        
        xdata <- RateDataModAvg[RateDataModAvg$time > 0 & RateDataModAvg$state4 %in% c('HardFreeLiving', 'HardParasitic'),]
        RateState4ModAvg01 <- lm(net.div ~ state4 + BAcut + net.div0, data = xdata)
            Yadj <- RateState4ModAvg01$coefficients['net.div0'] * xdata$net.div0
            RateState4ModAvg01 <- lm(net.div - Yadj ~ state4 + BAcut, data = xdata)
        RateState4ModAvg02 <- lm(turnover ~ state4 + BAcut + turnover0, data = xdata)
            Yadj <- RateState4ModAvg02$coefficients['turnover0'] * xdata$turnover0
            RateState4ModAvg02 <- lm(turnover - Yadj ~ state4 + BAcut, data = xdata)
        RateState4ModAvg03 <- lm(speciation ~ state4 + BAcut + speciation0, data = xdata)
            Yadj <- RateState4ModAvg03$coefficients['speciation0'] * xdata$speciation0
            RateState4ModAvg03 <- lm(speciation - Yadj ~ state4 + BAcut, data = xdata)
        RateState4ModAvg04 <- lm(extinction ~ state4 + BAcut + extinction0, data = xdata)
            Yadj <- RateState4ModAvg04$coefficients['extinction0'] * xdata$extinction0
            RateState4ModAvg04 <- lm(extinction - Yadj ~ state4 + BAcut, data = xdata)
        
        sink(paste('./ANCOVA report for State vs Rate on Nodes (', tag, '-model average).txt', sep = ''))
            cat('===== ANCOVA Report for All Nodes =====\n')
            cat('============ Liang Lü ============\n')
            cat('====',date(),'====\n\n')
            cat(paste('=====ANCOVA for Rate Difference between Before and After ', as.character(cut.point), ' Ma on Nodes=====\n', sep = ''))
            cat('## Net Diversification ##\n')
            print(car::Anova(RateModAvg01, type=3))
            # As long as the main effect of BAcut is significant, the LSD test is equivalent and not necessary for two-group comparison, but we still perform it for consistency with the following multiple comparisons
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateModAvg01, 'BAcut', group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            ## The Holm-adjusted LSD method is more conservative than the unadjusted LSD method, and the HSD method is more conservative than the LSD method with or without p-value adjustment, but all three methods are equivalent for two-group comparison, so we do not need to perform them here.
            # cat('\n# LSD Multiple Comparison\n')
            # RateNodeLSD <- LSD.test(RateModAvg01, 'BAcut', group=TRUE, p.adj = 'holm')
			# print(RateNodeLSD)
            # cat('\n# HSD Multiple Comparison\n')
            # RateNodeHSD <- HSD.test(RateModAvg01, 'BAcut', group=TRUE, unbalanced = TRUE)
            # print(RateNodeHSD)
            Effects <- rownames(RateNodeLSD$means)
            RateNodeTable1 <- cbind(PhyModel=strsplit(tag,'_')[[1]][2], AmbiLifeStyle=strsplit(tag,'_')[[1]][1],Model='Averaged',Effects,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Turnover ##\n')
            print(car::Anova(RateModAvg02, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateModAvg02, 'BAcut', group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable1 <- cbind(RateNodeTable1,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Speciation ##\n')
            print(car::Anova(RateModAvg03, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateModAvg03, 'BAcut', group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable1 <- cbind(RateNodeTable1,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Extinction ##\n')
            print(car::Anova(RateModAvg04, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateModAvg04, 'BAcut', group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable1 <- cbind(RateNodeTable1,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            
            cat(paste('=====ANCOVA for Rate Difference between Before and After ', as.character(cut.point), ' Ma on FreeLiving and Parasitic Nodes=====\n', sep = ''))
            cat('## Net Diversification ##\n')
            print(car::Anova(RateState2ModAvg01, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState2ModAvg01, c('state2','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            ## The Holm-adjusted LSD method is more conservative than the unadjusted LSD method, and the HSD method is more conservative than the LSD method. Here we performed the LSD test without p-value adjustment for multiple comparisons, as it is more powerful and less conservative than the other two methods. What we really expected is not simply a significant difference in any comparison, but a significant difference between post-90Ma parasitic nodes and all the other three groups (pre-90Ma parasitic nodes, pre-90Ma free-living nodes, and post-90Ma free-living nodes), where these three groups are not significantly different. We are willing to test whether our expectation remains robust under a slightly higher false positive rate for the sake of higher power to detect this specific pattern. However, we also provide reference codes that perform the LSD test with Holm p-value adjustment and the HSD test. We encourage the readers to try these two more conservative methods and copy the codes for the multiple comparisons in the following sections by themselves. In our case, the LSD test with Holm p-value adjustment and the HSD test showed little changes in significance, so the choice of multiple comparison method does not affect our conclusions.
            # cat('\n# LSD Multiple Comparison\n')
            # RateNodeLSD <- LSD.test(RateState2ModAvg01, c('state2','BAcut'), group=TRUE, p.adj = 'holm')
            # print(RateNodeLSD)
            # cat('\n# HSD Multiple Comparison\n')
            # RateNodeHSD <- HSD.test(RateState2ModAvg01, c('state2','BAcut'), group=TRUE, unbalanced = TRUE)
            # print(RateNodeHSD)
            Effects <- rownames(RateNodeLSD$means)
            RateNodeTable2 <- cbind(PhyModel=strsplit(tag,'_')[[1]][2], AmbiLifeStyle=strsplit(tag,'_')[[1]][1],Model='Averaged',Effects,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Turnover ##\n')
            print(car::Anova(RateState2ModAvg02, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState2ModAvg02, c('state2','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable2 <- cbind(RateNodeTable2,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Speciation ##\n')
            print(car::Anova(RateState2ModAvg03, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState2ModAvg03, c('state2','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable2 <- cbind(RateNodeTable2,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Extinction ##\n')
            print(car::Anova(RateState2ModAvg04, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState2ModAvg04, c('state2','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable2 <- cbind(RateNodeTable2,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            
            cat(paste('=====ANCOVA for Rate Difference between Before and After ', as.character(cut.point), ' Ma on HardFreeLiving and HardParasitic Nodes=====\n', sep = ''))
            cat('## Net Diversification ##\n')
            print(car::Anova(RateState4ModAvg01, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState4ModAvg01, c('state4','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            Effects <- rownames(RateNodeLSD$means)
            RateNodeTable3 <- cbind(PhyModel=strsplit(tag,'_')[[1]][2], AmbiLifeStyle=strsplit(tag,'_')[[1]][1],Model='Averaged',Effects,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Turnover ##\n')
            print(car::Anova(RateState4ModAvg02, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState4ModAvg02, c('state4','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable3 <- cbind(RateNodeTable3,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Speciation ##\n')
            print(car::Anova(RateState4ModAvg03, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState4ModAvg03, c('state4','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable3 <- cbind(RateNodeTable3,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Extinction ##\n')
            print(car::Anova(RateState4ModAvg04, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState4ModAvg04, c('state4','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable3 <- cbind(RateNodeTable3,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
        sink()
        
        xdata <- RateDataBest[KeepNodeIDs,]
        RateBest01 <- lm(net.div ~ BAcut + net.div0, data = xdata)
            Yadj <- RateBest01$coefficients['net.div0'] * xdata$net.div0
            RateBest01 <- lm(net.div - Yadj ~ BAcut, data = xdata)
        RateBest02 <- lm(turnover ~ BAcut + turnover0, data = xdata)
            Yadj <- RateBest02$coefficients['turnover0'] * xdata$turnover0
            RateBest02 <- lm(turnover - Yadj ~ BAcut, data = xdata)
        RateBest03 <- lm(speciation ~ BAcut + speciation0, data = xdata)
            Yadj <- RateBest03$coefficients['speciation0'] * xdata$speciation0
            RateBest03 <- lm(speciation - Yadj ~ BAcut, data = xdata)
        RateBest04 <- lm(extinction ~ BAcut + extinction0, data = xdata)
            Yadj <- RateBest04$coefficients['extinction0'] * xdata$extinction0
            RateBest04 <- lm(extinction - Yadj ~ BAcut, data = xdata)
        
        xdata <- RateDataBest[KeepNodeIDs,]
        RateState2Best01 <- lm(net.div ~ state2 + BAcut + net.div0, data = xdata)
            Yadj <- RateState2Best01$coefficients['net.div0'] * xdata$net.div0
            RateState2Best01 <- lm(net.div - Yadj ~ state2 + BAcut, data = xdata)
        RateState2Best02 <- lm(turnover ~ state2 + BAcut + turnover0, data = xdata)
            Yadj <- RateState2Best02$coefficients['turnover0'] * xdata$turnover0
            RateState2Best02 <- lm(turnover - Yadj ~ state2 + BAcut, data = xdata)
        RateState2Best03 <- lm(speciation ~ state2 + BAcut + speciation0, data = xdata)
            Yadj <- RateState2Best03$coefficients['speciation0'] * xdata$speciation0
            RateState2Best03 <- lm(speciation - Yadj ~ state2 + BAcut, data = xdata)
        RateState2Best04 <- lm(extinction ~ state2 + BAcut + extinction0, data = xdata)
            Yadj <- RateState2Best04$coefficients['extinction0'] * xdata$extinction0
            RateState2Best04 <- lm(extinction - Yadj ~ state2 + BAcut, data = xdata)
        
        xdata <- RateDataBest[RateDataBest$time > 0 & RateDataBest$state4 %in% c('HardFreeLiving', 'HardParasitic'),]
        RateState4Best01 <- lm(net.div ~ state4 + BAcut + net.div0, data = xdata)
            Yadj <- RateState4Best01$coefficients['net.div0'] * xdata$net.div0
            RateState4Best01 <- lm(net.div - Yadj ~ state4 + BAcut, data = xdata)
        RateState4Best02 <- lm(turnover ~ state4 + BAcut + turnover0, data = xdata)
            Yadj <- RateState4Best02$coefficients['turnover0'] * xdata$turnover0
            RateState4Best02 <- lm(turnover - Yadj ~ state4 + BAcut, data = xdata)
        RateState4Best03 <- lm(speciation ~ state4 + BAcut + speciation0, data = xdata)
            Yadj <- RateState4Best03$coefficients['speciation0'] * xdata$speciation0
            RateState4Best03 <- lm(speciation - Yadj ~ state4 + BAcut, data = xdata)
        RateState4Best04 <- lm(extinction ~ state4 + BAcut + extinction0, data = xdata)
            Yadj <- RateState4Best04$coefficients['extinction0'] * xdata$extinction0
            RateState4Best04 <- lm(extinction - Yadj ~ state4 + BAcut, data = xdata)
        
        sink(paste('./ANCOVA report for State vs Rate on Nodes (', tag, '-best model).txt', sep = ''))
            cat('===== ANCOVA Report for All Nodes =====\n')
            cat('============ Liang Lü ============\n')
            cat('====',date(),'====\n\n')
            cat(paste('=====ANCOVA for Rate Difference between Before and After ', as.character(cut.point), ' Ma on Nodes=====\n', sep = ''))
            cat('## Net Diversification ##\n')
            print(car::Anova(RateBest01, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateBest01, 'BAcut', group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            Effects <- rownames(RateNodeLSD$means)
            RateNodeTable4 <- cbind(PhyModel=strsplit(tag,'_')[[1]][2], AmbiLifeStyle=strsplit(tag,'_')[[1]][1],Model='Best',Effects,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Turnover ##\n')
            print(car::Anova(RateBest02, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateBest02, 'BAcut', group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable4 <- cbind(RateNodeTable4,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Speciation ##\n')
            print(car::Anova(RateBest03, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateBest03, 'BAcut', group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable4 <- cbind(RateNodeTable4,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Extinction ##\n')
            print(car::Anova(RateBest04, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateBest04, 'BAcut', group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable4 <- cbind(RateNodeTable4,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            
            cat(paste('=====ANCOVA for Rate Difference between Before and After ', as.character(cut.point), ' Ma on FreeLiving and Parasitic Nodes=====\n', sep = ''))
            cat('## Net Diversification ##\n')
            print(car::Anova(RateState2Best01, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState2Best01, c('state2','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            Effects <- rownames(RateNodeLSD$means)
            RateNodeTable5 <- cbind(PhyModel=strsplit(tag,'_')[[1]][2], AmbiLifeStyle=strsplit(tag,'_')[[1]][1],Model='Best',Effects,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Turnover ##\n')
            print(car::Anova(RateState2Best02, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState2Best02, c('state2','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable5 <- cbind(RateNodeTable5,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Speciation ##\n')
            print(car::Anova(RateState2Best03, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState2Best03, c('state2','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable5 <- cbind(RateNodeTable5,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Extinction ##\n')
            print(car::Anova(RateState2Best04, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState2Best04, c('state2','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable5 <- cbind(RateNodeTable5,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            
            cat(paste('=====ANCOVA for Rate Difference between Before and After ', as.character(cut.point), ' Ma on HardFreeLiving and HardParasitic Nodes=====\n', sep = ''))
            cat('## Net Diversification ##\n')
            print(car::Anova(RateState4Best01, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState4Best01, c('state4','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            Effects <- rownames(RateNodeLSD$means)
            RateNodeTable6 <- cbind(PhyModel=strsplit(tag,'_')[[1]][2], AmbiLifeStyle=strsplit(tag,'_')[[1]][1],Model='Best',Effects,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Turnover ##\n')
            print(car::Anova(RateState4Best02, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState4Best02, c('state4','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable6 <- cbind(RateNodeTable6,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Speciation ##\n')
            print(car::Anova(RateState4Best03, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState4Best03, c('state4','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable6 <- cbind(RateNodeTable6,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
            cat('## Extinction ##\n')
            print(car::Anova(RateState4Best04, type=3))
            cat('\n# LSD Multiple Comparison\n')
            RateNodeLSD <- LSD.test(RateState4Best04, c('state4','BAcut'), group=TRUE, p.adj = 'none')
            print(RateNodeLSD)
            RateNodeTable6 <- cbind(RateNodeTable6,RateNodeLSD$means[Effects,1:2],groups=RateNodeLSD$groups[Effects,'groups'])
            cat('=========================\n\n')
        sink()
        
        RateNodeTable <- rbind(RateNodeTable1,
                               RateNodeTable2,
                               RateNodeTable3,
                               RateNodeTable4,
                               RateNodeTable5,
                               RateNodeTable6)
        
        write.csv(RateNodeTable, file = paste('./ANCOVA report for State vs Rate on Nodes (', tag, ').csv', sep = ''), row.names = F)
        
    }else{
        xdata <- RateDataModAvg[KeepNodeIDs,]
        RateState2ModAvg01 <- lm(net.div ~ state2 + net.div0, data = xdata)
            Yadj <- RateState2ModAvg01$coefficients['net.div0'] * xdata$net.div0
            RateState2ModAvg01 <- lm(net.div - Yadj ~ state2, data = xdata)
        RateState2ModAvg02 <- lm(turnover ~ state2 + turnover0, data = xdata)
            Yadj <- RateState2ModAvg02$coefficients['turnover0'] * xdata$turnover0
            RateState2ModAvg02 <- lm(turnover - Yadj ~ state2, data = xdata)
        RateState2ModAvg03 <- lm(speciation ~ state2 + speciation0, data = xdata)
            Yadj <- RateState2ModAvg03$coefficients['speciation0'] * xdata$speciation0
            RateState2ModAvg03 <- lm(speciation - Yadj ~ state2, data = xdata)
        RateState2ModAvg04 <- lm(extinction ~ state2 + extinction0, data = xdata)
            Yadj <- RateState2ModAvg04$coefficients['extinction0'] * xdata$extinction0
            RateState2ModAvg04 <- lm(extinction - Yadj ~ state2, data = xdata)
        
        xdata <- RateDataModAvg[RateDataModAvg$time > 0 & RateDataModAvg$state4 %in% c('HardFreeLiving', 'HardParasitic'),]
        RateState4ModAvg01 <- lm(net.div ~ state4 + net.div0, data = xdata)
            Yadj <- RateState4ModAvg01$coefficients['net.div0'] * xdata$net.div0
            RateState4ModAvg01 <- lm(net.div - Yadj ~ state4, data = xdata)
        RateState4ModAvg02 <- lm(turnover ~ state4 + turnover0, data = xdata)
            Yadj <- RateState4ModAvg02$coefficients['turnover0'] * xdata$turnover0
            RateState4ModAvg02 <- lm(turnover - Yadj ~ state4, data = xdata)
        RateState4ModAvg03 <- lm(speciation ~ state4 + speciation0, data = xdata)
            Yadj <- RateState4ModAvg03$coefficients['speciation0'] * xdata$speciation0
            RateState4ModAvg03 <- lm(speciation - Yadj ~ state4, data = xdata)
        RateState4ModAvg04 <- lm(extinction ~ state4 + extinction0, data = xdata)
            Yadj <- RateState4ModAvg04$coefficients['extinction0'] * xdata$extinction0
            RateState4ModAvg04 <- lm(extinction - Yadj ~ state4, data = xdata)
        
        sink(paste('./ANCOVA report for State vs Rate on Nodes (', tag, '-model average).txt', sep = ''))
            cat('===== ANCOVA Report for Nodes =====\n')
            cat('============ Liang Lü ============\n')
            cat('====',date(),'====\n\n')
            
            cat('=====ANCOVA for Rate Difference between FreeLiving and Parasitic Nodes=====\n')
            cat('## Net Diversification ##\n')
            print(car::Anova(RateState2ModAvg01, type=3))
            cat('=========================\n\n')
            cat('## Turnover ##\n')
            print(car::Anova(RateState2ModAvg02, type=3))
            cat('=========================\n\n')
            cat('## Speciation ##\n')
            print(car::Anova(RateState2ModAvg03, type=3))
            cat('=========================\n\n')
            cat('## Extinction ##\n')
            print(car::Anova(RateState2ModAvg04, type=3))
            cat('=========================\n\n')
            
            cat('=====ANCOVA for Rate Difference between HardFreeLiving and HardParasitic Nodes=====\n')
            cat('## Net Diversification ##\n')
            print(car::Anova(RateState4ModAvg01, type=3))
            cat('=========================\n\n')
            cat('## Turnover ##\n')
            print(car::Anova(RateState4ModAvg02, type=3))
            cat('=========================\n\n')
            cat('## Speciation ##\n')
            print(car::Anova(RateState4ModAvg03, type=3))
            cat('=========================\n\n')
            cat('## Extinction ##\n')
            print(car::Anova(RateState4ModAvg04, type=3))
            cat('=========================\n\n')
        sink()
        
        xdata <- RateDataBest[KeepNodeIDs,]
        RateState2Best01 <- lm(net.div ~ state2 + net.div0, data = xdata)
            Yadj <- RateState2Best01$coefficients['net.div0'] * xdata$net.div0
            RateState2Best01 <- lm(net.div - Yadj ~ state2, data = xdata)
        RateState2Best02 <- lm(turnover ~ state2 + turnover0, data = xdata)
            Yadj <- RateState2Best02$coefficients['turnover0'] * xdata$turnover0
            RateState2Best02 <- lm(turnover - Yadj ~ state2, data = xdata)
        RateState2Best03 <- lm(speciation ~ state2 + speciation0, data = xdata)
            Yadj <- RateState2Best03$coefficients['speciation0'] * xdata$speciation0
            RateState2Best03 <- lm(speciation - Yadj ~ state2, data = xdata)
        RateState2Best04 <- lm(extinction ~ state2 + extinction0, data = xdata)
            Yadj <- RateState2Best04$coefficients['extinction0'] * xdata$extinction0
            RateState2Best04 <- lm(extinction - Yadj ~ state2, data = xdata)
        
        xdata <- RateDataBest[RateDataBest$time > 0 & RateDataBest$state4 %in% c('HardFreeLiving', 'HardParasitic'),]
        RateState4Best01 <- lm(net.div ~ state4 + net.div0, data = xdata)
            Yadj <- RateState4Best01$coefficients['net.div0'] * xdata$net.div0
            RateState4Best01 <- lm(net.div - Yadj ~ state4, data = xdata)
        RateState4Best02 <- lm(turnover ~ state4 + turnover0, data = xdata)
            Yadj <- RateState4Best02$coefficients['turnover0'] * xdata$turnover0
            RateState4Best02 <- lm(turnover - Yadj ~ state4, data = xdata)
        RateState4Best03 <- lm(speciation ~ state4 + speciation0, data = xdata)
            Yadj <- RateState4Best03$coefficients['speciation0'] * xdata$speciation0
            RateState4Best03 <- lm(speciation - Yadj ~ state4, data = xdata)
        RateState4Best04 <- lm(extinction ~ state4 + extinction0, data = xdata)
            Yadj <- RateState4Best04$coefficients['extinction0'] * xdata$extinction0
            RateState4Best04 <- lm(extinction - Yadj ~ state4, data = xdata)
        
        sink(paste('./ANCOVA report for State vs Rate on Nodes (', tag, '-best model).txt', sep = ''))
            cat('===== ANCOVA Report for Nodes =====\n')
            cat('============ Liang Lü ============\n')
            cat('====',date(),'====\n\n')
            
            cat('=====ANCOVA for Rate Difference between FreeLiving and Parasitic Nodes=====\n')
            cat('## Net Diversification ##\n')
            print(car::Anova(RateState2Best01, type=3))
            cat('=========================\n\n')
            cat('## Turnover ##\n')
            print(car::Anova(RateState2Best02, type=3))
            cat('=========================\n\n')
            cat('## Speciation ##\n')
            print(car::Anova(RateState2Best03, type=3))
            cat('=========================\n\n')
            cat('## Extinction ##\n')
            print(car::Anova(RateState2Best04, type=3))
            cat('=========================\n\n')
            
            cat('=====ANCOVA for Rate Difference between HardFreeLiving and HardParasitic Nodes=====\n')
            cat('## Net Diversification ##\n')
            print(car::Anova(RateState4Best01, type=3))
            cat('=========================\n\n')
            cat('## Turnover ##\n')
            print(car::Anova(RateState4Best02, type=3))
            cat('=========================\n\n')
            cat('## Speciation ##\n')
            print(car::Anova(RateState4Best03, type=3))
            cat('=========================\n\n')
            cat('## Extinction ##\n')
            print(car::Anova(RateState4Best04, type=3))
            cat('=========================\n\n')
        sink()
    }
    
    # Perform phylogenetic ANOVA for tip rates ----
    # Test for differences in rates between free-living and parasitic tips accounting for phylogeny
    
    grp <<- as.factor(RateDataBest[KeepTipIDs, 'state2'])
    dat01 <<- RateDataModAvg[KeepTipIDs, 'net.div']
    dat02 <<- RateDataModAvg[KeepTipIDs, 'turnover']
    dat03 <<- RateDataModAvg[KeepTipIDs, 'speciation']
    dat04 <<- RateDataModAvg[KeepTipIDs, 'extinction']
    dat05 <<- RateDataBest[KeepTipIDs, 'net.div']
    dat06 <<- RateDataBest[KeepTipIDs, 'turnover']
    dat07 <<- RateDataBest[KeepTipIDs, 'speciation']
    dat08 <<- RateDataBest[KeepTipIDs, 'extinction']
    
    
    names(grp) <<- names(dat01) <<- names(dat02) <<- names(dat03) <<- names(dat04) <<- names(dat05) <<- names(dat06) <<- names(dat07) <<- names(dat08) <<- KeepTips
    
    sink(paste('./PhyloANOVA for Comparison of Tip Rates (', tag, ').txt', sep = ''))
        cat('===== ANOVA Report for Tips ======\n')
        cat('============ Liang Lü ============\n')
        cat('====',date(),'====\n\n')
        
        cat('====Averaged model====\n\n')
        cat('## Net Diversification ##\n')
        (RateStateTip <- attributes(geiger::aov.phylo(dat01 ~ grp, Phy, nsim = 1000, test = 'Wilks'))$summary)
        RateTipANOVA1 <- RateStateTip[,-2]
        cat('=========================\n\n')
        cat('## Turnover ##\n')
        (RateStateTip <- attributes(geiger::aov.phylo(dat02 ~ grp, Phy, nsim = 1000, test = 'Wilks'))$summary)
        RateTipANOVA1 <- rbind(RateTipANOVA1, RateStateTip[,-2])
        cat('=========================\n\n')
        cat('## Speciation ##\n')
        (RateStateTip <- attributes(geiger::aov.phylo(dat03 ~ grp, Phy, nsim = 1000, test = 'Wilks'))$summary)
        RateTipANOVA1 <- rbind(RateTipANOVA1, RateStateTip[,-2])
        cat('## Extinction ##\n')
        (RateStateTip <- attributes(geiger::aov.phylo(dat04 ~ grp, Phy, nsim = 1000, test = 'Wilks'))$summary)
        RateTipANOVA1 <- rbind(RateTipANOVA1, RateStateTip[,-2])
        cat('=========================\n\n')
        
        RateTipANOVA1 <- cbind(PhyModel=c(strsplit(tag,'_')[[1]][2], rep(NA, 7)), 
                               AmbiLifeStyle=c(strsplit(tag,'_')[[1]][1], rep(NA, 7)), 
                               Rate=c('net.div', NA, 'turnover', NA, 'speciation', NA, 'extinction', NA),
                               SoV=rep(c('Lifestyle','Residuals'), 4),
                               RateTipANOVA1)
        
        cat('====Best-fitting model====\n\n')
        cat('## Net Diversification ##\n')
        (RateStateTip <- attributes(geiger::aov.phylo(dat05 ~ grp, Phy, nsim = 1000, test = 'Wilks'))$summary)
        RateTipANOVA2 <- RateStateTip[,-1:-2]
        cat('=========================\n\n')
        cat('## Turnover ##\n')
        (RateStateTip <- attributes(geiger::aov.phylo(dat06 ~ grp, Phy, nsim = 1000, test = 'Wilks'))$summary)
        RateTipANOVA2 <- rbind(RateTipANOVA2, RateStateTip[,-1:-2])
        cat('=========================\n\n')
        cat('## Speciation ##\n')
        (RateStateTip <- attributes(geiger::aov.phylo(dat07 ~ grp, Phy, nsim = 1000, test = 'Wilks'))$summary)
        RateTipANOVA2 <- rbind(RateTipANOVA2, RateStateTip[,-1:-2])
        cat('=========================\n\n')
        cat('## Extinction ##\n')
        (RateStateTip <- attributes(geiger::aov.phylo(dat08 ~ grp, Phy, nsim = 1000, test = 'Wilks'))$summary)
        RateTipANOVA2 <- rbind(RateTipANOVA2, RateStateTip[,-1:-2])
        cat('=========================\n\n')
        
    sink()
    
    RateTipANOVA <- cbind(RateTipANOVA1, RateTipANOVA2)
    
    write.csv(RateTipANOVA, file = paste('./PhyloANOVA for Comparison of Tip Rates (', tag, ').csv', sep = ''), row.names = F)
    
    rm(grp, dat01, dat02, dat03, dat04, dat05, dat06, dat07, dat08, envir = .GlobalEnv)
    
    ## Plot node and tip rate distributions ----
    if (is.numeric(cut.point)){
        # data preparation
        StateVars <- c("All", "FreeLiving", "Parasitic", "HardFreeLiving", "SoftFreeLiving", "SoftParasitic", "HardParasitic")
        GroupVars <- c('stateX', 'Mod', 'BAcut')
        
        # Node rate distribution
        Data <- rbind(RateDataModAvg[KeepNodeIDs,c('BAcut',MeasureVars)], 
                      RateDataBest[KeepNodeIDs,c('BAcut',MeasureVars)])
        Data$stateX <- 'All'
        Data <- rbind(RateDataModAvg[KeepNodeIDs,c('BAcut',MeasureVars,'state2')], 
                      RateDataBest[KeepNodeIDs,c('BAcut',MeasureVars,'state2')]) %>% plyr::rename(., c('state2' = 'stateX')) %>% rbind(Data,.)
        Data <- rbind(RateDataModAvg[KeepNodeIDs,c('BAcut',MeasureVars,'state4')], 
                      RateDataBest[KeepNodeIDs,c('BAcut',MeasureVars,'state4')]) %>% plyr::rename(., c('state4' = 'stateX')) %>% rbind(Data,.)
        Data$Mod <- rep(c('ModelAvg', 'BestModel'), each = length(KeepNodeIDs),times = 3)
        # Data <- Data[Data$stateX!='SoftFreeLiving' & Data$stateX!='SoftParasitic',]
        
        XMin <- apply(Data[, MeasureVars], 2, min, na.rm=TRUE)
        XMax <- apply(Data[, MeasureVars], 2, max, na.rm=TRUE)
        
        # plot
        plt <- NULL
        
        for (sat in StateVars){
            for (var in MeasureVars) {
                DataX <- Data[Data$stateX==sat, c(GroupVars, var)]
                colnames(DataX)[4] <- 'MeaVal'
                p <- ggplot() + 
                    geom_histogram(data=filter(DataX, Mod=='ModelAvg'), 
                                   bins = 20, position = 'dodge', 
                                   aes(x=MeaVal, fill=Mod, alpha=BAcut), color='#ffffff') + 
                    geom_histogram(data=filter(DataX, Mod=='BestModel'), 
                                   bins = 20, position = 'dodge', 
                                   aes(x=MeaVal, y=-after_stat(count), fill=Mod, alpha=BAcut), color='#ffffff') + 
                    # scale_color_manual(values=c('#000080','#c21e56')) + 
                    scale_fill_manual(values=c('#000080','#c21e56'), breaks = c('BestModel', 'ModelAvg')) + 
                    scale_alpha_manual(values=c(1, 0.3), breaks=c('After', 'Before')) + 
                    geom_hline(yintercept = 0, color='gray30') + 
                    # facet_grid(BAcut ~.) + 
                    theme_bw() + 
                    theme(panel.grid = element_blank(), 
                          plot.title = element_text(size = 20), 
                          axis.title = element_text(size = 18), 
                          axis.text = element_text(size = 18)) + 
                    scale_x_continuous(limits = c(XMin[var]*0.8, XMax[var]*1.2)) + 
                    xlab(var) + 
                    ylab('Frequency') + 
                    ggtitle(paste(sat, var, sep = ': '))
                plt <- append(plt, list(p))
                rm(p)
            }
        }
        
        # Tip rate distribution
        Data <- rbind(RateDataModAvg[KeepTipIDs,c('BAcut',MeasureVars,'state2')], 
                      RateDataBest[KeepTipIDs,c('BAcut',MeasureVars,'state2')]) %>% plyr::rename(., c('state2' = 'stateX'))
        Data$Mod <- rep(c('ModelAvg', 'BestModel'), each = length(KeepTipIDs))
        
        XMin <- apply(Data[, MeasureVars], 2, min, na.rm=TRUE)
        XMax <- apply(Data[, MeasureVars], 2, max, na.rm=TRUE)
        
        for (var in MeasureVars) {
            DataX <- Data[, c(GroupVars, var)]
            colnames(DataX)[4] <- 'MeaVal'
            p <- ggplot() + 
                geom_histogram(data=filter(DataX, Mod=='ModelAvg'), 
                               bins = 20, position = 'dodge', 
                               aes(x=MeaVal, fill=Mod, alpha=stateX), color='#ffffff') + 
                geom_histogram(data=filter(DataX, Mod=='BestModel'), 
                               bins = 20, position = 'dodge', 
                               aes(x=MeaVal, y=-after_stat(count), fill=Mod, alpha=stateX), color='#ffffff') + 
                # scale_color_manual(values=c('#000080','#c21e56')) + 
                scale_fill_manual(values=c('#000080','#c21e56'), breaks = c('BestModel', 'ModelAvg')) + 
                scale_alpha_manual(values=c(1, 0.3), breaks=c('Parasitic', 'FreeLiving')) +
                geom_hline(yintercept = 0, color='gray30') + 
                # facet_grid(BAcut ~.) + 
                theme_bw() + 
                theme(panel.grid = element_blank(), 
                      plot.title = element_text(size = 20), 
                      axis.title = element_text(size = 18), 
                      axis.text = element_text(size = 18)) + 
                scale_x_continuous(limits = c(XMin[var]*0.8, XMax[var]*1.2)) + 
                xlab(var) + 
                ylab('Frequency') + 
                ggtitle(paste('Tips', var, sep = ': '))
            plt <- append(plt, list(p))
            rm(p)
        }

    }else{
        
        # data preparation
        StateVars <- c("All", "FreeLiving", "Parasitic", "HardFreeLiving", "SoftFreeLiving", "SoftParasitic", "HardParasitic")
        GroupVars <- c('stateX', 'Mod')
        
        # Node rate distribution
        Data <- rbind(RateDataModAvg[KeepNodeIDs,c(MeasureVars)], 
                      RateDataBest[KeepNodeIDs,c(MeasureVars)])
        Data$stateX <- 'All'
        Data <- rbind(RateDataModAvg[KeepNodeIDs,c(MeasureVars,'state2')], 
                      RateDataBest[KeepNodeIDs,c(MeasureVars,'state2')]) %>% plyr::rename(., c('state2' = 'stateX')) %>% rbind(Data,.)
        Data <- rbind(RateDataModAvg[KeepNodeIDs,c(MeasureVars,'state4')], 
                      RateDataBest[KeepNodeIDs,c(MeasureVars,'state4')]) %>% plyr::rename(., c('state4' = 'stateX')) %>% rbind(Data,.)
        Data$Mod <- rep(c('ModelAvg', 'BestModel'), each = length(KeepNodeIDs),times = 3)
        
        XMin <- apply(Data[, MeasureVars], 2, min, na.rm=TRUE)
        XMax <- apply(Data[, MeasureVars], 2, max, na.rm=TRUE)
        
        # plot
        plt <- NULL
        
        for (sat in StateVars){
            for (var in MeasureVars) {
                DataX <- Data[Data$stateX==sat, c(GroupVars, var)]
                colnames(DataX)[3] <- 'MeaVal'
                p <- ggplot() + 
                    geom_histogram(data=filter(DataX, Mod=='ModelAvg'), 
                                   bins = 20, 
                                   aes(x=MeaVal, fill=Mod), color='#ffffff') + 
                    geom_histogram(data=filter(DataX, Mod=='BestModel'), 
                                   bins = 20, 
                                   aes(x=MeaVal, y=-after_stat(count), fill=Mod), color='#ffffff') + 
                    # scale_color_manual(values=c('#000080','#c21e56')) + 
                    scale_fill_manual(values=c('#000080','#c21e56'), breaks = c('BestModel', 'ModelAvg')) + 
                    geom_hline(yintercept = 0, color='gray30') + 
                    # facet_grid(BAcut ~.) + 
                    theme_bw() + 
                    theme(panel.grid = element_blank(), 
                          plot.title = element_text(size = 20), 
                          axis.title = element_text(size = 18), 
                          axis.text = element_text(size = 18)) + 
                    scale_x_continuous(limits = c(XMin[var]*0.8, XMax[var]*1.2)) + 
                    xlab(var) + 
                    ylab('Frequency') + 
                    ggtitle(paste(sat, var, sep = ': '))
                plt <- append(plt, list(p))
                rm(p)
            }
        }
        
        # Tip rate distribution
        Data <- rbind(RateDataModAvg[KeepTipIDs,c(MeasureVars,'state2')], 
                      RateDataBest[KeepTipIDs,c(MeasureVars,'state2')]) %>% plyr::rename(., c('state2' = 'stateX'))
        Data$Mod <- rep(c('ModelAvg', 'BestModel'), each = length(KeepTipIDs))
        
        XMin <- apply(Data[, MeasureVars], 2, min, na.rm=TRUE)
        XMax <- apply(Data[, MeasureVars], 2, max, na.rm=TRUE)
        
        for (var in MeasureVars) {
            DataX <- Data[, c(GroupVars, var)]
            colnames(DataX)[3] <- 'MeaVal'
            p <- ggplot() + 
                geom_histogram(data=filter(DataX, stateX=='Parasitic'), 
                               bins = 20, position = 'dodge', 
                               aes(x=MeaVal, fill=Mod), color='#ffffff') + 
                geom_histogram(data=filter(DataX, stateX=='FreeLiving'), 
                               bins = 20, position = 'dodge', 
                               aes(x=MeaVal, y=-after_stat(count), fill=Mod), color='#ffffff') + 
                # scale_color_manual(values=c('#000080','#c21e56')) + 
                scale_fill_manual(values=c('#000080','#c21e56'), breaks = c('BestModel', 'ModelAvg')) + 
                # scale_alpha_manual(values=c(1, 0.6), breaks=c('Parasitic', 'FreeLiving')) + 
                geom_hline(yintercept = 0, color='gray30') + 
                # facet_grid(BAcut ~.) + 
                theme_bw() + 
                theme(panel.grid = element_blank(), 
                      plot.title = element_text(size = 20), 
                      axis.title = element_text(size = 18), 
                      axis.text = element_text(size = 18)) + 
                scale_x_continuous(limits = c(XMin[var]*0.8, XMax[var]*1.2)) + 
                xlab(var) + 
                ylab('Frequency') + 
                ggtitle(paste('Tips', var, sep = ': '))
            plt <- append(plt, list(p))
            rm(p)
        }
    }
    
    # Combine all plots and save to file ----
    if (length(plt) == 32){
        
        mplt<-(plt[[1]]+plt[[2]]+plt[[3]]+plt[[4]]+
               plt[[5]]+plt[[6]]+plt[[7]]+plt[[8]]+
               plt[[9]]+plt[[10]]+plt[[11]]+plt[[12]]+
               plt[[13]]+plt[[14]]+plt[[15]]+plt[[16]]+
               plt[[17]]+plt[[18]]+plt[[19]]+plt[[20]]+
               plt[[21]]+plt[[22]]+plt[[23]]+plt[[24]]+
               plt[[25]]+plt[[26]]+plt[[27]]+plt[[28]]+
               plt[[29]]+plt[[30]]+plt[[31]]+plt[[32]])+
            patchwork::plot_layout(ncol = 4)
        
        ggsave(file = paste('./Rate Distributions Histogram for Nodes and Tips - ', tag, '.svg', sep = ''), plot = mplt, width = 60, height = 64, units = 'in', limitsize = FALSE)
        
    }else{
        
        print(paste('Insufficient plt elements: ', as.character(length(plt)), sep=''))
        
    }
    
    # Reset warning level
    options(warn = 0)
}
