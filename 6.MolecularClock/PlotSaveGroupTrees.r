# PlotSaveGroupTrees.r
# R script for visualizing and saving grouped phylogenetic trees from MCMC molecular clock analyses.
# This function processes clustered MCMC trees, computes summary statistics, and generates
# publication-quality plots with uncertainty intervals and taxonomic annotations.
#
# Function: PlotSaveGroupTrees
# Purpose: Generate averaged and centroid trees for each cluster, compute HPD intervals,
#          and create SVG plots with geological time scale and taxonomic coloring.
# Parameters:
#   good_trees: List of treedata objects from MCMCtree output
#   cluster: Vector assigning each tree to a cluster
#   centpoints: Indices of centroid trees for each cluster
#   fileprefix: Prefix for output file names

PlotSaveGroupTrees <- function(good_trees, cluster, centpoints, fileprefix){
    # Extract branching times from all trees for summary statistics
    BrTime <- Reduce(rbind, lapply(good_trees, function(x){branching.times(x@phylo)}))

    # Compute mean of MCMC-means for each cluster and node
    MuMedians <- t(sapply(1:max(cluster),
                          FUN=function(x){
                              apply(BrTime[cluster==x,], MARGIN = 2, mean)
                          }))
    colnames(MuMedians) <- as.character(seq(1, good_trees[[1]]@phylo$Nnode)+Ntip(good_trees[[1]]@phylo))

    # Compute standard error of MCMC-means for each cluster and node
    SEMedians <- t(sapply(1:max(cluster),
                          FUN=function(x){
                              apply(BrTime[cluster==x,], MARGIN = 2, sd)
                          }))
    colnames(SEMedians) <- as.character(seq(1, good_trees[[1]]@phylo$Nnode)+Ntip(good_trees[[1]]@phylo))

    # Helper function to convert branching times back to edge lengths
    time2edgelength <- function(branchingtimes, edges, rooted = F){
        if (rooted){
            if (length(branchingtimes) == 0.5 * (nrow(edges) - 1)){
                branchingtimes <- c(branchingtimes[1], branchingtimes)
            }
            times <- c(rep(0, 0.5 * (nrow(edges) + 1)), branchingtimes)
        }else{
            times <- c(rep(0, 0.5 * (nrow(edges) + 2)), branchingtimes)
        }
        new.edges <- apply(edges, 1, FUN = function(x){times[x[1]] - times[x[2]]})
        return(new.edges)
    }

    # Initialize lists to store averaged and centroid trees for each cluster
    GrpAvgTrees <- NULL
    GrpCntTrees <- NULL

    # Process each cluster (skip cluster 1 as it may be outliers)
    for (i in sort(unique(cluster))[-1]){
        # Extract trees belonging to current cluster
        GrpTrees <- good_trees[cluster==i]

        # Compute 95% HPD intervals across all trees in cluster
        # Lower bounds
        HPD95.L <- Reduce(rbind, lapply(GrpTrees, function(x){Reduce(rbind, x@data$`0.95HPD`)[,1]}))
        HPD95.L.min <- apply(HPD95.L, 2, min)
        HPD95.L.max <- apply(HPD95.L, 2, max)
        # Upper bounds
        HPD95.U <- Reduce(rbind, lapply(GrpTrees, function(x){Reduce(rbind, x@data$`0.95HPD`)[,2]}))
        HPD95.U.min <- apply(HPD95.U, 2, min)
        HPD95.U.max <- apply(HPD95.U, 2, max)

        # Create strict and relaxed HPD intervals
        HPD95.str <- apply(cbind(HPD95.L.max, HPD95.U.min), 1, list)
        HPD95.str <- lapply(HPD95.str, function(x){unlist(x, use.names = F)})
        HPD95.rlx <- apply(cbind(HPD95.L.min, HPD95.U.max), 1, list)
        HPD95.rlx <- lapply(HPD95.rlx, function(x){unlist(x, use.names = F)})

        # Similarly for 95% ETI (Equal-tailed intervals)
        ETI95.L <- Reduce(rbind, lapply(GrpTrees, function(x){Reduce(rbind, x@data$`0.95ETI`)[,1]}))
        ETI95.L.min <- apply(ETI95.L, 2, min)
        ETI95.L.max <- apply(ETI95.L, 2, max)
        ETI95.U <- Reduce(rbind, lapply(GrpTrees, function(x){Reduce(rbind, x@data$`0.95ETI`)[,2]}))
        ETI95.U.min <- apply(ETI95.U, 2, min)
        ETI95.U.max <- apply(ETI95.U, 2, max)
        ETI95.str <- apply(cbind(ETI95.L.max, ETI95.U.min), 1, list)
        ETI95.str <- lapply(ETI95.str, function(x){unlist(x, use.names = F)})
        ETI95.rlx <- apply(cbind(ETI95.L.min, ETI95.U.max), 1, list)
        ETI95.rlx <- lapply(ETI95.rlx, function(x){unlist(x, use.names = F)})

        # Create averaged tree using centroid tree as template
        tredat <- GrpTrees[[centpoints[i]]]
        tredat@phylo$edge.length <- time2edgelength(MuMedians[i,], tredat@phylo$edge, rooted = T)
        # Clear existing data columns and add new summary statistics
        tredat@data[,1:11] <- NULL
        tredat@data$node <- GrpTrees[[centpoints[i]]]@data$node
        tredat@data$EPSILON <- SEMedians[i,-1][GrpTrees[[centpoints[i]]]@data$node]
        tredat@data$`0.95HPD_STR` <- HPD95.str
        tredat@data$`0.95HPD_RLX` <- HPD95.rlx
        tredat@data$`0.95ETI_STR` <- ETI95.str
        tredat@data$`0.95ETI_RLX` <- ETI95.rlx
        tredat@file <- paste('AverageTree_Grp_',i,'.tre',sep = '')
        tredat@treetext <- ''
        GrpAvgTrees[[i]] <- tredat

        # Write averaged tree to BEAST format
        write.beast(treeio::as.treedata(tibble::as_tibble(tredat)), file = paste(fileprefix,'-AverageTree_Grp_',i,'.tre',sep = ''), translate = F, tree.name = paste('AvgGrp_',i,sep = ''))

        # Create centroid tree with summary statistics
        GrpTrees[[centpoints[i]]]@data$MU <- MuMedians[i,-1][GrpTrees[[centpoints[i]]]@data$node]
        GrpTrees[[centpoints[i]]]@data$EPSILON <- SEMedians[i,-1][GrpTrees[[centpoints[i]]]@data$node]
        GrpTrees[[centpoints[i]]]@data$`0.95HPD_STR` <- HPD95.str
        GrpTrees[[centpoints[i]]]@data$`0.95HPD_RLX` <- HPD95.rlx
        GrpTrees[[centpoints[i]]]@data$`0.95ETI_STR` <- ETI95.str
        GrpTrees[[centpoints[i]]]@data$`0.95ETI_RLX` <- ETI95.rlx
        GrpCntTrees[[i]] <- GrpTrees[[centpoints[i]]]

        # Write centroid tree to BEAST format
        write.beast(treeio::as.treedata(tibble::as_tibble(GrpTrees[[centpoints[i]]])), file = paste(fileprefix,'-CentroidTree_Grp_',i,'.tre',sep = ''), translate = F, tree.name = paste('CentroidGrp_',i,sep = ''))
    }
    
    # Define color palette for taxonomic groups (used in edge coloring)
    ColorX <- c('#cc8cbc','#95a7a4','#948f6c','#609c70','#f05d4a','#aedbe1','#90447c','#002468','#194f1d','#5d1d07','#edc720','#dccaa8','#b59af2','#b56710','#dc6626','#42e476','#d0aca2','#6fc8d8','#a68a8d','#32a04d','#d7c1c0','#9d0000','#7c664e','#1d4d66','#8b4c39','#2982a2','#ba8547','#846671','#5f6b6b','#584f31','#4525b1','#cdcd00','#c4bda5','#ebd261','#ac9f77','#627f89','#9d5d8c','#aea88f','#146d1e','#71542b','#bfbfb4','#8b7a50','#1cbcd5','#bab590','#553e2d','#92310f','#00868b','#591618','#978466','#a4915c','#e29873','#36648b','#693e5a','#327570','#8c571b','#ba6452','#a787a1','#b5a49e','#ba3000','#646d4f','#d21111','#7c9d9e','#3c3520','#be9e1e','#c380d9','#579261','#90988e','#7f7f66','#ffdcca','#917a63','#b6cd2b','#bda683') 
    
    # Function to assign colors to tree edges based on taxonomic classification
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
        m <- mrca(tree)
        for (i in 1:ncat1){
            ed <- unique(as.vector(m[tips[cat==CAT1[i,2]],tips[cat==CAT1[i,2]]]))
            Colors[tree$edge[,2] %in% ed] <- CAT1[i,1]
            BarTxtPos[i,4] <- CAT1[i,1]
            BarTxtPos[i,3] <- max(SppYPos[tips[cat==CAT1[i,2]]]) 
            BarTxtPos[i,2] <- BarTxtPos[i,3] + 0.1
            BarTxtPos[i,1] <- min(SppYPos[tips[cat==CAT1[i,2]]]) - 0.1
        }
        return(list(Colors, BarTxtPos, CAT0, CAT1[,2]))
    }
    
    # Calculate maximum branching times for alignment across plots
    xxx <- unlist(lapply(GrpCntTrees, function(x){max(branching.times(x@phylo))}))
    MaxBrTimeCompan <- max(xxx)-xxx
    xxx <- unlist(lapply(GrpAvgTrees, function(x){max(branching.times(x@phylo))}))
    MaxBrTimeCompan <- list(MaxBrTimeCompan, max(xxx)-xxx)
    rm(xxx)
    
    # Generate SVG plot for centroid trees of all clusters
    svg(filename = paste(fileprefix,'-MarkedGroupCentroidTrees.svg',sep=''), width = 40, height = 80, pointsize = 12)
    layout(mat = matrix(c(1:length(GrpCntTrees)),nrow = length(GrpCntTrees)))
    for (i in 1:length(GrpCntTrees)){
        # Ladderize tree for consistent orientation
        lad.tree <- ladderize(GrpCntTrees[[i]]@phylo, right = T)
        if (i == 1){
            EdgeBarTxtColors <- EdgeColors(lad.tree, ClassShort)
        }
        # Calculate positions for HPD bars
        nodebar_pos <- max(branching.times(lad.tree)) - Reduce(rbind, GrpCntTrees[[i]]@data$`0.95HPD`[order(as.numeric(GrpCntTrees[[i]]@data$node))])
        nodebar_pos <- cbind(nodebar_pos, node.height(lad.tree)[-c(1:(Ntip(lad.tree)+1))])
        avgnodes <- GrpCntTrees[[i]]@data[order(as.numeric(GrpCntTrees[[i]]@data$node)),c('MU','EPSILON')]
        hpd95 <- unlist(GrpCntTrees[[i]]@data[order(as.numeric(GrpCntTrees[[i]]@data$node)),'0.95HPD']) %>% matrix(ncol = 2, byrow = T)
        
        # Plot the tree with custom styling
        plot(lad.tree, edge.color = EdgeBarTxtColors[[1]], show.tip.label = FALSE, x.lim = c(0,750)-MaxBrTimeCompan[[1]][i], main = paste('Centroid Tree of Cluster ', as.character(i), ': ', GrpCntTrees[[i]]@file, sep = ''), cex.main = 3, edge.width = 4)
        # Add geological time scale grid lines
        abline(v=max(branching.times(lad.tree))-seq(0,600,100), lty = 'solid', lend = 'butt', col = "gray30", lwd = 3)
        abline(v=max(branching.times(lad.tree))-seq(50,650,100), lty = 'dotdash', lend = 'butt', col = "gray60", lwd = 3)
        # Add HPD bars
        segments(x0=nodebar_pos[,1], y0=nodebar_pos[,3], x1=nodebar_pos[,2], y1=nodebar_pos[,3], lwd = 4, col = 'royalblue', lend = 'butt')
        # Add mean age points with color coding based on HPD containment
        points(x=max(branching.times(lad.tree)) - avgnodes$MU, y=nodebar_pos[,3], pch = ifelse(avgnodes$MU >= hpd95[,1] & avgnodes$MU <= hpd95[,2], 18, 23), cex = 3, col = ifelse(branching.times(lad.tree)[-1] >= avgnodes$MU, 'darkolivegreen','orange'))
        # Mark calibrated nodes
        nodelabels(node = as.numeric(unlist(strsplit(lad.tree$node.label[grep('@',lad.tree$node.label)], '@'))) + ifelse(is.rooted(lad.tree), 1, 0), frame = 'none', pch = 8, adj = c(1.5,1), cex = 3)
        # Add taxonomic color legend
        segments(x0 = max(branching.times(lad.tree)) + 7, y0 = EdgeBarTxtColors[[2]]$by0, x1 = max(branching.times(lad.tree)) + 7, y1 = EdgeBarTxtColors[[2]]$by1, lwd = 8, col = EdgeBarTxtColors[[2]]$btc, lend = 'butt')
        if (i == 1){
            text(x = max(branching.times(lad.tree)) + 14, y = EdgeBarTxtColors[[2]]$ty, labels = EdgeBarTxtColors[[4]], cex = 4, col = EdgeBarTxtColors[[2]]$btc, adj = c(0,0.5))
            # Add geological period color bars
            rect(max(branching.times(lad.tree))-c(66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8, 650), 185, max(branching.times(lad.tree))-c(0, 66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8), 191, border = NA, col = c('#8ecfc9','#ffbe7a','#fa7f6f','#82b0d2cc','#beb8dc','#e7dad2','#999900','#c497b2','#1450908a','#a9b8c6','#8e8bfe'))
            text(x=max(branching.times(lad.tree))-c(33, 105.5, 173.2, 226.651, 275.401, 328.9, 389.05, 431.5, 464.6, 512.1, 570), y=188, labels = c('Cz', 'K', 'J', 'T', 'P', 'Cn', 'D', 'S', 'O', 'Cm', 'PreCm'), cex = 2)
        }
        if (i == length(GrpCntTrees)){
            # Add time axis and legend for the last plot
            axisPhylo(cex.axis=3, lwd = 3, padj = 1)
            mtext('Ma', side=1, at = 640, cex = 2, padj = 1.5)
            legend(x = 'bottomleft', legend = c('Cluster averaged median (<= node age, dropped in 95%HPD)', 'Cluster averaged median (<= node age, out of 95%HPD)', 'Cluster averaged median (> node age, dropped in 95%HPD)', 'Cluster averaged median (> node age, out of 95%HPD)', 'Calibrated nodes', EdgeBarTxtColors[[3]][,2]), pch = c(18, 23, 18, 23, 8, rep(15, nrow(EdgeBarTxtColors[[3]]))), col = c('darkolivegreen', 'darkolivegreen', 'orange', 'orange', 'black', EdgeBarTxtColors[[3]][,1]), text.col = c('darkolivegreen', 'darkolivegreen', 'orange', 'orange', 'black', EdgeBarTxtColors[[3]][,1]), cex = 2, bg = 'white', box.lwd =3)
        }
    }
    dev.off()
    
    # Generate SVG plot for averaged trees of all clusters
    svg(filename = paste(fileprefix,'-MarkedGroupAverageTrees.svg',sep=''), width = 40, height = 80, pointsize = 12)
    layout(mat = matrix(c(1:length(GrpAvgTrees)),nrow = length(GrpAvgTrees)))
    for (i in 1:length(GrpAvgTrees)){
        lad.tree <- ladderize(GrpAvgTrees[[i]]@phylo, right = T)
        if (i == 1){
            EdgeBarTxtColors <- EdgeColors(lad.tree, ClassShort)
        }
        # Calculate positions for strict and relaxed HPD bars
        nodebar_pos <- cbind(max(branching.times(lad.tree)) - Reduce(rbind, GrpAvgTrees[[i]]@data$`0.95HPD_RLX`[order(as.numeric(GrpAvgTrees[[i]]@data$node))]), 
                             max(branching.times(lad.tree)) - Reduce(rbind, GrpAvgTrees[[i]]@data$`0.95HPD_STR`[order(as.numeric(GrpAvgTrees[[i]]@data$node))]), 
                             node.height(lad.tree)[-c(1:(Ntip(lad.tree)+1))])
        avgnodes <- GrpAvgTrees[[i]]@data[order(as.numeric(GrpAvgTrees[[i]]@data$node)),'EPSILON']
        
        plot(lad.tree, edge.color = EdgeBarTxtColors[[1]], show.tip.label = FALSE, x.lim = c(0,750)-MaxBrTimeCompan[[2]][i], main = paste('Average Tree of Cluster ', as.character(i), ': ', GrpAvgTrees[[i]]@file, sep = ''), cex.main = 3, edge.width = 4)
        abline(v=max(branching.times(lad.tree))-seq(0,600,100), lty = 'solid', lend = 'butt', col = "gray30", lwd = 3)
        abline(v=max(branching.times(lad.tree))-seq(50,650,100), lty = 'dotdash', lend = 'butt', col = "gray60", lwd = 3)
        # Add relaxed HPD bars (thicker)
        segments(x0=nodebar_pos[,1], y0=nodebar_pos[,5], x1=nodebar_pos[,2], y1=nodebar_pos[,5], lwd = 6, col = 'royalblue', lend = 'butt')
        # Add strict HPD bars (thinner, overlaid)
        segments(x0=nodebar_pos[,3], y0=nodebar_pos[,5], x1=nodebar_pos[,4], y1=nodebar_pos[,5], lwd = 4, col = 'coral', lend = 'butt')
        nodelabels(node = as.numeric(unlist(strsplit(lad.tree$node.label[grep('@',lad.tree$node.label)], '@'))) + ifelse(is.rooted(lad.tree), 1, 0), frame = 'none', pch = 8, adj = c(1.5,1), cex = 3)
        segments(x0 = max(branching.times(lad.tree)) + 7, y0 = EdgeBarTxtColors[[2]]$by0, x1 = max(branching.times(lad.tree)) + 7, y1 = EdgeBarTxtColors[[2]]$by1, lwd = 8, col = EdgeBarTxtColors[[2]]$btc, lend = 'butt')
        if (i == 1){
            text(x = max(branching.times(lad.tree)) + 14, y = EdgeBarTxtColors[[2]]$ty, labels = EdgeBarTxtColors[[4]], cex = 4, col = EdgeBarTxtColors[[2]]$btc, adj = c(0,0.5))
            rect(max(branching.times(lad.tree))-c(66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8, 650), 185, max(branching.times(lad.tree))-c(0, 66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8), 191, border = NA, col = c('#8ecfc9','#ffbe7a','#fa7f6f','#82b0d2cc','#beb8dc','#e7dad2','#999900','#c497b2','#1450908a','#a9b8c6','#8e8bfe'))
            text(x=max(branching.times(lad.tree))-c(33, 105.5, 173.2, 226.651, 275.401, 328.9, 389.05, 431.5, 464.6, 512.1, 570), y=188, labels = c('Cz', 'K', 'J', 'T', 'P', 'Cn', 'D', 'S', 'O', 'Cm', 'PreCm'), cex = 2)
        }
        if (i == length(GrpAvgTrees)){
            axisPhylo(cex.axis=3, lwd = 3, padj = 1)
            mtext('Ma', side=1, at = 640, cex = 2, padj = 1.5)
            legend(x = 'bottomleft', legend = c('Relaxed 95% HPD', 'Strict 95% HPD','Calibrated nodes', EdgeBarTxtColors[[3]][,2]), pch = c(22, 22, 8, rep(15, nrow(EdgeBarTxtColors[[3]]))), col = c('royalblue','coral', 'black', EdgeBarTxtColors[[3]][,1]), text.col = c('royalblue','coral', 'black', EdgeBarTxtColors[[3]][,1]), cex = 2, bg = 'white', box.lwd =3)
        }
    }
    dev.off()
    
    # Generate combined plots showing both centroid and averaged trees for each cluster
    for (i in 1:length(GrpCntTrees)){
        svg(filename = paste(fileprefix,'-MarkedGroupCentroid&AverageTrees-Group-',as.character(i),'.svg', sep = ''), width = 100, height = 120, pointsize = 20)
        layout(mat = matrix(c(1:2), nrow = 1, ncol = 2))
        
        # Plot centroid tree (left panel)
        lad.tree1 <- ladderize(GrpCntTrees[[i]]@phylo, right = T)
        EdgeBarTxtColors <- EdgeColors(lad.tree1, ClassLong)
        nodebar_pos <- max(branching.times(lad.tree1)) - Reduce(rbind, GrpCntTrees[[i]]@data$`0.95HPD`[order(as.numeric(GrpCntTrees[[i]]@data$node))])
        nodebar_pos <- cbind(nodebar_pos, node.height(lad.tree1)[-c(1:(Ntip(lad.tree1)+1))])
        avgnodes <- GrpCntTrees[[i]]@data[order(as.numeric(GrpCntTrees[[i]]@data$node)),c('MU','EPSILON')]
        hpd95 <- unlist(GrpCntTrees[[i]]@data[order(as.numeric(GrpCntTrees[[i]]@data$node)),'0.95HPD']) %>% matrix(ncol = 2, byrow = T)
        plot(lad.tree1, edge.color = EdgeBarTxtColors[[1]], show.tip.label = TRUE, x.lim = c(0,800)-MaxBrTimeCompan[[1]][i], main = paste('Centroid Tree of Cluster ', as.character(i), ': ', GrpCntTrees[[i]]@file, sep = ''), cex.main = 3, edge.width = 4)
        abline(v=max(branching.times(lad.tree1))-seq(0,600,100), lty = 'solid', lend = 'butt', col = "gray30", lwd = 3)
        abline(v=max(branching.times(lad.tree1))-seq(50,650,100), lty = 'dotdash', lend = 'butt', col = "gray60", lwd = 3)
        segments(x0=nodebar_pos[,1], y0=nodebar_pos[,3], x1=nodebar_pos[,2], y1=nodebar_pos[,3], lwd = 14, col = 'royalblue', lend = 'butt')
        points(x=max(branching.times(lad.tree1)) - avgnodes$MU, y=nodebar_pos[,3], pch = ifelse(avgnodes$MU >= hpd95[,1] & avgnodes$MU <= hpd95[,2], 18, 23), cex = 3, col = ifelse(branching.times(lad.tree1)[-1] >= avgnodes$MU, 'darkolivegreen','orange'))
        nodelabels(node = as.numeric(unlist(strsplit(lad.tree1$node.label[grep('@',lad.tree1$node.label)], '@'))) + ifelse(is.rooted(lad.tree1), 1, 0), frame = 'none', pch = 8, adj = c(1.5,1), cex = 2)
        segments(x0 = max(branching.times(lad.tree1)) + 75, y0 = EdgeBarTxtColors[[2]]$by0, x1 = max(branching.times(lad.tree1)) + 75, y1 = EdgeBarTxtColors[[2]]$by1, lwd = 16, col = EdgeBarTxtColors[[2]]$btc, lend = 'butt')
        text(x = max(branching.times(lad.tree1)) + 80, y = EdgeBarTxtColors[[2]]$ty, labels = EdgeBarTxtColors[[4]], cex = 2, col = EdgeBarTxtColors[[2]]$btc, adj = c(0,0.5))
        rect(max(branching.times(lad.tree1))-c(66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8, 650), 186, max(branching.times(lad.tree1))-c(0, 66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8), 190, border = NA, col = c('#8ecfc9','#ffbe7a','#fa7f6f','#82b0d2cc','#beb8dc','#e7dad2','#999900','#c497b2','#1450908a','#a9b8c6','#8e8bfe'))
        text(x=max(branching.times(lad.tree1))-c(33, 105.5, 173.2, 226.651, 275.401, 328.9, 389.05, 431.5, 464.6, 512.1, 570), y=188, labels = c('Cz', 'K', 'J', 'T', 'P', 'Cn', 'D', 'S', 'O', 'Cm', 'PreCm'), cex = 2)
        axisPhylo(cex.axis=3, lwd = 3, padj = 1)
        mtext('Ma', side=1, at = 640, cex = 2, padj = 1.5)
        legend(x = 'bottomleft', legend = c('Cluster averaged median (<= node age, dropped in 95%HPD)', 'Cluster averaged median (<= node age, out of 95%HPD)', 'Cluster averaged median (> node age, dropped in 95%HPD)', 'Cluster averaged median (> node age, out of 95%HPD)', 'Calibrated nodes', EdgeBarTxtColors[[3]][,2]), pch = c(18, 23, 18, 23, 8, rep(15, nrow(EdgeBarTxtColors[[3]]))), col = c('darkolivegreen', 'darkolivegreen', 'orange', 'orange', 'black', EdgeBarTxtColors[[3]][,1]), text.col = c('darkolivegreen', 'darkolivegreen', 'orange', 'orange', 'black', EdgeBarTxtColors[[3]][,1]), cex = 2, bg = 'white', box.lwd =3)
        
        # Plot averaged tree (right panel)
        lad.tree2 <- ladderize(GrpAvgTrees[[i]]@phylo, right = T)
        nodebar_pos <- cbind(max(branching.times(lad.tree2)) - Reduce(rbind, GrpAvgTrees[[i]]@data$`0.95HPD_RLX`[order(as.numeric(GrpAvgTrees[[i]]@data$node))]), 
                             max(branching.times(lad.tree2)) - Reduce(rbind, GrpAvgTrees[[i]]@data$`0.95HPD_STR`[order(as.numeric(GrpAvgTrees[[i]]@data$node))]), 
                             node.height(lad.tree2)[-c(1:(Ntip(lad.tree2)+1))])
        avgnodes <- GrpAvgTrees[[i]]@data[order(as.numeric(GrpAvgTrees[[i]]@data$node)),'EPSILON']
        plot(lad.tree2, edge.color = EdgeBarTxtColors[[1]], show.tip.label = TRUE, x.lim = c(0,800)-MaxBrTimeCompan[[2]][i], main = paste(fileprefix,'-Average Tree of Cluster ', as.character(i), sep = ''), cex.main = 3, edge.width = 4)
        abline(v=max(branching.times(lad.tree2))-seq(0,600,100), lty = 'solid', lend = 'butt', col = "gray30", lwd = 3)
        abline(v=max(branching.times(lad.tree2))-seq(50,650,100), lty = 'dotdash', lend = 'butt', col = "gray60", lwd = 3)
        segments(x0=nodebar_pos[,1], y0=nodebar_pos[,5], x1=nodebar_pos[,2], y1=nodebar_pos[,5], lwd = 14, col = 'royalblue', lend = 'butt')
        segments(x0=nodebar_pos[,3], y0=nodebar_pos[,5], x1=nodebar_pos[,4], y1=nodebar_pos[,5], lwd = 6, col = 'coral', lend = 'butt')
        nodelabels(node = as.numeric(unlist(strsplit(lad.tree2$node.label[grep('@',lad.tree2$node.label)], '@'))) + ifelse(is.rooted(lad.tree2), 1, 0), frame = 'none', pch = 8, adj = c(1.5,1), cex = 2)
        segments(x0 = max(branching.times(lad.tree2)) + 75, y0 = EdgeBarTxtColors[[2]]$by0, x1 = max(branching.times(lad.tree2)) + 75, y1 = EdgeBarTxtColors[[2]]$by1, lwd = 16, col = EdgeBarTxtColors[[2]]$btc, lend = 'butt')
        text(x = max(branching.times(lad.tree2)) + 80, y = EdgeBarTxtColors[[2]]$ty, labels = EdgeBarTxtColors[[4]], cex = 2, col = EdgeBarTxtColors[[2]]$btc, adj = c(0,0.5))
        rect(max(branching.times(lad.tree2))-c(66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8, 650), 186, max(branching.times(lad.tree2))-c(0, 66, 145, 201.4, 251.902, 298.9, 358.9, 419.2, 443.8, 485.4, 538.8), 190, border = NA, col = c('#8ecfc9','#ffbe7a','#fa7f6f','#82b0d2cc','#beb8dc','#e7dad2','#999900','#c497b2','#1450908a','#a9b8c6','#8e8bfe'))
        text(x=max(branching.times(lad.tree2))-c(33, 105.5, 173.2, 226.651, 275.401, 328.9, 389.05, 431.5, 464.6, 512.1, 570), y=188, labels = c('Cz', 'K', 'J', 'T', 'P', 'Cn', 'D', 'S', 'O', 'Cm', 'PreCm'), cex = 2)
        axisPhylo(cex.axis=3, lwd = 3, padj = 1)
        mtext('Ma', side=1, at = 640, cex = 2, padj = 1.5)
        legend(x = 'bottomleft', legend = c('Relaxed 95% HPD', 'Strict 95% HPD','Calibrated nodes', EdgeBarTxtColors[[3]][,2]), pch = c(22, 22, 8, rep(15, nrow(EdgeBarTxtColors[[3]]))), col = c('royalblue','coral', 'black', EdgeBarTxtColors[[3]][,1]), text.col = c('royalblue','coral', 'black', EdgeBarTxtColors[[3]][,1]), cex = 2, bg = 'white', box.lwd =3)
        dev.off()
    }
    return(list(GrpAvgTrees,GrpCntTrees))
}
