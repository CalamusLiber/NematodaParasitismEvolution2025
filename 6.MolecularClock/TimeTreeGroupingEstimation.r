## Best time estimation based on density classification----
# ==============================================================================
# TimeTreeGroupingEstimation.r
# Purpose: Group and classify time-tree estimates from MCMCtree analyses using
#   multiple dimensionality reduction and clustering approaches (PCA-DBSCAN,
#   tSNE-DBSCAN, PRD-DBSCAN, AEC-DBSCAN) to identify robust time estimates.
# ==============================================================================

# Loading libraries----
library(ape)
library(ggtree);library(ggplot2);library(ggalt);library(patchwork)
library(treeio);library(pracma)
library(fpc);library(factoextra);library(dbscan);library(Rtsne)
library(torch);library(brulee)
library(doParallel)
library(igraph)
source('/To/Your/Directory/nematoda/timing/PlotSaveGroupTrees.r')

options(max.print = 10000)

setwd('/To/Your/Directory/nematoda/timing/Res_posttrees')

## PCA of Branching times----
# Perform principal component analysis on branching times of converged trees
# to reduce dimensionality and identify major patterns in time estimates.
pca_res <- prcomp(BrTime[mcmc$cvg=='G', -1], center = T, scale. = T)
pca_res.sum <- summary(pca_res)
# Importance of components:
#                            PC1    PC2     PC3     PC4     PC5     PC6     PC7    PC8     PC9    PC10    PC11    PC12
# Standard deviation     12.1731 4.7262 2.61302 1.19306 1.14120 0.94644 0.79470 0.7161 0.47327 0.39887 0.38000 0.29992
# Proportion of Variance  0.8097 0.1221 0.03731 0.00778 0.00712 0.00489 0.00345 0.0028 0.00122 0.00087 0.00079 0.00049
# Cumulative Proportion   0.8097 0.9318 0.96912 0.97690 0.98402 0.98891 0.99236 0.9952 0.99639 0.99726 0.99805 0.99854

# Using the first 3 PCs
BrTimePC1to3 <- pca_res$x[, 1:3]

# Define a function to find the knee/elbow point in a curve for optimal parameter selection.
# This identifies the point where the curve transitions from steep to gradual slope.
findKnee <- function(Cv, method='maxPerpend'){
    if (method=='maxPerpend'){ # it is actually to maximize the parallelepiped areas instead of the perpendicular lines.
        CvVect1 <- rbind(Cv[,1]-Cv,0)
        CvVect2 <- rbind(Cv[,ncol(Cv)]-Cv,0)
        crsprod <- pracma::cross(CvVect1,CvVect2)
        knee_id <- which.max(apply(crsprod, 2, FUN = function(x){sqrt(sum(x^2))}))
        return(c(knee_id, Cv[2,knee_id]))
    }else if (method=='paraTangent'){
        CvTan <- c(diff(Cv[2,], lag=1)/diff(Cv[1,], lag=1), Inf)
        blTan <- (Cv[2,ncol(Cv)]-Cv[2,1])/(Cv[1,ncol(Cv)]-Cv[1,1])
        knee_id <- which.min(abs(CvTan-blTan))
        return(c(knee_id, Cv[2,knee_id]))
    }
}

# Calculate k-nearest neighbor distances to determine optimal epsilon for DBSCAN.
# minPts = 2 * dim(dataset) = 6 (Sander et al., 1998)
knn <- sort(kNNdist(dist(BrTimePC1to3, 'euclidean'), k=2 * ncol(BrTimePC1to3) - 1))
eps_opt_id <- findKnee(Cv=rbind(point=1:length(knn),knn), method='maxPerpend')[1]
eps_opt_pca <- knn[eps_opt_id]

# Plot k-NN distances to visualize the optimal epsilon selection.
svg(filename = 'k-Nearest Neighbor Distances plot for PCA-DBSCAN.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
kNNdistplot(dist(BrTimePC1to3, 'euclidean'), k=2 * ncol(BrTimePC1to3) - 1, minPts = 2 * ncol(BrTimePC1to3))
lines(cbind(c(1,length(knn)),knn[c(1,length(knn))]), lty = 'dotted', lwd=2, col = 'blue')
abline(h=eps_opt_pca, col='red', lty = 'dashed')
abline(v=eps_opt_id, col='red', lty = 'dashed')
points(x=eps_opt_id, y=eps_opt_pca, col='red', pch=16)
title('k-Nearest Neighbor Distances Plot for PCA-DBSCAN')
dev.off()

# PCA-DBSCAN-clustering using optimal epsilon----
# Perform density-based clustering on PCA-reduced branching times.
pca_db <- dbscan::dbscan(dist(BrTimePC1to3, 'euclidean'), minPts = 2 * ncol(BrTimePC1to3), eps = eps_opt_pca)

# Assign cluster labels, with -1 indicating noise/outliers.
pca_group <- rep(-1, nrow(mcmc))
pca_group[mcmc$cvg=='G'] <- pca_db$cluster
pca_table <- addmargins(table(pca_group, mcmc$cst, mcmc$fp, mcmc$part, mcmc$corr, dnn = c('group', 'cst', 'fp', 'part', 'corr')))

# Visualize clustering results in 2D projections of the PCA space.
svg(filename = 'DBSCAN-clustering using PCA-axes 1-3 (shown on axis 1-2).svg', width = 8, height = 6, bg = 'white', pointsize = 12)
fviz_cluster(pca_db, BrTimePC1to3[,c(1,2)], show.clust.cent = F, stand = F, ellipse = T, ellipse.type = 'norm', labelsize = 2, shape=16, geom = 'point', main = 'DBSCAN-clustering using PCA-axes 1:3 (shown on axis 1 & 2)', palette = scales::hue_pal()(max(pca_group)))
dev.off()
svg(filename = 'DBSCAN-clustering using PCA-axes 1-3 (shown on axis 1-3).svg', width = 8, height = 6, bg = 'white', pointsize = 12)
fviz_cluster(pca_db, BrTimePC1to3[,c(1,3)], show.clust.cent = F, stand = F, ellipse = T, ellipse.type = 'norm', labelsize = 2, shape=16, geom = 'point', main = 'DBSCAN-clustering using PCA-axes 1:3 (shown on axis 1 & 3)', palette = scales::hue_pal()(max(pca_group)))
dev.off()
svg(filename = 'DBSCAN-clustering using PCA-axes 1-3 (shown on axis 2-3).svg', width = 8, height = 6, bg = 'white', pointsize = 12)
fviz_cluster(pca_db, BrTimePC1to3[,c(2,3)], show.clust.cent = F, stand = F, ellipse = T, ellipse.type = 'norm', labelsize = 2, shape=16, geom = 'point', main = 'DBSCAN-clustering using PCA-axes 1:3 (shown on axis 2 & 3)', palette = scales::hue_pal()(max(pca_group)))
dev.off()

# Output frequency tables summarizing clustering results and parameter associations.
sink('FreqTables.PCA-DBSCAN.txt')
cat('====PCA-DBSCAN parameters====\n')
pca_db

cat('\n====PCA grouping vs. Calibration Strategy vs. Calibration Coverage====\n')
pca_table

cat('\n====BirthDeath vs. RgeneRate vs. Convergence====\n')
addmargins(table(mcmc$sf, mcmc$rg, mcmc$cvg, dnn = c('sf', 'rg', 'cvg')))

cat('\n====Calibration Strategy vs. Calibration Coverage vs. Convergence====\n')
addmargins(table(mcmc$cst, mcmc$fp, mcmc$cvg, dnn = c('cst', 'fp', 'cvg')))

cat('\n====Data Partition vs. Clock Model vs. Convergence====\n')
addmargins(table(mcmc$part, mcmc$corr, mcmc$cvg, dnn = c('part', 'corr', 'cvg')))
sink()

# t-SNE of Branching times----
# Apply t-distributed Stochastic Neighbor Embedding to reduce dimensionality
# while preserving local structure in the branching time data.
tsne_res <- Rtsne(BrTime[mcmc$cvg=='G', -1], dims = 3, theta = 0, pca = TRUE, pca_center = TRUE, pca_scale = TRUE)
tsne_res.sum <- summary(tsne_res)
BrTimeTSNE3d <- tsne_res$Y
colnames(BrTimeTSNE3d) <- c('tSNE1', 'tSNE2', 'tSNE3')

# Calculate k-NN distances for tSNE space to find optimal epsilon.
knn <- sort(kNNdist(dist(BrTimeTSNE3d, 'euclidean'), k=2 * ncol(BrTimeTSNE3d) - 1))
eps_opt_id <- findKnee(Cv=rbind(point=1:length(knn),knn), method='maxPerpend')[1]
eps_opt_tsne <- knn[eps_opt_id]

# Plot k-NN distances for tSNE-DBSCAN.
svg(filename = 'k-Nearest Neighbor Distances plot for tSNE-DBSCAN.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
kNNdistplot(dist(BrTimeTSNE3d, 'euclidean'), k=2 * ncol(BrTimeTSNE3d) - 1, minPts = 2 * ncol(BrTimeTSNE3d))
lines(cbind(c(1,length(knn)),knn[c(1,length(knn))]), lty = 'dotted', lwd=2, col = 'blue')
abline(h=eps_opt_tsne, col='red', lty = 'dashed')
abline(v=eps_opt_id, col='red', lty = 'dashed')
points(x=eps_opt_id, y=eps_opt_tsne, col='red', pch=16)
title('k-Nearest Neighbor Distances Plot for tSNE-DBSCAN')
dev.off()

## tSNE-DBSCAN-clustering using optimal epsilon----
# Perform DBSCAN clustering on tSNE-reduced data.
tsne_db <- dbscan::dbscan(dist(BrTimeTSNE3d, 'euclidean'), minPts = 2 * ncol(BrTimeTSNE3d), eps = eps_opt_tsne)

tsne_group <- rep(-1, nrow(mcmc))
tsne_group[mcmc$cvg=='G'] <- tsne_db$cluster
tsne_table <- addmargins(table(tsne_group, mcmc$cst, mcmc$fp, mcmc$part, mcmc$corr, dnn = c('group', 'cst', 'fp', 'part', 'corr')))

# Visualize tSNE clustering in 2D projections.
svg(filename = 'DBSCAN-clustering using tSNE-axes 1-3 (shown on axis 1-2).svg', width = 8, height = 6, bg = 'white', pointsize = 12)
fviz_cluster(tsne_db, BrTimeTSNE3d[,c(1,2)], show.clust.cent = F, stand = F, ellipse = T, ellipse.type = 'norm', labelsize = 2, shape=16, xlim=c(-12, 10), ylim=c(-10, 15), geom = 'point', main = 'DBSCAN-clustering using tSNE-axes 1:3 (shown on axis 1 & 2)', palette = scales::hue_pal()(max(tsne_group)))
dev.off()
svg(filename = 'DBSCAN-clustering using tSNE-axes 1-3 (shown on axis 1-3).svg', width = 8, height = 6, bg = 'white', pointsize = 12)
fviz_cluster(tsne_db, BrTimeTSNE3d[,c(1,3)], show.clust.cent = F, stand = F, ellipse = T, ellipse.type = 'norm', labelsize = 2, shape=16, xlim=c(-12, 10), ylim=c(-25, 15), geom = 'point', main = 'DBSCAN-clustering using tSNE-axes 1:3 (shown on axis 1 & 3)', palette = scales::hue_pal()(max(tsne_group)))
dev.off()
svg(filename = 'DBSCAN-clustering using tSNE-axes 1-3 (shown on axis 2-3).svg', width = 8, height = 6, bg = 'white', pointsize = 12)
fviz_cluster(tsne_db, BrTimeTSNE3d[,c(2,3)], show.clust.cent = F, stand = F, ellipse = T, ellipse.type = 'norm', labelsize = 2, shape=16, xlim=c(-10, 15), ylim=c(-25, 15), geom = 'point', main = 'DBSCAN-clustering using tSNE-axes 1:3 (shown on axis 2 & 3)', palette = scales::hue_pal()(max(tsne_group)))
dev.off()

# Output frequency tables for tSNE clustering.
sink('FreqTables.tSNE-DBSCAN.txt')
cat('====tSNE-DBSCAN parameters====\n')
tsne_db

cat('\n====tSNE grouping vs. Calibration Strategy vs. Calibration Coverage====\n')
tsne_table

cat('\n====BirthDeath vs. RgeneRate vs. Convergence====\n')
addmargins(table(mcmc$sf, mcmc$rg, mcmc$cvg, dnn = c('sf', 'rg', 'cvg')))

cat('\n====Calibration Strategy vs. Calibration Coverage vs. Convergence====\n')
addmargins(table(mcmc$cst, mcmc$fp, mcmc$cvg, dnn = c('cst', 'fp', 'cvg')))

cat('\n====Data Partition vs. Clock Model vs. Convergence====\n')
addmargins(table(mcmc$part, mcmc$corr, mcmc$cvg, dnn = c('part', 'corr', 'cvg')))

cat('\n====The End====\n')
sink()

## DBSCAN-clustering using customized range distance matrix between paired matrices (PRD)----
# Define a custom distance function for comparing time range matrices.
# This calculates distances based on HPD intervals and node weights.
GoodTimeRange3 <- TimeRange3[mcmc$cvg=='G']

RangeDist <- function(mat1, mat2, w=1, dtype='t', byrow=TRUE){
    if (! byrow){
        mat1 <- t(mat1)
        mat2 <- t(mat2)
    }
    mat1 <- apply(mat1, 2, sort)
    mat2 <- apply(mat2, 2, sort)
    if (length(mat1)==3){
        mat1[3,] <- apply(mat1, 2, function(x){x[3]*diff(x)[1]/diff(x)[2]})
        mat2[3,] <- apply(mat2, 2, function(x){x[3]*diff(x)[1]/diff(x)[2]})
        mat1 <- mat1[-2,]
        mat2 <- mat2[-2,]
    }
    if (length(w)==ncol(mat1)){
        w <- w
    }else if (length(w)==1){
        w <- rep(w, ncol(mat1))
    }else{
        break
    }
    if (dtype=='t'){
        d <- abs(mat1[2,] - mat2[2,]) + abs(mat1[1,] - mat2[1,])
    }else if (dtype=='r'){
        d <- abs(mat1[2,] - mat2[1,]) / abs(mat1[1,] - mat2[2,])
        d[d<1] <- 1 / d[d<1]
    }else{
        break
    }
    return(sum(d * w)/sum(w))
}

# Calculate node weights based on subtree sizes for weighted distance computation.
nodewight <- unlist(lapply(subtrees(AllTrees[[1]]@phylo),Nnode))[-1]/unlist(lapply(subtrees(AllTrees[[1]]@phylo),Nnode))[1]

# Compute pairwise distance matrix using parallel processing.
TreeDistMatrix <- matrix(NA, nrow = length(TimeRange3[mcmc$cvg=='G']), ncol = length(TimeRange3[mcmc$cvg=='G']))

cl <- makeCluster(0.5 * detectCores())  # 0.5 * detectCores() or user-specified
registerDoParallel(cl)

TreeDistMatrix <- foreach (i = 1:length(GoodTimeRange3), .combine = rbind) %dopar% {
    Reduce(c, lapply(GoodTimeRange3, function(x){RangeDist(GoodTimeRange3[[i]], x, w=nodewight, dtype = 't')}))
}

stopCluster(cl)

# Find optimal epsilon for PRD distance matrix.
knn <- sort(kNNdist(TreeDistMatrix, k=2 * 3 - 1))
eps_opt_id <- findKnee(Cv=rbind(point=1:length(knn),knn), method='maxPerpend')[1]
eps_opt_prd <- knn[eps_opt_id]

# Plot k-NN distances for PRD-DBSCAN.
svg(filename = 'k-Nearest Neighbor Distances plot for PRD-DBSCAN.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
kNNdistplot(TreeDistMatrix, k=2 * 3 - 1, minPts = 2 * 3)
lines(cbind(c(1,length(knn)),knn[c(1,length(knn))]), lty = 'dotted', lwd=2, col = 'blue')
abline(h=eps_opt_prd, col='red', lty = 'dashed')
abline(v=eps_opt_id, col='red', lty = 'dashed')
points(x=eps_opt_id, y=eps_opt_prd, col='red', pch=16)
title('k-Nearest Neighbor Distances Plot for PRD-DBSCAN')
dev.off()

## PRD-DBSCAN-clustering using optimal epsilon----
# Perform DBSCAN on the custom range distance matrix.
prd_db <- dbscan::dbscan(TreeDistMatrix, minPts = 2 * 3, eps = eps_opt_prd)

prd_group <- rep(-1, nrow(mcmc))
prd_group[mcmc$cvg=='G'] <- prd_db$cluster
prd_table <- addmargins(table(prd_group, mcmc$cst, mcmc$fp, mcmc$part, mcmc$corr, dnn = c('group', 'cst', 'fp', 'part', 'corr')))

# Visualize PRD clustering in tSNE projections.
svg(filename = 'PRD-DBSCAN-clustering using tSNE-axes 1-3 (shown on axis 1-2).svg', width = 8, height = 6, bg = 'white', pointsize = 12)
fviz_cluster(prd_db, BrTimeTSNE3d[,c(1,2)], show.clust.cent = F, stand = F, ellipse = T, ellipse.type = 'norm', labelsize = 2, shape=16, xlim=c(-12, 10), ylim=c(-10, 15), geom = 'point', main = 'PRD-DBSCAN-clustering using tSNE-axes 1:3 (shown on axis 1 & 2)', palette = scales::hue_pal()(max(prd_group)))
dev.off()
svg(filename = 'PRD-DBSCAN-clustering using tSNE-axes 1-3 (shown on axis 1-3).svg', width = 8, height = 6, bg = 'white', pointsize = 12)
fviz_cluster(prd_db, BrTimeTSNE3d[,c(1,3)], show.clust.cent = F, stand = F, ellipse = T, ellipse.type = 'norm', labelsize = 2, shape=16, xlim=c(-12, 10), ylim=c(-25, 15), geom = 'point', main = 'PRD-DBSCAN-clustering using tSNE-axes 1:3 (shown on axis 1 & 3)', palette = scales::hue_pal()(max(prd_group)))
dev.off()
svg(filename = 'PRD-DBSCAN-clustering using tSNE-axes 1-3 (shown on axis 2-3).svg', width = 8, height = 6, bg = 'white', pointsize = 12)
fviz_cluster(prd_db, BrTimeTSNE3d[,c(2,3)], show.clust.cent = F, stand = F, ellipse = T, ellipse.type = 'norm', labelsize = 2, shape=16, xlim=c(-10, 15), ylim=c(-25, 15), geom = 'point', main = 'PRD-DBSCAN-clustering using tSNE-axes 1:3 (shown on axis 2 & 3)', palette = scales::hue_pal()(max(prd_group)))
dev.off()

# Output frequency tables for PRD clustering.
sink('FreqTables.PRD-DBSCAN.txt')
cat('====PRD-DBSCAN parameters====\n')
prd_db

cat('\n====PRD grouping vs. Calibration Strategy vs. Calibration Coverage====\n')
prd_table

cat('\n====BirthDeath vs. RgeneRate vs. Convergence====\n')
addmargins(table(mcmc$sf, mcmc$rg, mcmc$cvg, dnn = c('sf', 'rg', 'cvg')))

cat('\n====Calibration Strategy vs. Calibration Coverage vs. Convergence====\n')
addmargins(table(mcmc$cst, mcmc$fp, mcmc$cvg, dnn = c('cst', 'fp', 'cvg')))

cat('\n====Data Partition vs. Clock Model vs. Convergence====\n')
addmargins(table(mcmc$part, mcmc$corr, mcmc$cvg, dnn = c('part', 'corr', 'cvg')))

cat('\n====The End====\n')
sink()

## DBSCAN-clustering using AutoEncoder clustering (AEC) with 95HPD distance matrix----
# Define an autoencoder function to learn compressed representations of time range matrices.
MatEncoder <- function(mat, outdim=1, epochs=50){
    torch_manual_seed('20250410')
    data <- torch_tensor(mat, dtype = torch_float())
    autoencoder <- nn_module(
        initialize = function() {
            self$encoder <- nn_sequential(
                nn_linear(ncol(mat), outdim),  # compressed to 1d
                nn_relu()
            )
            self$decoder <- nn_sequential(
                nn_linear(outdim, ncol(mat)),  # recovered to original dimensionality
                nn_sigmoid()
            )
        },
        forward = function(x) {
            encoded <- self$encoder(x)  # encode
            decoded <- self$decoder(encoded)  # decode
            list(encoded, decoded)  
        }
    )
    model <- autoencoder()
    # loss function and optimizer
    criterion <- nn_mse_loss()  # MSE loss
    optimizer <- optim_adam(model$parameters, lr = 0.01)
    # model training
    batch_size <- 1
    num_batches <- nrow(data) / batch_size
    cat('\n==== AutoEncoder Training ====\n')
    for (epoch in 1:epochs) {
        model$train()
        total_loss <- 0
        for (i in seq(1, nrow(data), by = batch_size)) {
            batch <- data[i:(i + batch_size - 1), ]
            optimizer$zero_grad()
            outputs <- model(batch)
            encoded <- outputs[[1]]  # extracting potential features
            decoded <- outputs[[2]]
            loss <- criterion(decoded, batch)
            loss$backward()
            optimizer$step()
            total_loss <- total_loss + loss$item()
        }
        cat(sprintf("Epoch %d, Loss: %.4f\n", epoch, total_loss / num_batches))
    }
    # extracting potential features by decoder
    return(as.numeric(model$encoder(data)))
}

# Train autoencoders on each time range matrix to extract latent features.
sink('AutoEncoder-training.log')
aec_res <- Reduce(rbind, lapply(GoodTimeRange3, MatEncoder, epochs=500))
sink()

colnames(aec_res) <- c('AEC1', 'AEC2', 'AEC3')

# Find optimal epsilon for AEC features.
knn <- sort(kNNdist(dist(aec_res, 'euclidean'), k=2 * ncol(aec_res) - 1))
eps_opt_id <- findKnee(Cv=rbind(point=1:length(knn),knn), method='maxPerpend')[1]
eps_opt_aec <- knn[eps_opt_id]

# Plot k-NN distances for AEC-DBSCAN.
svg(filename = 'k-Nearest Neighbor Distances plot for AEC-DBSCAN.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
kNNdistplot(dist(aec_res, 'euclidean'), k=2 * ncol(aec_res) - 1, minPts = 2 * ncol(aec_res))
lines(cbind(c(1,length(knn)),knn[c(1,length(knn))]), lty = 'dotted', lwd=2, col = 'blue')
abline(h=eps_opt_aec, col='red', lty = 'dashed')
abline(v=eps_opt_id, col='red', lty = 'dashed')
points(x=eps_opt_id, y=eps_opt_aec, col='red', pch=16)
title('k-Nearest Neighbor Distances Plot for AEC-DBSCAN')
dev.off()

## AEC-DBSCAN-clustering using optimal epsilon----
# Perform DBSCAN on autoencoder features.
aec_db <- dbscan::dbscan(dist(aec_res, 'euclidean'), minPts = 2 * ncol(aec_res), eps = eps_opt_aec)

aec_group <- rep(-1, nrow(mcmc))
aec_group[mcmc$cvg=='G'] <- aec_db$cluster
aec_table <- addmargins(table(aec_group, mcmc$cst, mcmc$fp, mcmc$part, mcmc$corr, dnn = c('group', 'cst', 'fp', 'part', 'corr')))

svg(filename = 'AEC-DBSCAN-clustering using tSNE-axes 1-3 (shown on axis 1-2).svg', width = 8, height = 6, bg = 'white', pointsize = 12)
fviz_cluster(aec_db, BrTimeTSNE3d[,c(1,2)], show.clust.cent = F, stand = F, ellipse = T, ellipse.type = 'norm', labelsize = 2, shape=16, xlim=c(-12, 10), ylim=c(-10, 15), geom = 'point', main = 'AEC-DBSCAN-clustering using tSNE-axes 1:3 (shown on axis 1 & 2)', palette = scales::hue_pal()(max(aec_group)))
dev.off()
svg(filename = 'AEC-DBSCAN-clustering using tSNE-axes 1-3 (shown on axis 1-3).svg', width = 8, height = 6, bg = 'white', pointsize = 12)
fviz_cluster(aec_db, BrTimeTSNE3d[,c(1,3)], show.clust.cent = F, stand = F, ellipse = T, ellipse.type = 'norm', labelsize = 2, shape=16, xlim=c(-12, 10), ylim=c(-25, 15), geom = 'point', main = 'AEC-DBSCAN-clustering using tSNE-axes 1:3 (shown on axis 1 & 3)', palette = scales::hue_pal()(max(aec_group)))
dev.off()
svg(filename = 'AEC-DBSCAN-clustering using tSNE-axes 1-3 (shown on axis 2-3).svg', width = 8, height = 6, bg = 'white', pointsize = 12)
fviz_cluster(aec_db, BrTimeTSNE3d[,c(2,3)], show.clust.cent = F, stand = F, ellipse = T, ellipse.type = 'norm', labelsize = 2, shape=16, xlim=c(-10, 15), ylim=c(-25, 15), geom = 'point', main = 'AEC-DBSCAN-clustering using tSNE-axes 1:3 (shown on axis 2 & 3)', palette = scales::hue_pal()(max(aec_group)))
dev.off()

sink('FreqTables.AEC-DBSCAN.txt')
cat('====AEC-DBSCAN parameters====\n')
aec_db

cat('\n====AEC grouping vs. Calibration Strategy vs. Calibration Coverage====\n')
aec_table

cat('\n====BirthDeath vs. RgeneRate vs. Convergence====\n')
addmargins(table(mcmc$sf, mcmc$rg, mcmc$cvg, dnn = c('sf', 'rg', 'cvg')))

cat('\n====Calibration Strategy vs. Calibration Coverage vs. Convergence====\n')
addmargins(table(mcmc$cst, mcmc$fp, mcmc$cvg, dnn = c('cst', 'fp', 'cvg')))

cat('\n====Data Partition vs. Clock Model vs. Convergence====\n')
addmargins(table(mcmc$part, mcmc$corr, mcmc$cvg, dnn = c('part', 'corr', 'cvg')))

cat('\n====The End====\n')
sink()

## find the center point of the highest density of points of all groups----
pca_centpoints.ind <- as.vector(sapply(1:max(pca_db$cluster), 
                                       FUN=function(x){
                                           which.min(apply(as.matrix(dist(BrTimePC1to3[pca_db$cluster==x,], diag = T, upper = T)), MARGIN = 1, sum))
                                       }, 
                                       simplify = T, USE.NAMES = F))
pca_centpoints.tre <- apply(cbind(1:max(pca_db$cluster),pca_centpoints.ind), 1, 
                            FUN=function(x){
                                mcmc[pca_group==x[1],][x[2],]
                            })
pca_centpoints.tre <- Reduce(rbind,pca_centpoints.tre)

tsne_centpoints.ind <- as.vector(sapply(1:max(tsne_db$cluster), 
                                        FUN=function(x){
                                            which.min(apply(as.matrix(dist(BrTimeTSNE3d[tsne_db$cluster==x,], diag = T, upper = T)), MARGIN = 1, sum))
                                        }, 
                                        simplify = T, USE.NAMES = F))
tsne_centpoints.tre <- apply(cbind(1:max(tsne_db$cluster),tsne_centpoints.ind), 1, 
                             FUN=function(x){
                                 mcmc[tsne_group==x[1],][x[2],]
                             })
tsne_centpoints.tre <- Reduce(rbind,tsne_centpoints.tre)

prd_centpoints.ind <- as.vector(sapply(1:max(prd_db$cluster), 
                                       FUN=function(x){
                                           which.min(apply(TreeDistMatrix[prd_db$cluster==x,prd_db$cluster==x], MARGIN = 1, sum))
                                       }, 
                                       simplify = T, USE.NAMES = F))
prd_centpoints.tre <- apply(cbind(1:max(prd_db$cluster),prd_centpoints.ind), 1, 
                            FUN=function(x){
                                mcmc[prd_group==x[1],][x[2],]
                            })
prd_centpoints.tre <- Reduce(rbind,prd_centpoints.tre)

aec_centpoints.ind <- as.vector(sapply(1:max(aec_db$cluster), 
                                       FUN=function(x){
                                           which.min(apply(as.matrix(dist(aec_res[aec_db$cluster==x,], diag = T, upper = T)), MARGIN = 1, sum))
                                       }, 
                                       simplify = T, USE.NAMES = F))
aec_centpoints.tre <- apply(cbind(1:max(aec_db$cluster),aec_centpoints.ind), 1, 
                            FUN=function(x){
                                mcmc[aec_group==x[1],][x[2],]
                            })
aec_centpoints.tre <- Reduce(rbind,aec_centpoints.tre)

mcmc.class.report <- mcmc
mcmc.class.report$pca.group <- pca_group
mcmc.class.report$pca.group[mcmc.class.report$pca.group == -1] <- NA
mcmc.class.report$pca.centroid <- rep(NA, nrow(mcmc.class.report))
mcmc.class.report$pca.centroid[as.numeric(row.names(pca_centpoints.tre))] <- '*'
mcmc.class.report$tsne.group <- tsne_group
mcmc.class.report$tsne.group[mcmc.class.report$tsne.group == -1] <- NA
mcmc.class.report$tsne.centroid <- rep(NA, nrow(mcmc.class.report))
mcmc.class.report$tsne.centroid[as.numeric(row.names(tsne_centpoints.tre))] <- '*'
mcmc.class.report$prd.group <- prd_group
mcmc.class.report$prd.group[mcmc.class.report$prd.group == -1] <- NA
mcmc.class.report$prd.centroid <- rep(NA, nrow(mcmc.class.report))
mcmc.class.report$prd.centroid[as.numeric(row.names(prd_centpoints.tre))] <- '*'
mcmc.class.report$aec.group <- aec_group
mcmc.class.report$aec.group[mcmc.class.report$aec.group == -1] <- NA
mcmc.class.report$aec.centroid <- rep(NA, nrow(mcmc.class.report))
mcmc.class.report$aec.centroid[as.numeric(row.names(aec_centpoints.tre))] <- '*'
write.csv(mcmc.class.report, file = 'mcmc.filtering.and.classification.report.csv', na = '')
rm(knn, eps_opt_id)

## compute the mean and SE of median ages and HPD95% trapezoids for PCA-DBSCAN trees (good tree only)----
xxx <- PlotSaveGroupTrees(good_trees=AllTrees[mcmc$cvg=='G'], cluster=pca_db$cluster, centpoints=pca_centpoints.ind, fileprefix='PCA-DBSCAN')

GrpAvgTrees_PCA <- xxx[[1]]
GrpCntTrees_PCA <- xxx[[2]]

## compute the mean and SE of median ages and HPD95% trapezoids for tSNE-DBSCAN trees (good tree only)----
xxx <- PlotSaveGroupTrees(good_trees=AllTrees[mcmc$cvg=='G'], cluster=tsne_db$cluster, centpoints=tsne_centpoints.ind, fileprefix='tSNE-DBSCAN')

GrpAvgTrees_tSNE <- xxx[[1]]
GrpCntTrees_tSNE <- xxx[[2]]

## compute the mean and SE of median ages and HPD95% trapezoids for PRD-DBSCAN trees (good tree only)----
xxx <- PlotSaveGroupTrees(good_trees=AllTrees[mcmc$cvg=='G'], cluster=prd_db$cluster, centpoints=prd_centpoints.ind, fileprefix='PRD-DBSCAN')

GrpAvgTrees_PRD <- xxx[[1]]
GrpCntTrees_PRD <- xxx[[2]]

## compute the mean and SE of median ages and HPD95% trapezoids for AEC-DBSCAN trees (good tree only)----
xxx <- PlotSaveGroupTrees(good_trees=AllTrees[mcmc$cvg=='G'], cluster=aec_db$cluster, centpoints=aec_centpoints.ind, fileprefix='AEC-DBSCAN')

GrpAvgTrees_AEC <- xxx[[1]]
GrpCntTrees_AEC <- xxx[[2]]

rm(xxx)

# Save group-level summary trees for later inspection and downstream analysis.
saveRDS(GrpCntTrees_PCA, './GrpCntTrees_PCA.rds')
saveRDS(GrpAvgTrees_PCA, './GrpAvgTrees_PCA.rds')
saveRDS(GrpCntTrees_tSNE, './GrpCntTrees_tSNE.rds')
saveRDS(GrpAvgTrees_tSNE, './GrpAvgTrees_tSNE.rds')
saveRDS(GrpCntTrees_PRD, './GrpCntTrees_PRD.rds')
saveRDS(GrpAvgTrees_PRD, './GrpAvgTrees_PRD.rds')
saveRDS(GrpCntTrees_AEC, './GrpCntTrees_AEC.rds')
saveRDS(GrpAvgTrees_AEC, './GrpAvgTrees_AEC.rds')

# divergence time of key nodes----
# Define a fixed set of focal nodes/clades for reporting divergence times.
keyclades <- c('root', 'EuMetazoa', 'Ecdysozoa', 'Arthropoda', 'Tardigrada', 'Nematomorpha', 'Nematoda', 'Enoplia (Clade II)', 'Dorylaimia (Clade I)', 'Chromadoria', 'Rhabditida','Spirurina (Clade III)', 'Tylenchina (Clade IV)', 'Rhabditina (Clade V)', 'Trichinellida', 'Tylenchomorpha', 'Diplogastromorpha', 'Strongyloididae', 'Ascarididae', 'Onchocercidae', 'Strongyloidea')

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
        # Return formatted age strings with 95% HPD ranges for reporting.
        if (tre=='cnt'){
            hdi <- unlist(Map(function(x){paste(round(unlist(phy@data[phy@data$node==x, '0.95HPD']), digit=2), collapse=',')}, keynodes))
        }else if (tre=='avg'){
            hdi <- unlist(Map(function(x){paste(round(unlist(phy@data[phy@data$node==x, '0.95HPD_RLX']), digit=2), collapse=',')}, keynodes))
        }
        bt <- round(branching.times(phy@phylo)[keynodes - Ntip(phy@phylo)],2)
        
        res <- paste(bt, '[', hdi, ']', sep='')
        names(res) <- keyclades
    }else if (out == 'num'){
        # Return numeric age values and HPD bounds for downstream analysis.
        if (tre=='cnt'){
            res <- cbind(bt = branching.times(phy@phylo)[keynodes - Ntip(phy@phylo)], t(sapply(Map(function(x){unlist(phy@data[phy@data$node==x, '0.95HPD'])}, keynodes), c)))
        }else if (tre=='avg'){
            res <- cbind(bt = branching.times(phy@phylo)[keynodes - Ntip(phy@phylo)], t(sapply(Map(function(x){unlist(phy@data[phy@data$node==x, '0.95HPD_RLX'])}, keynodes), c)))
        }
        rownames(res) <- keyclades
    }
    return(res)
}

write.csv(cbind(Cluster=c(1:length(GrpCntTrees_PCA)), Reduce(rbind, lapply(GrpCntTrees_PCA, GetTime, tre='cnt')), Treefile=unlist(Map(function(x){GrpCntTrees_PCA[[x]]@file}, 1:length(GrpCntTrees_PCA)))), file = './KeyCladeTime_GrpCntTrees_PCA.csv', row.names = F)
write.csv(cbind(Cluster=c(1:length(GrpAvgTrees_PCA)), Reduce(rbind, lapply(GrpAvgTrees_PCA, GetTime, tre='avg')), Treefile=unlist(Map(function(x){GrpAvgTrees_PCA[[x]]@file}, 1:length(GrpAvgTrees_PCA)))), file = './KeyCladeTime_GrpAvgTrees_PCA.csv', row.names = F)
write.csv(cbind(Cluster=c(1:length(GrpCntTrees_tSNE)), Reduce(rbind, lapply(GrpCntTrees_tSNE, GetTime, tre='cnt')), Treefile=unlist(Map(function(x){GrpCntTrees_tSNE[[x]]@file}, 1:length(GrpCntTrees_tSNE)))), file = './KeyCladeTime_GrpCntTrees_tSNE.csv', row.names = F)
write.csv(cbind(Cluster=c(1:length(GrpAvgTrees_tSNE)), Reduce(rbind, lapply(GrpAvgTrees_tSNE, GetTime, tre='avg')), Treefile=unlist(Map(function(x){GrpAvgTrees_tSNE[[x]]@file}, 1:length(GrpAvgTrees_tSNE)))), file = './KeyCladeTime_GrpAvgTrees_tSNE.csv', row.names = F)
write.csv(cbind(Cluster=c(1:length(GrpCntTrees_PRD)), Reduce(rbind, lapply(GrpCntTrees_PRD, GetTime, tre='cnt')), Treefile=unlist(Map(function(x){GrpCntTrees_PRD[[x]]@file}, 1:length(GrpCntTrees_PRD)))), file = './KeyCladeTime_GrpCntTrees_PRD.csv', row.names = F)
write.csv(cbind(Cluster=c(1:length(GrpAvgTrees_PRD)), Reduce(rbind, lapply(GrpAvgTrees_PRD, GetTime, tre='avg')), Treefile=unlist(Map(function(x){GrpAvgTrees_PRD[[x]]@file}, 1:length(GrpAvgTrees_PRD)))), file = './KeyCladeTime_GrpAvgTrees_PRD.csv', row.names = F)
write.csv(cbind(Cluster=c(1:length(GrpCntTrees_AEC)), Reduce(rbind, lapply(GrpCntTrees_AEC, GetTime, tre='cnt')), Treefile=unlist(Map(function(x){GrpCntTrees_AEC[[x]]@file}, 1:length(GrpCntTrees_AEC)))), file = './KeyCladeTime_GrpCntTrees_AEC.csv', row.names = F)
write.csv(cbind(Cluster=c(1:length(GrpAvgTrees_AEC)), Reduce(rbind, lapply(GrpAvgTrees_AEC, GetTime, tre='avg')), Treefile=unlist(Map(function(x){GrpAvgTrees_AEC[[x]]@file}, 1:length(GrpAvgTrees_AEC)))), file = './KeyCladeTime_GrpAvgTrees_AEC.csv', row.names = F)

# Remove redundant or empty annotations inserted by treeio::write.nexus() on tip labels.
# This is a workaround for known treeio NEXUS export behavior.
system("treefiles=$(ls ./*.tre); for tr in $treefiles; do sed -i 's#\\[&0\\.95HPD_STR={},0\\.95HPD_RLX={},0\\.95ETI_STR={},0\\.95ETI_RLX={}\\]##g;s#\\[&0\\.95ETI={},0\\.95HPD={},MAD_RANGE={},SD_RANGE={},0\\.95HPD_STR={},0\\.95HPD_RLX={},0\\.95ETI_STR={},0\\.95ETI_RLX={}\\]##g' $tr; done")

# or the following equivalent but runnable only for R 4.0+:
# system(r"(treefiles=$(ls ./*.tre); for tr in $treefiles; do sed -i 's#\[&0\.95HPD_STR={},0\.95HPD_RLX={},0\.95ETI_STR={},0\.95ETI_RLX={}\]##g;s#\[&0\.95ETI={},0\.95HPD={},MAD_RANGE={},SD_RANGE={},0\.95HPD_STR={},0\.95HPD_RLX={},0\.95ETI_STR={},0\.95ETI_RLX={}\]##g' $tr; done)")

# or you should run:
# treefiles=$(ls ./*.tre); for tr in $treefiles; do sed -i 's#\[&0\.95HPD_STR={},0\.95HPD_RLX={},0\.95ETI_STR={},0\.95ETI_RLX={}\]##g;s#\[&0\.95ETI={},0\.95HPD={},MAD_RANGE={},SD_RANGE={},0\.95HPD_STR={},0\.95HPD_RLX={},0\.95ETI_STR={},0\.95ETI_RLX={}\]##g' $tr; done
# in your shell prompt

# generate a PDF copy for each SVG file
system('bash ./svg2pdf.sh r')

# Save the entire R workspace so downstream analyses can reload the full session.
save.image('/To/Your/Directory/nematoda/timing/AllPosteriorTimeTree.cLGPMSF.RData')
