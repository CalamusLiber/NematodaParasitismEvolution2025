# ===============================================================================
# VariableAssessmentRDA.r
# Purpose: Evaluate which MCMC specification and calibration variables explain
#   time-tree variation using redundancy analysis (RDA) across full and converged
#   data subsets, with model selection, permutation testing, and ordination plots.
# ===============================================================================

## Loading libraries----
# vegan : ecological ordination and redundancy analysis functions
# ggplot2: plotting results and biplots
library(vegan)
library(ggplot2)
# library(CCP)

# Set the working directory for output SVG files and R data saving.
setwd('/To/Your/Directory/nematoda/timing/Res_posttrees')

# Perform detrended correspondence analysis (DCA) on branching times.
# This is used to assess whether the data are better suited for RDA or CCA.
# Rule of thumb: if the first four axis lengths are all less than 3, RDA is preferred.
decorana(veg=BrTime)
# Detrended correspondence analysis with 26 segments.
# Rescaling of axes with 4 iterations.
# Total inertia (scaled Chi-square): 0.0368 
# 
#                         DCA1     DCA2     DCA3     DCA4
# Eigenvalues          0.02076 0.005420 0.004294 0.004967
# Additive Eigenvalues 0.02076 0.005346 0.004100 0.001459
# Decorana values      0.02139 0.003622 0.001671 0.000752
# Axis lengths         0.40583 0.331044 0.345910 0.432101
#  Use RDA, if the first 4 axis lengths all less than 3, or otherwise, use CCA.

# Determine available CPU cores for parallel permutation tests.
ncrs <- parallel::detectCores(all.tests = FALSE, logical = TRUE)

## redundancy analysis (full data)----
rdaFL_res <- rda(BrTime ~ ., data=mcmc.dump, scale=TRUE)
rdaFL_res.scaling1 <- summary(rdaFL_res, scaling = 1)

# Ordination diagram (Distance biplot)----
rdaFL_p1 <- ggplot() + 
    geom_line(aes(x=c(1:length(rdaFL_res$CCA$eig)), y=as.vector(rdaFL_res$CCA$eig)), linetype="dashed", linewidth=1.5, color="darkgrey") + 
    geom_point(aes(x=c(1:length(rdaFL_res$CCA$eig)), y=as.vector(rdaFL_res$CCA$eig)), size=3, color="magenta") +
    scale_x_discrete(name = "Ordination axes") + 
    ylab("Inertia") + 
    ggtitle('RDA loading variation (full data)') + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaFL_p1_loading.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaFL_p1
dev.off()

rdaFL_p2.12 <- ggplot() +
    geom_point(aes(x=rdaFL_res$CCA$wa[mcmc.dump$cvg=='B',1], y=rdaFL_res$CCA$wa[mcmc.dump$cvg=='B',2]), col = "gray86") +
    geom_point(aes(x=rdaFL_res$CCA$wa[mcmc.dump$cvg=='G',1], y=rdaFL_res$CCA$wa[mcmc.dump$cvg=='G',2]), col = c('black', scales::hue_pal()(max(prd_group)))[prd_group[mcmc.dump$cvg=='G']+1]) +
    geom_segment(aes(xend=rdaFL_res$CCA$biplot[,1]*0.15, yend=rdaFL_res$CCA$biplot[,2]*0.15, x=0, y=0), colour="black", linewidth=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(aes(x=0.17*rdaFL_res$CCA$biplot[,1], y=0.17*rdaFL_res$CCA$biplot[,2], label = rownames(rdaFL_res$CCA$biplot))) +
    labs(x=paste('RDA 1:',round(rdaFL_res.scaling1$concont$importance[2,1] * 100, 2), '%'), y=paste('RDA 2:',round(rdaFL_res.scaling1$concont$importance[2,2] * 100, 2), '%')) +
    ggtitle('RDA (full data): time ~ mcmc specification (RDA 1 & 2)') + 
    theme_bw() +
    theme(legend.position='right', plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaFL_p2.12.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaFL_p2.12
dev.off()

rdaFL_p2.13 <- ggplot() +
    geom_point(aes(x=rdaFL_res$CCA$wa[mcmc.dump$cvg=='B',1], y=rdaFL_res$CCA$wa[mcmc.dump$cvg=='B',3]), col = "gray86") +
    geom_point(aes(x=rdaFL_res$CCA$wa[mcmc.dump$cvg=='G',1], y=rdaFL_res$CCA$wa[mcmc.dump$cvg=='G',3]), col = c('black', scales::hue_pal()(max(prd_group)))[prd_group[mcmc.dump$cvg=='G']+1]) +
    geom_segment(aes(xend=rdaFL_res$CCA$biplot[,1]*0.15, yend=rdaFL_res$CCA$biplot[,3]*0.15, x=0, y=0), colour="black", linewidth=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(aes(x=0.17*rdaFL_res$CCA$biplot[,1], y=0.17*rdaFL_res$CCA$biplot[,3], label = rownames(rdaFL_res$CCA$biplot))) +
    labs(x=paste('RDA 1:',round(rdaFL_res.scaling1$concont$importance[2,1] * 100, 2), '%'), y=paste('RDA 3:',round(rdaFL_res.scaling1$concont$importance[2,3] * 100, 2), '%')) +
    ggtitle('RDA (full data): time ~ mcmc specification (RDA 1 & 3)') + 
    theme_bw() + 
    theme(legend.position="right", plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaFL_p2.13.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaFL_p2.13
dev.off()

rdaFL_p2.23 <- ggplot() +
    geom_point(aes(x=rdaFL_res$CCA$wa[mcmc.dump$cvg=='B',2], y=rdaFL_res$CCA$wa[mcmc.dump$cvg=='B',3]), col = "gray86") +
    geom_point(aes(x=rdaFL_res$CCA$wa[mcmc.dump$cvg=='G',2], y=rdaFL_res$CCA$wa[mcmc.dump$cvg=='G',3]), col = c('black', scales::hue_pal()(max(prd_group)))[prd_group[mcmc.dump$cvg=='G']+1]) +
    geom_segment(aes(xend=rdaFL_res$CCA$biplot[,2]*0.15, yend=rdaFL_res$CCA$biplot[,3]*0.15, x=0, y=0), colour="black", linewidth=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(aes(x=0.17*rdaFL_res$CCA$biplot[,2], y=0.17*rdaFL_res$CCA$biplot[,3], label = rownames(rdaFL_res$CCA$biplot))) +
    labs(x=paste('RDA 2:',round(rdaFL_res.scaling1$concont$importance[2,2] * 100, 2), '%'), y=paste('RDA 3:',round(rdaFL_res.scaling1$concont$importance[2,3] * 100, 2), '%')) +
    ggtitle('RDA (full data): time ~ mcmc specification (RDA 2 & 3)') + 
    theme_bw() + 
    theme(legend.position="right", plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaFL_p2.23.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaFL_p2.23
dev.off()

# Permutation test and model selection for the full dataset----
# Use permutation ANOVA to test the significance of the full RDA model,
# of each canonical axis, and of each constrained predictor.
rdaFL_perm <- anova(rdaFL_res, step=2000, parallel = ncrs)
rdaFL_perm.axis <- anova(rdaFL_res, by="axis", step=2000, parallel = ncrs)
rdaFL_perm.margin <- anova(rdaFL_res, by="margin", step=2000, parallel = ncrs)

# Record variance partitioning and predictor importance for the full model.
sink(file = 'RDAFL-margin-ANOVA-and-bidirection-step-factor-selection.txt')
cat('====rda partitioning of variance (full data)====\n')
cat(paste('Total:', rdaFL_res.scaling1$tot.chi, '1.0','\n', sep = '\t'))
cat(paste('Constrained (R.squared):', rdaFL_res.scaling1$constr.chi, rdaFL_res.scaling1$constr.chi / rdaFL_res.scaling1$tot.chi,'\n', sep = '\t'))
cat(paste('Unconstrained:', rdaFL_res.scaling1$unconst.chi, rdaFL_res.scaling1$unconst.chi / rdaFL_res.scaling1$tot.chi,'\n', sep = '\t'))
cat(paste('Adj.R.squared:', RsquareAdj(rdaFL_res)$adj.r.squared,'\n', sep = '\t'))
cat('\n\n====rda importance (full data)====\n')
rdaFL_res.scaling1$concont$importance
cat('\n\n====rda perm margin anova (full data)====\n')
rdaFL_perm.margin
cat('\n\n====rda perm bidirection-step (full data)====\n')
# Stepwise selection using both forward and backward moves, with permutation-based p-values.
rdaFL_step.both <- ordistep(rda(BrTime ~ 1, data=mcmc.dump, scale=TRUE), scope=formula(rdaFL_res), direction="both", pstep=2000, parallel = ncrs)
sink()

# RDA refined (full data)----
rdaFL_res_simple <- eval(rdaFL_step.both$call)
rdaFL_res_simple.scaling1 <- summary(rdaFL_res_simple, scaling = 1)
rdaFL_perm_simple <- anova(rdaFL_res_simple, step=2000, parallel = ncrs)
rdaFL_perm_simple.axis <- anova(rdaFL_res_simple, by="axis", step=2000, parallel = ncrs)
rdaFL_perm_simple.margin <- anova(rdaFL_res_simple, by="margin", step=2000, parallel = ncrs)
sink(file = 'RDAFL-margin-ANOVA-and-bidirection-step-factor-selection.txt', append = T)
cat('====rda partitioning of variance (simple, full data)====\n')
cat(paste('Total:', rdaFL_res_simple.scaling1$tot.chi, '1.0','\n', sep = '\t'))
cat(paste('Constrained (R.squared):', rdaFL_res_simple.scaling1$constr.chi, rdaFL_res_simple.scaling1$constr.chi / rdaFL_res_simple.scaling1$tot.chi,'\n', sep = '\t'))
cat(paste('Unconstrained:', rdaFL_res_simple.scaling1$unconst.chi, rdaFL_res_simple.scaling1$unconst.chi / rdaFL_res_simple.scaling1$tot.chi,'\n', sep = '\t'))
cat(paste('Adj.R.squared:', RsquareAdj(rdaFL_res_simple)$adj.r.squared,'\n', sep = '\t'))
cat('====rda importance (simple, full data)====\n')
rdaFL_res_simple.scaling1$concont$importance
cat('\n\n====rda perm margin anova (simple, full data)====\n')
rdaFL_perm_simple.margin
sink()

# Ordination diagram (Distance biplot)----
rdaFL_p3 <- ggplot() +
    geom_line(aes(x=c(1:length(rdaFL_res_simple$CCA$eig)), y=as.vector(rdaFL_res_simple$CCA$eig)), linetype="dashed", linewidth=1.5, color="darkgrey") +
    geom_point(aes(x=c(1:length(rdaFL_res_simple$CCA$eig)), y=as.vector(rdaFL_res_simple$CCA$eig)), size=3, color="magenta") +
    scale_x_discrete(name = "Ordination axes") +
    ylab("Inertia") +
    ggtitle('RDA loading variation (simple, full data)') + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaFL_simple_p3_loading.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaFL_p3
dev.off()

rdaFL_p4.12 <- ggplot() +
    geom_point(aes(x=rdaFL_res_simple$CCA$wa[mcmc.dump$cvg=='B',1], y=rdaFL_res_simple$CCA$wa[mcmc.dump$cvg=='B',2]), col = "gray86") +
    geom_point(aes(x=rdaFL_res_simple$CCA$wa[mcmc.dump$cvg=='G',1], y=rdaFL_res_simple$CCA$wa[mcmc.dump$cvg=='G',2]), col = c('black', scales::hue_pal()(max(prd_group)))[prd_group[mcmc.dump$cvg=='G']+1]) +
    geom_segment(aes(xend=rdaFL_res_simple$CCA$biplot[,1]*0.15, yend=rdaFL_res_simple$CCA$biplot[,2]*0.15, x=0, y=0), colour="black", linewidth=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(aes(x=0.17*rdaFL_res_simple$CCA$biplot[,1], y=0.17*rdaFL_res_simple$CCA$biplot[,2], label = rownames(rdaFL_res_simple$CCA$biplot))) +
    labs(x=paste('RDA 1:',round(rdaFL_res_simple.scaling1$concont$importance[2,1] * 100, 2), '%'), y=paste('RDA 2:',round(rdaFL_res_simple.scaling1$concont$importance[2,2] * 100, 2), '%')) +
    ggtitle('RDA (simple, full data): time ~ mcmc specification (RDA 1 & 2)') + 
    theme_bw() +
    theme(legend.position="right", plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaFL_p4.12.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaFL_p4.12
dev.off()

rdaFL_p4.13 <- ggplot() +
    geom_point(aes(x=rdaFL_res_simple$CCA$wa[mcmc.dump$cvg=='B',1], y=rdaFL_res_simple$CCA$wa[mcmc.dump$cvg=='B',3]), col = "gray86") +
    geom_point(aes(x=rdaFL_res_simple$CCA$wa[mcmc.dump$cvg=='G',1], y=rdaFL_res_simple$CCA$wa[mcmc.dump$cvg=='G',3]), col = c('black', scales::hue_pal()(max(prd_group)))[prd_group[mcmc.dump$cvg=='G']+1]) +
    geom_segment(aes(xend=rdaFL_res_simple$CCA$biplot[,1]*0.15, yend=rdaFL_res_simple$CCA$biplot[,3]*0.15, x=0, y=0), colour="black", linewidth=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(aes(x=0.17*rdaFL_res_simple$CCA$biplot[,1], y=0.17*rdaFL_res_simple$CCA$biplot[,3], label = rownames(rdaFL_res_simple$CCA$biplot))) +
    labs(x=paste('RDA 1:',round(rdaFL_res_simple.scaling1$concont$importance[2,1] * 100, 2), '%'), y=paste('RDA 3:',round(rdaFL_res_simple.scaling1$concont$importance[2,3] * 100, 2), '%')) +
    ggtitle('RDA (simple, full data): time ~ mcmc specification (RDA 1 & 3)') + 
    theme_bw() +
    theme(legend.position="right", plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaFL_p4.13.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaFL_p4.13
dev.off()

rdaFL_p4.23 <- ggplot() +
    geom_point(aes(x=rdaFL_res_simple$CCA$wa[mcmc.dump$cvg=='B',2], y=rdaFL_res_simple$CCA$wa[mcmc.dump$cvg=='B',3]), col = "gray86") +
    geom_point(aes(x=rdaFL_res_simple$CCA$wa[mcmc.dump$cvg=='G',2], y=rdaFL_res_simple$CCA$wa[mcmc.dump$cvg=='G',3]), col = c('black', scales::hue_pal()(max(prd_group)))[prd_group[mcmc.dump$cvg=='G']+1]) +
    geom_segment(aes(xend=rdaFL_res_simple$CCA$biplot[,2]*0.15, yend=rdaFL_res_simple$CCA$biplot[,3]*0.15, x=0, y=0), colour="black", linewidth=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(aes(x=0.17*rdaFL_res_simple$CCA$biplot[,2], y=0.17*rdaFL_res_simple$CCA$biplot[,3], label = rownames(rdaFL_res_simple$CCA$biplot))) +
    labs(x=paste('RDA 2:',round(rdaFL_res_simple.scaling1$concont$importance[2,2] * 100, 2), '%'), y=paste('RDA 3:',round(rdaFL_res_simple.scaling1$concont$importance[2,3] * 100, 2), '%')) +
    ggtitle('RDA (simple, full data): time ~ mcmc specification (RDA 2 & 3)') + 
    theme_bw() +
    theme(legend.position="right", plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaFL_p4.23.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaFL_p4.23
dev.off()

# RDA partitioning of variance (simple, full data)----
X1 <- mcmc.dump[,names(rdaFL_step.both$terminfo$xlev)[1]]
X2 <- mcmc.dump[,names(rdaFL_step.both$terminfo$xlev)[2]]
X3 <- mcmc.dump[,names(rdaFL_step.both$terminfo$xlev)[3]]
X4 <- mcmc.dump[,names(rdaFL_step.both$terminfo$xlev)[4]]

rdaFL_res_simple_var_part <- varpart(BrTime, X1, X2, X3, X4, data = mcmc.dump, scale=TRUE)
svg(filename = 'rdaFL_res_simple_variance_partitioning.svg', width = 8, height = 8, bg = 'transparent', pointsize = 12)
plot(rdaFL_res_simple_var_part, Xnames = names(rdaFL_step.both$terminfo$xlev)[1:4], bg = 2:5)
dev.off()

## redundancy analysis (good conv data)----
rdaST_res <- rda(BrTime[mcmc.dump$cvg=='G',] ~ ., data=mcmc.dump[mcmc.dump$cvg=='G',-1], scale=TRUE)
rdaST_res.scaling1 <- summary(rdaST_res, scaling = 1)

# Ordination diagram (Distance biplot)----
rdaST_p1 <- ggplot() + 
    geom_line(aes(x=c(1:length(rdaST_res$CCA$eig)), y=as.vector(rdaST_res$CCA$eig)), linetype="dashed", linewidth=1.5, color="darkgrey") + 
    geom_point(aes(x=c(1:length(rdaST_res$CCA$eig)), y=as.vector(rdaST_res$CCA$eig)), size=3, color="magenta") +
    scale_x_discrete(name = "Ordination axes") + 
    ylab("Inertia") + 
    ggtitle('RDA loading variation (good conv data)') + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaST_p1_loading.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaST_p1
dev.off()

rdaST_p2.12 <- ggplot() +
    geom_point(aes(x=rdaST_res$CCA$wa[,1], y=rdaST_res$CCA$wa[,2]), col = c('black', scales::hue_pal()(max(prd_group)))[prd_group[mcmc.dump$cvg=='G']+1]) +
    geom_segment(aes(xend=rdaST_res$CCA$biplot[,1]*0.15, yend=rdaST_res$CCA$biplot[,2]*0.15, x=0, y=0), colour="black", linewidth=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(aes(x=0.17*rdaST_res$CCA$biplot[,1], y=0.17*rdaST_res$CCA$biplot[,2], label = rownames(rdaST_res$CCA$biplot))) +
    labs(x=paste('RDA 1:',round(rdaST_res.scaling1$concont$importance[2,1] * 100, 2), '%'), y=paste('RDA 2:',round(rdaST_res.scaling1$concont$importance[2,2] * 100, 2), '%')) +
    ggtitle('RDA (good conv data): time ~ mcmc specification (RDA 1 & 2)') + 
    theme_bw() + 
    theme(legend.position="right", plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaST_p2.12.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaST_p2.12
dev.off()

rdaST_p2.13 <- ggplot() +
    geom_point(aes(x=rdaST_res$CCA$wa[,1], y=rdaST_res$CCA$wa[,3]), col = c('black', scales::hue_pal()(max(prd_group)))[prd_group[mcmc.dump$cvg=='G']+1]) +
    geom_segment(aes(xend=rdaST_res$CCA$biplot[,1]*0.15, yend=rdaST_res$CCA$biplot[,3]*0.15, x=0, y=0), colour="black", linewidth=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(aes(x=0.17*rdaST_res$CCA$biplot[,1], y=0.17*rdaST_res$CCA$biplot[,3], label = rownames(rdaST_res$CCA$biplot))) +
    labs(x=paste('RDA 1:',round(rdaST_res.scaling1$concont$importance[2,1] * 100, 2), '%'), y=paste('RDA 3:',round(rdaST_res.scaling1$concont$importance[2,3] * 100, 2), '%')) +
    ggtitle('RDA (good conv data): time ~ mcmc specification (RDA 1 & 3)') + 
    theme_bw() + 
    theme(legend.position="right", plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaST_p2.13.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaST_p2.13
dev.off()

rdaST_p2.23 <- ggplot() +
    geom_point(aes(x=rdaST_res$CCA$wa[,2], y=rdaST_res$CCA$wa[,3]), col = c('black', scales::hue_pal()(max(prd_group)))[prd_group[mcmc.dump$cvg=='G']+1]) +
    geom_segment(aes(xend=rdaST_res$CCA$biplot[,2]*0.15, yend=rdaST_res$CCA$biplot[,3]*0.15, x=0, y=0), colour="black", linewidth=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(aes(x=0.17*rdaST_res$CCA$biplot[,2], y=0.17*rdaST_res$CCA$biplot[,3], label = rownames(rdaST_res$CCA$biplot))) +
    labs(x=paste('RDA 2:',round(rdaST_res.scaling1$concont$importance[2,2] * 100, 2), '%'), y=paste('RDA 3:',round(rdaST_res.scaling1$concont$importance[2,3] * 100, 2), '%')) +
    ggtitle('RDA (good conv data): time ~ mcmc specification (RDA 2 & 3)') + 
    theme_bw() + 
    theme(legend.position="right", plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaST_p2.23.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaST_p2.23
dev.off()

# Permutation test (good conv data)----
rdaST_perm <- anova(rdaST_res, step=2000, parallel = ncrs)
rdaST_perm.axis <- anova(rdaST_res, by="axis", step=2000, parallel = ncrs)
rdaST_perm.margin <- anova(rdaST_res, by="margin", step=2000, parallel = ncrs)
sink(file = 'RDAST-margin-ANOVA-and-bidirection-step-factor-selection.txt')
cat('====rda partitioning of variance (good conv data)====\n')
cat(paste('Total:', rdaST_res.scaling1$tot.chi, '1.0','\n', sep = '\t'))
cat(paste('Constrained (R.squared):', rdaST_res.scaling1$constr.chi, rdaST_res.scaling1$constr.chi / rdaST_res.scaling1$tot.chi,'\n', sep = '\t'))
cat(paste('Unconstrained:', rdaST_res.scaling1$unconst.chi, rdaST_res.scaling1$unconst.chi / rdaST_res.scaling1$tot.chi,'\n', sep = '\t'))
cat(paste('Adj.R.squared:', RsquareAdj(rdaST_res)$adj.r.squared,'\n', sep = '\t'))
cat('====rda importance (good conv data)====\n')
rdaST_res.scaling1$concont$importance
cat('\n\n====rda perm margin anova (good conv data)====\n')
rdaST_perm.margin
cat('\n\n====rda perm bidirection-step (good conv data)====\n')
rdaST_step.both <- ordistep(rda(BrTime[mcmc.dump$cvg=='G',] ~ 1, data=mcmc.dump[mcmc.dump$cvg=='G',-1], scale = FALSE), scope=formula(rdaST_res), direction="both", pstep=2000, parallel = ncrs)
sink()

# RDA refined (good conv data)----
rdaST_res_simple <- eval(rdaST_step.both$call)
rdaST_res_simple.scaling1 <- summary(rdaST_res_simple, scaling = 1)
rdaST_perm_simple <- anova(rdaST_res_simple, step=2000, parallel = ncrs)
rdaST_perm_simple.axis <- anova(rdaST_res_simple, by="axis", step=2000, parallel = ncrs)
rdaST_perm_simple.margin <- anova(rdaST_res_simple, by="margin", step=2000, parallel = ncrs)
sink(file = 'RDAST-margin-ANOVA-and-bidirection-step-factor-selection.txt', append = T)
cat('====rda partitioning of variance (simple, good conv data)====\n')
cat(paste('Total:', rdaST_res_simple.scaling1$tot.chi, '1.0','\n', sep = '\t'))
cat(paste('Constrained (R.squared):', rdaST_res_simple.scaling1$constr.chi, rdaST_res_simple.scaling1$constr.chi / rdaST_res_simple.scaling1$tot.chi,'\n', sep = '\t'))
cat(paste('Unconstrained:', rdaST_res_simple.scaling1$unconst.chi, rdaST_res_simple.scaling1$unconst.chi / rdaST_res_simple.scaling1$tot.chi,'\n', sep = '\t'))
cat(paste('Adj.R.squared:', RsquareAdj(rdaST_res_simple)$adj.r.squared,'\n', sep = '\t'))
cat('====rda importance (simple, good conv data)====\n')
rdaST_res_simple.scaling1$concont$importance
cat('\n\n====rda perm margin anova (simple, good conv data)====\n')
rdaST_perm_simple.margin
sink()

# Ordination diagram (Distance biplot)----
rdaST_p3 <- ggplot() +
    geom_line(aes(x=c(1:length(rdaST_res_simple$CCA$eig)), y=as.vector(rdaST_res_simple$CCA$eig)), linetype="dashed", linewidth=1.5, color="darkgrey") +
    geom_point(aes(x=c(1:length(rdaST_res_simple$CCA$eig)), y=as.vector(rdaST_res_simple$CCA$eig)), size=3, color="magenta") +
    scale_x_discrete(name = "Ordination axes") +
    ylab("Inertia") +
    ggtitle('RDA loading variation (simple, good conv data)') + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaST_simple_p3_loading.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaST_p3
dev.off()

rdaST_p4.12 <- ggplot() +
    geom_point(aes(x=rdaST_res_simple$CCA$wa[,1], y=rdaST_res_simple$CCA$wa[,2]), col = c('black', scales::hue_pal()(max(prd_group)))[prd_group[mcmc.dump$cvg=='G']+1]) +
    geom_segment(aes(xend=rdaST_res_simple$CCA$biplot[,1]*0.15, yend=rdaST_res_simple$CCA$biplot[,2]*0.15, x=0, y=0), colour="black", linewidth=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(aes(x=0.17*rdaST_res_simple$CCA$biplot[,1], y=0.17*rdaST_res_simple$CCA$biplot[,2], label = rownames(rdaST_res_simple$CCA$biplot))) +
    labs(x=paste('RDA 1:',round(rdaST_res_simple.scaling1$concont$importance[2,1] * 100, 2), '%'), y=paste('RDA 2:',round(rdaST_res_simple.scaling1$concont$importance[2,2] * 100, 2), '%')) +
    ggtitle('RDA (simple, good conv data): time ~ mcmc specification (RDA 1 & 2)') + 
    theme_bw() +
    theme(legend.position="right", plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaST_p4.12.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaST_p4.12
dev.off()

rdaST_p4.13 <- ggplot() +
    geom_point(aes(x=rdaST_res_simple$CCA$wa[,1], y=rdaST_res_simple$CCA$wa[,3]), col = c('black', scales::hue_pal()(max(prd_group)))[prd_group[mcmc.dump$cvg=='G']+1]) +
    geom_segment(aes(xend=rdaST_res_simple$CCA$biplot[,1]*0.15, yend=rdaST_res_simple$CCA$biplot[,3]*0.15, x=0, y=0), colour="black", linewidth=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(aes(x=0.17*rdaST_res_simple$CCA$biplot[,1], y=0.17*rdaST_res_simple$CCA$biplot[,3], label = rownames(rdaST_res_simple$CCA$biplot))) +
    labs(x=paste('RDA 1:',round(rdaST_res_simple.scaling1$concont$importance[2,1] * 100, 2), '%'), y=paste('RDA 3:',round(rdaST_res_simple.scaling1$concont$importance[2,3] * 100, 2), '%')) +
    ggtitle('RDA (simple, good conv data): time ~ mcmc specification (RDA 1 & 3)') + 
    theme_bw() +
    theme(legend.position="right", plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaST_p4.13.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaST_p4.13
dev.off()

rdaST_p4.23 <- ggplot() +
    geom_point(aes(x=rdaST_res_simple$CCA$wa[,2], y=rdaST_res_simple$CCA$wa[,3]), col = c('black', scales::hue_pal()(max(prd_group)))[prd_group[mcmc.dump$cvg=='G']+1]) +
    geom_segment(aes(xend=rdaST_res_simple$CCA$biplot[,2]*0.15, yend=rdaST_res_simple$CCA$biplot[,3]*0.15, x=0, y=0), colour="black", linewidth=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(aes(x=0.17*rdaST_res_simple$CCA$biplot[,2], y=0.17*rdaST_res_simple$CCA$biplot[,3], label = rownames(rdaST_res_simple$CCA$biplot))) +
    labs(x=paste('RDA 2:',round(rdaST_res_simple.scaling1$concont$importance[2,2] * 100, 2), '%'), y=paste('RDA 3:',round(rdaST_res_simple.scaling1$concont$importance[2,3] * 100, 2), '%')) +
    ggtitle('RDA (simple, good conv data): time ~ mcmc specification (RDA 2 & 3)') + 
    theme_bw() +
    theme(legend.position="right", plot.title = element_text(hjust = 0.5))
svg(filename = 'rdaST_p4.23.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
rdaST_p4.23
dev.off()

# RDA partitioning of variance (simple, good conv data)----
# Decompose explained variance among the selected predictors from stepwise RDA.
# This clarifies how much variation is uniquely attributable to each predictor,
# and how much is shared across predictors.
X1 <- mcmc.dump[,names(rdaST_step.both$terminfo$xlev)[1]]
X2 <- mcmc.dump[,names(rdaST_step.both$terminfo$xlev)[2]]
X3 <- mcmc.dump[,names(rdaST_step.both$terminfo$xlev)[3]]
X4 <- mcmc.dump[,names(rdaST_step.both$terminfo$xlev)[4]]

rdaST_res_simple_var_part <- varpart(BrTime, X1, X2, X3, X4, data = mcmc.dump, scale = FALSE)
svg(filename = 'rdaST_res_simple_variance_partitioning.svg', width = 8, height = 8, bg = 'transparent', pointsize = 12)
plot(rdaST_res_simple_var_part, Xnames = names(rdaST_step.both$terminfo$xlev)[1:4], bg = 2:5)
dev.off()

rm(X1,X2,X3,X4)

# Convert generated SVG figures to PDF for downstream reporting.
system('bash ./svg2pdf.sh r')

# Save the full workspace so these RDA results can be reloaded later.
save.image('/To/Your/Directory/nematoda/timing/AllPosteriorTimeTree.cLGPMSF.RData')
