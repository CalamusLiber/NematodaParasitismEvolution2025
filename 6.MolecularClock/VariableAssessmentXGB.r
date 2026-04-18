# ===============================================================================
# VariableAssessmentXGB.r
# Purpose: Build and interpret an XGBoost classifier for convergence status
#   (good vs bad convergence), then explain feature contributions using SHAP.
# ===============================================================================

## Loading libraries----
# xgboost: gradient boosting model fitting
# tidyverse: data wrangling and plotting helpers
# DataExplorer: exploratory data report utilities
# pROC: ROC curve and AUC evaluation
# shapr / iml / SHAPforxgboost: SHAP explanation tools
# tidymodels: tuning workflow and resampling
# DiagrammeR: optional model/graph visualization
# doParallel: parallel backend control for caret/tidymodels
# svglite: write vector SVG plots reliably
library(xgboost)
library(tidyverse)
library(DataExplorer)
library(pROC)
library(shapr)
library(iml)
library(tidymodels)
library(DiagrammeR)
library(SHAPforxgboost)
library(doParallel)
library(svglite)

# Set working directory for plot and log outputs.
setwd('/To/Your/Directory/nematoda/timing/Res_posttrees')

# Custom helper: compute SHAP contributions from a fitted XGBoost model.
# Returns raw per-sample SHAP values, mean absolute SHAP scores, and bias term.
shap.value.xgb <- function (xgb_model, X_train, target_class = NULL) {
    shap_contrib <- predict(xgb_model, (X_train), predcontrib = TRUE)
    if (is.list(shap_contrib)) {
        shap_contrib <- if (!is.null(target_class) && !is.na(target_class))
            shap_contrib[[target_class + 1]]
        else Reduce("+", lapply(shap_contrib, abs))
    }
    if (is.null(colnames(shap_contrib))) {
        colnames(shap_contrib) <- c(colnames(X_train), "BIAS")
    } else {
        colnames(shap_contrib)[colnames(shap_contrib) == '(Intercept)'] <- "BIAS"
    }
    shap_contrib <- data.table::as.data.table(shap_contrib)
    BIAS0 <- shap_contrib[, ncol(shap_contrib), with = FALSE][1]
    shap_contrib[, `:=`(BIAS, NULL)]
    imp <- colMeans(abs(shap_contrib))
    mean_shap_score <- imp[order(imp, decreasing = TRUE)]
    return(list(shap_score = shap_contrib, mean_shap_score = mean_shap_score,
                BIAS0 = BIAS0))
}

## Factors that determine the convergence (XGBoost classification)----
# Split the binary convergence dataset into training and validation partitions.
# We use an 80/20 partition and preserve the class distribution of the target.
set.seed(20260404)
XGB_train2 <- as.vector(createDataPartition(y = mcmc.dump.bin[, 1], p = 0.8, list = FALSE))

# Create XGBoost DMatrix objects for faster training and evaluation.
XGB_data_train2 <- xgb.DMatrix(data = as.matrix(mcmc.dump.bin[XGB_train2, -1]), label = as.matrix(mcmc.dump.bin[XGB_train2, 1]))
XGB_data_valid2 <- xgb.DMatrix(data = as.matrix(mcmc.dump.bin[-XGB_train2, -1]), label = as.matrix(mcmc.dump.bin[-XGB_train2, 1]))

watchlist2 <- list(train = XGB_data_train2, test = XGB_data_valid2)

# Disable parallel backend for tidymodels tuning to avoid conflicts.
stopImplicitCluster()
registerDoSEQ()

# XGBoost hyperparameter tuning using tidymodels and repeated cross-validation.
# The function returns a native xgboost parameter list for xgb.train().
tune_xgb_params <- function(data, target_var, v_folds = 5, repeats = 3, grid_size = 5) {
    # Disable TODO ALTRET to avoid known xgboost pointer errors on some R builds.
    Sys.setenv("XGBOOST_ALTREP_OPTIMIZATION" = "FALSE")

    # Ensure the target class is a factor for tidymodels preprocessing.
    data[[target_var]] <- as.factor(data[[target_var]])

    # Build a preprocessing recipe for categorical encoding and zero-variance filtering.
    xgb_recipe <- recipe(as.formula(paste(target_var, "~ .")), data = data) %>%
        step_novel(all_nominal_predictors()) %>%
        step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
        step_zv(all_predictors())

    # Define a tunable XGBoost classification model.
    # All potential tuning parameters must be wrapped with tune().
    xgb_spec <- boost_tree(
        trees = 2000,
        tree_depth = tune(),
        min_n = tune(),
        loss_reduction = tune(),
        sample_size = tune(),
        mtry = tune(),
        learn_rate = tune()
    ) %>%
        set_engine("xgboost", nthread = 1) %>%
        set_mode("classification")

    # Create repeated stratified folds for robust tuning.
    set.seed(20260404)
    folds <- vfold_cv(
        data,
        v = v_folds,
        repeats = repeats,
        strata = !!sym(target_var)
    )

    # Generate a space-filling tuning grid for the chosen hyperparameters.
    xgb_grid <- grid_space_filling(
        tree_depth(),
        min_n(),
        loss_reduction(),
        sample_size = sample_prop(),
        finalize(mtry(), data),
        learn_rate(),
        size = grid_size
    )

    message(paste0("Starting Repeated CV: ", v_folds, "-fold, ", repeats, " repeats..."))

    xgb_workflow <- workflow() %>%
        add_recipe(xgb_recipe) %>%
        add_model(xgb_spec)

    tune_results <- tune_grid(
        xgb_workflow,
        resamples = folds,
        grid = xgb_grid,
        metrics = metric_set(accuracy, roc_auc),
        control = control_grid(verbose = FALSE, allow_par = FALSE)
    )

    # Choose the best tuning candidate by classification accuracy.
    best_params_raw <- select_best(tune_results, metric = "accuracy")

    # Convert tidymodels parameter names to native xgboost parameter names.
    num_features <- ncol(data) - 1
    native_params <- list(
        booster          = "gbtree",
        objective        = "binary:logistic",
        eta              = best_params_raw$learn_rate,
        max_depth        = best_params_raw$tree_depth,
        min_child_weight = best_params_raw$min_n,
        gamma            = best_params_raw$loss_reduction,
        subsample        = best_params_raw$sample_size,
        colsample_bytree = min(best_params_raw$mtry / num_features, 1)
    )
    return(native_params)
}

bestTune <- tune_xgb_params(data = mcmc.dump.bin, target_var = 'cvg', v_folds = 5, repeats = 3, grid_size = 20)

# Model training using the tuned hyperparameters from tidymodels.
params <- xgb.params(
    booster = 'gbtree',
    eta = bestTune$eta,
    gamma = bestTune$gamma,
    max_depth = bestTune$max_depth,
    min_child_weight = bestTune$min_child_weight,
    subsample = bestTune$subsample,
    colsample_bytree = bestTune$colsample_bytree,
    objective = 'binary:logistic',
    seed = 20260404
)

sink(file = 'XGB-model-2-class-training-logging.txt')
cat("===== custom train model =====\n\n")
XGB_mod2 <- xgb.train(
    data = XGB_data_train2,
    params = params,
    nrounds = 2000,
    evals = watchlist2,
    verbose = 1,
    print_every_n = 100
)
XGB_mod2
cat("\n\n===== best Tune by tidymodels =====\n\n")
bestTune
sink()

# Save variable importance scores to text and plot the top features.
sink(file = 'XGB-model-2-class-variable-importance-Gain.txt')
VarImp_xgb <- xgb.importance(model = XGB_mod2)
print(VarImp_xgb)
sink()

# Plot the top 8 important variables by Gain.
svg(filename = 'XGB-model-2-class-variable-importance-Gain.svg', width = 8, height = 16, bg = 'white', pointsize = 12)
xgb.plot.importance(importance_matrix = VarImp_xgb, measure = 'Gain', top_n = 8)
dev.off()

# Create a long-format dataset for multi-metric importance plotting.
VarImp_xgb_all4 <- VarImp_xgb %>% 
    select(Feature, Gain, Cover, Frequency) %>% 
    pivot_longer(cols = 2:4, names_to = 'Measure', values_to = 'Value') %>% 
    mutate(
        Angle = -360 * (as.numeric(factor(Feature)) - 0.5) / nrow(VarImp_xgb),
        Angle = ifelse(Angle < -90 & Angle > -270, Angle + 180, Angle)
    )

# Radial bar plot showing multiple importance measures for each feature.
imp <- ggplot(VarImp_xgb_all4, aes(Feature, Value, fill = Feature)) + 
    geom_col(width = 1, position = "identity") + 
    geom_text(aes(label = Feature, y = 1, angle = Angle),
              color = "black", size = 6) + 
    geom_text(aes(label = round(Value, 4), y = Value + 0.1, angle = Angle),
              color = "black", size = 4) + 
    coord_radial(start = 0, theta = "x", clip = "off", r.axis.inside = FALSE, inner.radius = 0.2) + 
    scale_y_continuous(limits = NULL) +
    facet_wrap(~Measure, scales = "free_y") + 
    scale_fill_manual(values = c("#00BFFFFF", "#B22222FF", "#227733FF", "#DAA520FF", "#ED3B4FFF", "#6BBF23FF", "#9932CCFF",  "#708090FF")) + 
    theme_void() + 
    theme(axis.text.x = element_blank(),
          strip.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom") + 
    ggtitle("XGB-model Variable Importance")

svg(filename = 'XGB-model-2-class-variable-importance-all4.svg', width = 18, height = 6, bg = 'white', pointsize = 12)
imp
dev.off()

# Compute SHAP values to explain prediction contributions for each feature.
shap_values_2 <- shap.value.xgb(xgb_model = XGB_mod2, X_train = as.matrix(mcmc.dump.bin[XGB_train2, -1]))

svg(filename = 'XGB-model-2-class-SHAP-dependence.svg', width = 8, height = 16, bg = 'white', pointsize = 12)
xgb.plot.shap(data = as.matrix(mcmc.dump.bin[XGB_train2, -1]), model = XGB_mod2, top_n = 8, n_col = 1, pch = 16, cex = 0.8)
dev.off()

# Use svglite for pure SVG output and avoid bitmap fallback issues.
svglite(filename = 'XGB-model-2-class-SHAP-summary-1.svg', width = 8, height = 6, bg = 'white', pointsize = 12)
xgb.ggplot.shap.summary(data = as.matrix(mcmc.dump.bin[XGB_train2, -1]), model = XGB_mod2, top_n = 8)
dev.off()

# Generate a second SHAP summary plot using the shap.prep pipeline.
svglite(filename = 'XGB-model-2-class-SHAP-summary-2.svg', width = 8, height = 6, bg = 'transparent', pointsize = 12)
shap_long_2 <- shap.prep(xgb_model = XGB_mod2, X_train = as.matrix(mcmc.dump.bin[XGB_train2, -1]), shap_contrib = shap_values_2$shap_score)
shap.plot.summary(shap_long_2, kind = 'sina')
dev.off()

# Optional SHAP dependence and force plots are commented out for now.
# plot_list_2 <- lapply(c('part', 'corr', 'cstLS', 'fpFull', 'cstLM', 'fpOI', 'bd', 'rg'), shap.plot.dependence, data_long = shap_long_2)
# gridExtra::grid.arrange(grobs = plot_list_2, ncol = 1)
# plot_data_2 <- shap.prep.stack.data(shap_contrib = shap_values_2$shap_score, top_n = 4, n_groups = 4)
# shap.plot.force_plot(plot_data_2, y_parent_limit = c(-8,8))
# shap.plot.force_plot_bygroup(plot_data_2)

# Generate predicted probabilities on the training set and compute ROC metrics.
XGB_pred2 <- data.frame(Prob = predict(XGB_mod2, newdata = XGB_data_train2))
train_roc2 <- roc(response = mcmc.dump.bin[XGB_train2, 1],
                  predictor = XGB_pred2[,1])
# Area under the curve: 0.9459
bestp <- train_roc2$thresholds[which.max(train_roc2$sensitivities + train_roc2$specificities - 1)]
XGB_pred2$Conv <- factor(ifelse(XGB_pred2$Prob > bestp, 'G', 'B'))

svg(filename = 'XGB-model-2-class-ROC-curve.svg', width = 8, height = 8, bg = 'white', pointsize = 12)
plot(train_roc2, print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, max.auc.polygon = TRUE, auc.polygon.col = 'skyblue', print.thres = TRUE, legacy.axes = TRUE, bty = 'l')
dev.off()

# Write classification performance metrics to text files.
sink(file = 'XGB-model-2-class-Confusion-Matrix.txt')
confusionMatrix(data = XGB_pred2$Conv,
                reference = factor(mcmc.dump[XGB_train2, 1]),
                positive = 'G',
                mode = 'everything')
sink()

sink(file = 'XGB-model-2-class-Binary-Summary.txt')
defaultSummary(
    data.frame(obs = factor(mcmc.dump[XGB_train2, 1]),
               pred = XGB_pred2$Conv),
    lev = levels(factor(mcmc.dump[XGB_train2, 1]))
)
sink()

# Convert all generated SVG output to PDF for reporting.
system('bash ./svg2pdf.sh r')

# Save the workspace so the full XGBoost model and SHAP results can be restored.
save.image('/To/Your/Directory/nematoda/timing/AllPosteriorTimeTree.cLGPMSF.RData')
