# ===============================================================================
# HiSSEModels.r
# Purpose: Custom function to fit a suite of HiSSE (Hidden State Speciation and Extinction)
#   models for ancestral state reconstruction, including model selection, pruning,
#   and marginal ancestral state estimation.
# ===============================================================================

# HiSSEModels function: Fits multiple HiSSE model variants, selects the best,
# performs ancestral reconstruction, and returns model-averaged rates.
# Arguments:
#   phy: phylogenetic tree
#   data: trait data with tip labels
#   f: sampling fractions
#   keep.all.models: whether to keep all fitted models
#   criterion: 'AIC' or 'AICc' for model selection
#   multiprocessors: number of cores for parallel computation
HiSSEModels <- function(phy, data, f, keep.all.models=FALSE, criterion='AIC', multiprocessors=8){
    # Nested function: PruneRedundantModelsX
    # Removes models with identical log-likelihoods (within precision) to avoid redundancy.
    # Useful for HiSSE, MuHiSSE, GeoHiSSE, or MiSSE model lists.
    PruneRedundantModelsX <- function(..., precision = 1e-05)
    {
        models <- list(...)
        if (!inherits(models[[1]], what = c("geohisse.fit", "hisse.fit", "misse.fit", "muhisse.fit"))) {
            models <- models[[1]]
        }
        mod.class.geohisse <- sapply(models, function(x) inherits(x, what = "geohisse.fit"))
        mod.class.hisse <- sapply(models, function(x) inherits(x, what = "hisse.fit"))
        mod.class.misse <- sapply(models, function(x) inherits(x, what = "misse.fit"))
        mod.class.muhisse <- sapply(models, function(x) inherits(x, what = "muhisse.fit"))
        if (all(mod.class.geohisse) & all(mod.class.hisse) & all(mod.class.misse) & 
            all(mod.class.muhisse)) {
            stop("list of models need to be only HiSSE, MuHiSSE, GeoHiSSE, or MiSSE fits.")
        }
        if (!all(mod.class.geohisse) & !all(mod.class.hisse) & !all(mod.class.misse) & 
            !all(mod.class.muhisse)) {
            stop("list of models need to be only HiSSE, MuHiSSE, GeoHiSSE, or MiSSE fits.")
        }
        mod.nparameters <- simplify2array(lapply(lapply(models, "[[", "starting.vals"), length))
        models <- models[order(mod.nparameters, decreasing = FALSE)]
        mod.loglik <- simplify2array(lapply(models, "[[", "loglik"))
        models_to_delete <- c()
        isTrueAllEqual <- function(...) {
            return(isTRUE(all.equal(...)))
        }
        if (length(models) > 1) {
            for (i in 2:(length(models))) {
                if (any(sapply(mod.loglik[1:(i - 1)], isTrueAllEqual, 
                               mod.loglik[i], tolerance = precision))) {
                    models_to_delete <- c(models_to_delete, i)
                }
            }
        }
        return(models_to_delete)
    }
    
    MODELS <- NULL
    Conf <- NULL

    # ===============================================================================
    # Model Configurations
    # Defines a set of HiSSE model variants with different parameter constraints.
    # Types:
    # - BiSSE: Binary State Speciation and Extinction (no hidden states)
    # - CID2/CID4: Character-Independent Diversification (2 or 4 hidden states)
    # - HiSSE: Hidden State Speciation and Extinction (2 hidden states)
    # Each config specifies turnover (tau), extinction fraction (eps), transition rates (tr),
    # hidden states flag (HiState), and description (Acnt).
    # ===============================================================================

	# mod01: Yule BiSSE model (t0=t1, e0=e1=0, q's equal, np=2)
	tau <- c(1,1)
	eps <- c(0,0)
	tr <- TransMatMakerHiSSE(hidden.traits=0)
	tr[!is.na(tr)] <- 1
	HiState <- FALSE
	Acnt <- "Yule BiSSE (t0=t1, e0=e1=0, q's equal)"
	Conf[[1]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
    # mod02: Yule BiSSE model (t0=t1, e0=e1=0, q's free, np=3)
	tau <- c(1,1)
	eps <- c(0,0)
	tr <- TransMatMakerHiSSE(hidden.traits=0)
	HiState <- FALSE
	Acnt <- "Yule BiSSE (t0=t1, e0=e1=0, q's free)"
	Conf[[2]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod03: Yule BiSSE model (t0, t1, e0=e1=0, q's equal, np=3)
	tau <- c(1,2)
	eps <- c(0,0)
	tr <- TransMatMakerHiSSE(hidden.traits=0)
	tr[!is.na(tr)] <- 1
	HiState <- FALSE
	Acnt <- "Yule BiSSE (e0=e1=0, q's equal)"
	Conf[[3]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
    # mod04: Yule BiSSE model (t0, t1, e0=e1=0, q's free, np=4)
	tau <- c(1,2)
	eps <- c(0,0)
	tr <- TransMatMakerHiSSE(hidden.traits=0)
	HiState <- FALSE
	Acnt <- "Yule BiSSE (e0=e1=0, q's free)"
	Conf[[4]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)

    # mod05: normal BiSSE model (t0=t1, e0=e1, q's equal, np=3)
	tau <- c(1,1)
	eps <- c(1,1)
	tr <- TransMatMakerHiSSE(hidden.traits=0)
	tr[!is.na(tr)] <- 1
	HiState <- FALSE
	Acnt <- "BiSSE (t0=t1, e0=e1, q's equal)"
	Conf[[5]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
    # mod06: normal BiSSE model (t0=t1, e0=e1, q's free, np=4)
	tau <- c(1,1)
	eps <- c(1,1)
	tr <- TransMatMakerHiSSE(hidden.traits=0)
	HiState <- FALSE
	Acnt <- "BiSSE (t0=t1, e0=e1, q's free)"
	Conf[[6]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod07: normal BiSSE model (t0, t1, e0=e1, q's equal, np=4)
	tau <- c(1,2)
	eps <- c(1,1)
	tr <- TransMatMakerHiSSE(hidden.traits=0)
	tr[!is.na(tr)] <- 1
	HiState <- FALSE
	Acnt <- "BiSSE (e0=e1, q's equal)"
	Conf[[7]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
    # mod08: normal BiSSE model (t0, t1, e0=e1, q's free, np=5)
	tau <- c(1,2)
	eps <- c(1,1)
	tr <- TransMatMakerHiSSE(hidden.traits=0)
	HiState <- FALSE
	Acnt <- "BiSSE (e0=e1, q's free)"
	Conf[[8]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod09: normal BiSSE model (t0, t1, e0, e1, q's equal, np=5)
	tau <- c(1,2)
	eps <- c(1,2)
	tr <- TransMatMakerHiSSE(hidden.traits=0)
	tr[!is.na(tr)] <- 1
	HiState <- FALSE
	Acnt <- "BiSSE (q's equal)"
	Conf[[9]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod10: normal BiSSE model (t0, t1, e0, e1, q's free, np=6)
	tau <- c(1,2)
	eps <- c(1,2)
	tr <- TransMatMakerHiSSE(hidden.traits=0)
	HiState <- FALSE
	Acnt <- "BiSSE (q's free)"
	Conf[[10]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod11: CID2 HiSSE model (t0A=t1A, t0B=t1B, e's equal, q's equal, np=4)
	tau <- c(1,1,2,2)
	eps <- c(1,1,1,1)
	tr <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "CID2 (t0A=t1A, t0B=t1B, e's equal, q's equal)"
	Conf[[11]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod12: CID2 HiSSE model (t0A=t1A, t0B=t1B, e's equal, q0A0B=q0B0A=q1A1B=q1B1A, np=6)
	tau <- c(1,1,2,2)
	eps <- c(1,1,1,1)
	tr <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
	HiState <- TRUE
	Acnt <- "CID2 (t0A=t1A, t0B=t1B, e's equal, q0A0B=q0B0A=q1A1B=q1B1A)"
	Conf[[12]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod13: CID2 HiSSE model (t0A=t1A, t0B=t1B, e0A=e1A, e0B=e1B, q's equal, np=5)
	tau <- c(1,1,2,2)
	eps <- c(1,1,2,2)
	tr <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "CID2 (t0A=t1A, t0B=t1B, e0A=e1A, e0B=e1B, q's equal)"
	Conf[[13]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod14: CID2 HiSSE model (t0A=t1A, t0B=t1B, e0A=e1A, e0B=e1B, q0A1A=q0B1B, q1A0A=q1B0B, all other q's equal, np=7)
	tau <- c(1,1,2,2)
	eps <- c(1,1,2,2)
	tr <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
	HiState <- TRUE
	Acnt <- "CID2 (t0A=t1A, t0B=t1B, e0A=e1A, e0B=e1B, q0A1A=q0B1B, q1A0A=q1B0B, all other q's equal)"
	Conf[[14]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod15: CID4 HiSSE model (t0A=t1A, t0B=t1B, t0C=t1C, t0D=t1D, e's equal, q's equal, np=6)
	tau <- c(1,1,2,2,3,3,4,4)
	eps <- rep(1, 8)
	tr <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "CID4 (t0A=t1A, t0B=t1B, t0C=t1C, t0D=t1D, e's equal, q's equal)"
	Conf[[15]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod16: CID4 HiSSE model (t0A=t1A, t0B=t1B, t0C=t1C, t0D=t1D, e's equal, q0A1A=q0B1B=q0C1C=q0D1D, q1A0A=q1B0B=q1C0C=q1D0D, all other q's equal, np=8)
	tau <- c(1,1,2,2,3,3,4,4)
	eps <- rep(1, 8)
	tr <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
	HiState <- TRUE
	Acnt <- "CID4 (t0A=t1A, t0B=t1B, t0C=t1C, t0D=t1D, e's equal, q0A1A=q0B1B=q0C1C=q0D1D, q1A0A=q1B0B=q1C0C=q1D0D, all other q's equal)"
	Conf[[16]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod17: CID4 HiSSE model (t0A=t1A, t0B=t1B, t0C=t1C, t0D=t1D, e0A=e1A, e0B=e1B, e0C=e1C, e0D=e1D, q's equal, np=9)
	tau <- c(1,1,2,2,3,3,4,4)
	eps <- c(1,1,2,2,3,3,4,4)
	tr <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "CID4 (t0A=t1A, t0B=t1B, t0C=t1C, t0D=t1D, e0A=e1A, e0B=e1B, e0C=e1C, e0D=e1D, q's equal)"
	Conf[[17]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod18: CID4 HiSSE model (t0A=t1A, t0B=t1B, t0C=t1C, t0D=t1D, e0A=e1A, e0B=e1B, e0C=e1C, e0D=e1D, q0A1A=q0B1B=q0C1C=q0D1D, q1A0A=q1B0B=q1C0C=q1D0D, all other q's equal, np=11)
	tau <- c(1,1,2,2,3,3,4,4)
	eps <- c(1,1,2,2,3,3,4,4)
	tr <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
	HiState <- TRUE
	Acnt <- "CID4 (t0A=t1A, t0B=t1B, t0C=t1C, t0D=t1D, e0A=e1A, e0B=e1B, e0C=e1C, e0D=e1D, q0A1A=q0B1B=q0C1C=q0D1D, q1A0A=q1B0B=q1C0C=q1D0D, all other q's equal)"
	Conf[[18]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod19: character-dependent HiSSE model (t's equal, e's equal, q's equal, np=3)	
	tau <- rep(1, 4)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t's equal, e's equal, q's equal)"
	Conf[[19]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod20: character-dependent HiSSE model (t's equal, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=7)	
	tau <- rep(1, 4)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (t's equal, e's equal, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[20]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod21: character-dependent HiSSE model (t's equal, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=6)	
	tau <- rep(1, 4)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (t's equal, e's equal, q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[21]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod22: character-dependent HiSSE model (t's free, e's equal, q's equal, np=6)	
	tau <- c(1:4)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (e's equal, q's equal)"
	Conf[[22]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)		
	
	# mod23: character-dependent HiSSE model (t's free, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=10)	
	tau <- c(1:4)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (e's equal, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[23]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)		
	
	# mod24: character-dependent HiSSE model (t's free, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=9)	
	tau <- c(1:4)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (e's equal, q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[24]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)		

	# mod25: character-dependent HiSSE model (t's free, e's equal, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A, np=8)	
	tau <- c(1:4)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	HiState <- TRUE
	Acnt <- "HiSSE (e's equal, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[25]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)		

	# mod26: character-dependent HiSSE model (t's free, e's equal, q0B1B=q1B0B=0, all other q's equal, np=6)	
	tau <- c(1:4)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	tr[tr > 0] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (e's equal, q0B1B=q1B0B=0, all other q's equal)"
	Conf[[26]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod27: character-dependent HiSSE model (t's free, e's free, q's equal, np=9)	
	tau <- c(1:4)
	eps <- c(1:4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (q's equal)"
	Conf[[27]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod28: character-dependent HiSSE model (t's free, e's free, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=13)	
	tau <- c(1:4)
	eps <- c(1:4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[28]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)		
	
	# mod29: character-dependent HiSSE model (t's free, e's free, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=12)	
	tau <- c(1:4)
	eps <- c(1:4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[29]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)		

	# mod30: character-dependent HiSSE model (t's free, e's free, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A, np=11)	
	tau <- c(1:4)
	eps <- c(1:4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	HiState <- TRUE
	Acnt <- "HiSSE (q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[30]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod31: character-dependent HiSSE model (t's free, e's free, q0B1B=q1B0B=0, all other q's equal, np=9)	
	tau <- c(1:4)
	eps <- c(1:4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	tr[tr > 0] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (q0B1B=q1B0B=0, all other q's equal)"
	Conf[[31]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)		

	# mod32: character-dependent HiSSE model (t0A=t1A=t0B, e0A=e1A=e0B, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=9)	
	tau <- c(1,1,1,2)
	eps <- c(1,1,1,2)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A=t0B, e0A=e1A=e0B, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[32]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)		
	
	# mod33: character-dependent HiSSE model (t0A=t1A=t0B, e0A=e1A=e0B, q's equal, np=5)	
	tau <- c(1,1,1,2)
	eps <- c(1,1,1,2)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A=t0B, e0A=e1A=e0B, q's equal)"
	Conf[[33]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod34: character-dependent HiSSE model (t0A=t1A=t0B, e0A=e1A=e0B, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=8)	
	tau <- c(1,1,1,2)
	eps <- c(1,1,1,2)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A=t0B, e0A=e1A=e0B, q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[34]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod35: character-dependent HiSSE model (t0A=t1A=t0B, e0A=e1A=e0B, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A, np=7)	
	tau <- c(1,1,1,2)
	eps <- c(1,1,1,2)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A=t0B, e0A=e1A=e0B, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[35]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)		
	
	# mod36: character-dependent HiSSE model (t0A=t1A=t0B, e0A=e1A=e0B, q0B1B=q1B0B=0, all other q's equal, np=5)	
	tau <- c(1,1,1,2)
	eps <- c(1,1,1,2)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	tr[tr > 0] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A=t0B, e0A=e1A=e0B, q0B1B=q1B0B=0, all other q's equal)"
	Conf[[36]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)		

	# mod37: character-dependent HiSSE model (t0A=t1A=t0B, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=8)	
	tau <- c(1,1,1,2)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A=t0B, e's equal, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[37]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod38: character-dependent HiSSE model (t0A=t1A=t0B, e's equal, q's equal, np=4)	
	tau <- c(1,1,1,2)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A=t0B, e's equal, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[38]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod39: character-dependent HiSSE model (t0A=t1A=t0B, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=7)	
	tau <- c(1,1,1,2)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A=t0B, e's equal, q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[39]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod40: character-dependent HiSSE model (t0A=t1A=t0B, e's equal, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A, np=6)	
	tau <- c(1,1,1,2)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A=t0B, e's equal, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[40]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)		
	
	# mod41: character-dependent HiSSE model (t0A=t1A=t0B, e's equal, q0B1B=q1B0B=0, all other q's equal, np=4)	
	tau <- c(1,1,1,2)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	tr[tr > 0] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A=t0B, e's equal, q0B1B=q1B0B=0, all other q's equal)"
	Conf[[41]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)		

	# mod42: character-dependent HiSSE model (t0A=t1A, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=9)	
	tau <- c(1,1,2,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A, e's equal, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[42]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod43: character-dependent HiSSE model (t0A=t1A, e's equal, q's equal, np=5)	
	tau <- c(1,1,2,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A, e's equal, q's equal)"
	Conf[[43]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod44: character-dependent HiSSE model (t0A=t1A, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=8)	
	tau <- c(1,1,2,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A, e's equal, q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[44]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod45: character-dependent HiSSE model (t0A=t1A, e's equal, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A, np=7)	
	tau <- c(1,1,2,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A, e's equal, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[45]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod46: character-dependent HiSSE model (t0A=t1A, e's equal, q0B1B=q1B0B=0, all other q's equal, np=5)	
	tau <- c(1,1,2,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	tr[tr > 0] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A, e's equal, q0B1B=q1B0B=0, all other q's equal)"
	Conf[[46]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	

	# mod47: character-dependent HiSSE model (t0A=t1A, e0A=e1A, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=11)	
	tau <- c(1,1,2,3)
	eps <- c(1,1,2,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A, e0A=e1A, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[47]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod48: character-dependent HiSSE model (t0A=t1A, e0A=e1A, q's equal, np=7)	
	tau <- c(1,1,2,3)
	eps <- c(1,1,2,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A, e0A=e1A, q's equal)"
	Conf[[48]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod49: character-dependent HiSSE model (t0A=t1A, e0A=e1A, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=10)	
	tau <- c(1,1,2,3)
	eps <- c(1,1,2,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A, e0A=e1A, q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[49]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod50: character-dependent HiSSE model (t0A=t1A, e0A=e1A, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A, np=9)	
	tau <- c(1,1,2,3)
	eps <- c(1,1,2,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A, e0A=e1A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[50]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod51: character-dependent HiSSE model (t0A=t1A, e0A=e1A, q0B1B=q1B0B=0, all other q's equal, np=7)	
	tau <- c(1,1,2,3)
	eps <- c(1,1,2,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	tr[tr > 0] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t1A, e0A=e1A, q0B1B=q1B0B=0, all other q's equal)"
	Conf[[51]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)

	# mod52: character-dependent HiSSE model (t0B=t1B, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=9)	
	tau <- c(1,2,3,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (t0B=t1B, e's equal, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[52]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod53: character-dependent HiSSE model (t0B=t1B, e's equal, q's equal, np=5)	
	tau <- c(1,2,3,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0B=t1B, e's equal, q's equal)"
	Conf[[53]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod54: character-dependent HiSSE model (t0B=t1B, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=8)	
	tau <- c(1,2,3,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (t0B=t1B, e's equal, q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[54]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod55: character-dependent HiSSE model (t0B=t1B, e's equal, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A, np=7)	
	tau <- c(1,2,3,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	HiState <- TRUE
	Acnt <- "HiSSE (t0B=t1B, e's equal, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[55]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod56: character-dependent HiSSE model (t0B=t1B, e's equal, q0B1B=q1B0B=0, all other q's equal, np=5)	
	tau <- c(1,2,3,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	tr[tr > 0] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0B=t1B, e's equal, q0B1B=q1B0B=0, all other q's equal)"
	Conf[[56]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)

	# mod57: character-dependent HiSSE model (t0B=t1B, e0B=e1B, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=11)	
	tau <- c(1,2,3,3)
	eps <- c(1,2,3,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (t0B=t1B, e0B=e1B, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[57]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod58: character-dependent HiSSE model (t0B=t1B, e0B=e1B, q's equal, np=7)	
	tau <- c(1,2,3,3)
	eps <- c(1,2,3,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0B=t1B, e0B=e1B, q's equal)"
	Conf[[58]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod59: character-dependent HiSSE model (t0B=t1B, e0B=e1B, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=10)	
	tau <- c(1,2,3,3)
	eps <- c(1,2,3,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (t0B=t1B, e0B=e1B, q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[59]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod60: character-dependent HiSSE model (t0B=t1B, e0B=e1B, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A, np=9)	
	tau <- c(1,2,3,3)
	eps <- c(1,2,3,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	HiState <- TRUE
	Acnt <- "HiSSE (t0B=t1B, e0B=e1B, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[60]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod61: character-dependent HiSSE model (t0B=t1B, e0B=e1B, q0B1B=q1B0B=0, all other q's equal, np=7)	
	tau <- c(1,2,3,3)
	eps <- c(1,2,3,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	tr[tr > 0] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0B=t1B, e0B=e1B, q0B1B=q1B0B=0, all other q's equal)"
	Conf[[61]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)

	# mod62: character-dependent HiSSE model (t0A=t0B, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=9)	
	tau <- c(1,2,1,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t0B, e's equal, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[62]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod63: character-dependent HiSSE model (t0A=t0B, e's equal, q's equal, np=5)	
	tau <- c(1,2,1,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t0B, e's equal, q's equal)"
	Conf[[63]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod64: character-dependent HiSSE model (t0A=t0B, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=8)	
	tau <- c(1,2,1,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t0B, e's equal, q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[64]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod65: character-dependent HiSSE model (t0A=t0B, e's equal, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A, np=7)	
	tau <- c(1,2,1,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t0B, e's equal, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[65]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod66: character-dependent HiSSE model (t0A=t0B, e's equal, q0B1B=q1B0B=0, all other q's equal, np=5)	
	tau <- c(1,2,1,3)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	tr[tr > 0] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t0B, e's equal, q0B1B=q1B0B=0, all other q's equal)"
	Conf[[66]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)

	# mod67: character-dependent HiSSE model (t0A=t0B, e0A=e0B, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=11)	
	tau <- c(1,2,1,3)
	eps <- c(1,2,1,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t0B, e0A=e0B, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[67]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod68: character-dependent HiSSE model (t0A=t0B, e0A=e0B, q's equal, np=7)	
	tau <- c(1,2,1,3)
	eps <- c(1,2,1,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t0B, e0A=e0B, q's equal)"
	Conf[[68]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod69: character-dependent HiSSE model (t0A=t0B, e0A=e0B, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=10)	
	tau <- c(1,2,1,3)
	eps <- c(1,2,1,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t0B, e0A=e0B, q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[69]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod70: character-dependent HiSSE model (t0A=t0B, e0A=e0B, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A, np=9)	
	tau <- c(1,2,1,3)
	eps <- c(1,2,1,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t0B, e0A=e0B, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[70]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod71: character-dependent HiSSE model (t0A=t0B, e0A=e0B, q0B1B=q1B0B=0, all other q's equal, np=7)	
	tau <- c(1,2,1,3)
	eps <- c(1,2,1,3)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	tr[tr > 0] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t0A=t0B, e0A=e0B, q0B1B=q1B0B=0, all other q's equal)"
	Conf[[71]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	

	# mod72: character-dependent HiSSE model (t1A=t1B, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=9)	
	tau <- c(1,2,3,2)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (t1A=t1B, e's equal, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[72]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod73: character-dependent HiSSE model (t1A=t1B, e's equal, q's equal, np=5)	
	tau <- c(1,2,3,2)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t1A=t1B, e's equal, q's equal)"
	Conf[[73]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod74: character-dependent HiSSE model (t1A=t1B, e's equal, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=8)	
	tau <- c(1,2,3,2)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (t1A=t1B, e's equal, q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[74]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod75: character-dependent HiSSE model (t1A=t1B, e's equal, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A, np=7)	
	tau <- c(1,2,3,2)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	HiState <- TRUE
	Acnt <- "HiSSE (t1A=t1B, e's equal, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[75]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod76: character-dependent HiSSE model (t1A=t1B, e's equal, q0B1B=q1B0B=0, all other q's equal, np=5)	
	tau <- c(1,2,3,2)
	eps <- rep(1, 4)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	tr[tr > 0] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t1A=t1B, e's equal, q0B1B=q1B0B=0, all other q's equal)"
	Conf[[76]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	

	# mod77: character-dependent HiSSE model (t1A=t1B, e1A=e1B, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A, np=11)	
	tau <- c(1,2,3,2)
	eps <- c(1,2,3,2)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	HiState <- TRUE
	Acnt <- "HiSSE (t1A=t1B, e1A=e1B, q0A0B=q1A1B=q0B0A=q1B1)"
	Conf[[77]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod78: character-dependent HiSSE model (t1A=t1B, e1A=e1B, q's equal, np=7)	
	tau <- c(1,2,3,2)
	eps <- c(1,2,3,2)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr[!is.na(tr)] <- 1
	HiState <- TRUE
	Acnt <- "HiSSE (t1A=t1B, e1A=e1B, q's equal)"
	Conf[[78]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod79: character-dependent HiSSE model (t1A=t1B, e1A=e1B, q0A1A, q1A0A, q0B1B, q1B0B, q0A0B=q1A1B=q0B0A=q1B1A=0, np=10)	
	tau <- c(1,2,3,2)
	eps <- c(1,2,3,2)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, 5)
	HiState <- TRUE
	Acnt <- "HiSSE (t1A=t1B, e1A=e1B, q0A0B=q1A1B=q0B0A=q1B1A=0)"
	Conf[[79]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)	
	
	# mod80: character-dependent HiSSE model (t1A=t1B, e1A=e1B, q0A1A, q1A0A, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A, np=9)	
	tau <- c(1,2,3,2)
	eps <- c(1,2,3,2)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	HiState <- TRUE
	Acnt <- "HiSSE (t1A=t1B, e1A=e1B, q0B1B=q1B0B=0, q0A0B=q1A1B=q0B0A=q1B1A)"
	Conf[[80]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# mod81: character-dependent HiSSE model (t1A=t1B, e1A=e1B, q0B1B=q1B0B=0, all other q's equal, np=7)	
	tau <- c(1,2,3,2)
	eps <- c(1,2,3,2)
	tr <- TransMatMakerHiSSE(hidden.traits=1)
	tr <- ParDrop(tr, c(3, 4))
	tr[tr > 0] <- 1
    HiState <- TRUE
	Acnt <- "HiSSE (t1A=t1B, e1A=e1B, q0B1B=q1B0B=0, all other q's equals)"
	Conf[[81]] <- list(tau=tau, eps=eps, trate=tr, HiState=HiState, Acnt=Acnt)
	
	# ===============================================================================
	# Model Fitting
	# Fit all 81 HiSSE model configurations in parallel using foreach.
	# ===============================================================================
	library(doParallel)  
    cl <- makeCluster(multiprocessors, type="PSOCK")  
    registerDoParallel(cl)  
    MODELS <- foreach(i=1:81, .combine = list, .multicombine = TRUE, .maxcombine = 200, .packages = c('ape', 'hisse')) %dopar% {
	    hisse(phy=phy, data=data, f=f, turnover=Conf[[i]]$tau, eps=Conf[[i]]$eps, hidden.states=Conf[[i]]$HiState, trans.rate=Conf[[i]]$trate)
	}
	
	stopCluster(cl)
	
	cat('HiSSE modelling finished.\n')
	
	# ===============================================================================
	# Model Selection and Ancestral Reconstruction
	# Compile model selection table, prune redundant models, and perform marginal
	# ancestral state estimation with model averaging.
	# ===============================================================================
	
	# Create a data frame summarizing all fitted models.
	ModelSelection <- data.frame(Model=rep(NA, 81), Account=rep(NA, 81), lnlik=rep(NA, 81), n.parameters=rep(NA, 81), AIC=rep(NA, 81), AICc=rep(NA, 81))

	for (j in 1:81){
	    ModelSelection[j, 'Model'] <- j
	    ModelSelection[j, 'Account'] <- Conf[[j]]$Acnt
	    ModelSelection[j, 'lnlik'] <- MODELS[[j]]$loglik
	    ModelSelection[j, 'n.parameters'] <- length(MODELS[[j]]$starting.vals)
	    ModelSelection[j, 'AIC'] <- MODELS[[j]]$AIC
	    ModelSelection[j, 'AICc'] <- MODELS[[j]]$AICc
	}
    
	# Prune models with identical log-likelihoods.
	models_to_delete <- PruneRedundantModelsX(MODELS)
	
    if (keep.all.models){
        # If keeping all models, sort by criterion and compute AIC weights.
        ModelSelectionSorted <- ModelSelection[with(ModelSelection,order(ModelSelection[,criterion])),]
        ModelSelectionSorted[,ifelse(criterion=='AICc', 'dAICc', 'dAIC')] <- ModelSelectionSorted[,criterion]-ModelSelectionSorted[1,criterion]
        ModelSelectionSorted[,'w'] <- exp(-0.5*ModelSelectionSorted[,ifelse(criterion=='AICc', 'dAICc', 'dAIC')])/sum(exp(-0.5*ModelSelectionSorted[,ifelse(criterion=='AICc', 'dAICc', 'dAIC')]))
        
        MODELS <- MODELS[ModelSelectionSorted[,'Model']]
        
        cat('HiSSE sorting finished.\n')
        
        # Perform marginal ancestral state reconstruction for all models.
        Pred.AllModels <- list()
        for (k in 1:length(MODELS)){
            Pred.AllModels[[k]] <- MarginReconHiSSE(phy=phy, data=data, f=f, pars=MODELS[[k]]$solution, hidden.states=nrow(MODELS[[k]]$trans.matrix)/2, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils = FALSE, k.samples = NULL, AIC=MODELS[[k]][[criterion]], get.tips.only=FALSE, verbose=TRUE, n.cores=multiprocessors, dt.threads=1)
        }
        
        cat('HiSSE marginal ancestral state estimation finished.\n')
        
        # Compute model-averaged rates across all models.
        ModelAveRates <- GetModelAveRates(Pred.AllModels, AIC.weights=ModelSelectionSorted[,'w'], type="nodes")
        
        cat('HiSSE all-model averaging finished.\n')
        
        # Compute rates for the best model only.
        ModelAveRates.BestModel <- GetModelAveRates(Pred.AllModels[[1]], AIC.weights=NULL, type="both")
        
        cat('HiSSE best-model averaging finished.\n')
        
        # Optional: Adaptive sampling of likelihood surface (commented out).
        # SupportRegion.BestModel <- SupportRegionHiSSE(MODELS[[1]])
        
        return(list(ModelSelection=ModelSelectionSorted, AllModels=MODELS, ModelAveRates=ModelAveRates, Pred.AllModels=Pred.AllModels, BestModel=MODELS[[1]], ModelAveRates.BestModel=ModelAveRates.BestModel))
        
    }else if (length(models_to_delete) > 0){
        # If pruning redundant models, remove them and proceed with refined set.
        ModelSelection <- ModelSelection[-models_to_delete, ]
        
        cat('HiSSE pruning finished.\n')
        
        # Sort remaining models and compute weights.
        ModelSelectionSorted <- ModelSelection[with(ModelSelection,order(ModelSelection[,criterion])),]
        ModelSelectionSorted[,ifelse(criterion=='AICc', 'dAICc', 'dAIC')] <- ModelSelectionSorted[,criterion]-ModelSelectionSorted[1,criterion]
        ModelSelectionSorted[,'w'] <- exp(-0.5*ModelSelectionSorted[,ifelse(criterion=='AICc', 'dAICc', 'dAIC')])/sum(exp(-0.5*ModelSelectionSorted[,ifelse(criterion=='AICc', 'dAICc', 'dAIC')]))       
        RefinedModels <- MODELS[ModelSelectionSorted[,'Model']]
        
        cat('HiSSE sorting finished.\n')
        
        # Ancestral reconstruction for refined models.
        Pred.RefinedModels <- list()
        for (k in 1:length(RefinedModels)){
            Pred.RefinedModels[[k]] <- MarginReconHiSSE(phy=phy, data=data, f=f, pars=RefinedModels[[k]]$solution, hidden.states=nrow(RefinedModels[[k]]$trans.matrix)/2, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils = FALSE, k.samples = NULL, AIC=RefinedModels[[k]][[criterion]], get.tips.only=FALSE, verbose=TRUE, n.cores=multiprocessors, dt.threads=1)
        }
        
        cat('HiSSE marginal ancestral state estimation finished.\n')
        
        # Model-averaged rates for refined models.
        ModelAveRates <- GetModelAveRates(Pred.RefinedModels, AIC.weights=ModelSelectionSorted[,criterion], type="both")
        
        cat('HiSSE all-model averaging finished.\n')
        
        # Best model rates.
        ModelAveRates.BestModel <- GetModelAveRates(Pred.RefinedModels[[1]], AIC.weights=NULL, type="both")
        
        cat('HiSSE best-model averaging finished.\n')
        
        # Optional: Likelihood surface sampling.
        # SupportRegion.BestModel <- SupportRegionHiSSE(RefinedModels[[1]])
        
        return(list(ModelSelection=ModelSelectionSorted, RefinedModels=RefinedModels, ModelAveRates=ModelAveRates, Pred.RefinedModels=Pred.RefinedModels, BestModel=RefinedModels[[1]], ModelAveRates.BestModel=ModelAveRates.BestModel))
	}
}
