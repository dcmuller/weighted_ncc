# p03_sim_weakmatch.r
# dcmuller
# Simulation study of weighted unconditional analysis of an NCC study
# with a matching variable strongly correlated with the exposure
# 
# Compares three potential analyses:
# 1) theoretical parametric analysis of the full cohort
# 2) conditional analysis of the nested case-control sample
# 3) weighted parametric analysis of the nested case-control sample

# survival package for weibull regression,
# devtools to install utilities from github
# randsurv to simulate survival times
# ncctools to generate ncc data with sampling weights
library(survival)
library(devtools)
library(data.table)
library(ggplot2)
version_randsurv <- "0.2.1"
version_ncctools <- "0.2"
if ("randsurv" %in% rownames(installed.packages()) == FALSE) {
  install_github("randsurv", 
                 username="dcmuller", 
                 ref=paste0("v", version_randsurv))
}
if (installed.packages()["randsurv", "Version"] != version_randsurv) {
  warning(paste("replacing randsurv with version", version_randsurv))
  install_github("randsurv", 
                 username="dcmuller", 
                 ref=paste0("v", version_randsurv))
}
if ("ncctools" %in% rownames(installed.packages()) == FALSE) {
  install_github("ncctools", 
                 username="dcmuller", 
                 ref=paste0("v", version_ncctools))
}
if (installed.packages()["ncctools", "Version"] != version_ncctools) {
  warning(paste("replacing ncctools with version", version_ncctools))
  install_github("ncctools", 
                 username="dcmuller", 
                 ref=paste0("v", version_ncctools))
}
library(randsurv)
library(ncctools)


##############################
## function: coef_se()
##############################
# function to return estimated parameters and their standard errors
# from a survreg object
coef_se <- function(mod) {
  if (class(mod) != "survreg") {
    stop("Fitted model must be of class survreg")
  }
  res <- as.data.frame(rbind(c(mod$coefficients, log(mod$scale)),
                             c(sqrt(diag(mod$var)))))
  colnames(res) <- c("intercept", 
                     names(mod$coefficients)[-1],
                     "log_scale")
  res$param <- c("coef", "se")
  return(res)
}

##############################
## function: dots()
##############################
# function to display dots to show progress of, and provide
# entertainment during, simulations
dots <- function(iter) {
  di_str <- ifelse(((iter %% 10) == 0) & !(iter %% 50 == 0),
                   paste0(iter),
                   ifelse((iter %% 50) == 0, paste(iter, "\n"), "."))
  return(di_str)
}

################################
## function: gen_binary_match()
################################
# function to generate a binary matching factor that is
# associated with the exposure
gen_binary_match <- function(exposure, exposure_or, prob) {
  logodds0 <- log(prob) - log(1-prob)
  logor <- log(exposure_or)
  xb <- logodds0 + logor*exposure
  pr <- exp(xb)/(1 + exp(xb))
  match <- rbinom(length(pr), 1, pr)
  match
}

###############################
## function: sim_strong_match()
###############################
# mega function to sample the cohort and the NCC study, and return
# the results of each analysis. -params- is a list of parameters
# for the simulation.
sim_strong_match <- function(params) {
  exposure <- rnorm(params$n)
  match <- gen_binary_match(exposure, 
                            exposure_or=params$matchvar_exposure_or,
                            prob=params$prev_matchvar)
  X <- cbind(exposure = exposure, matchvar = match)
  compet_exposure <- as.matrix(rnorm(params$n))
  cohort <- weibull_compet(lambda0 = params$baseline_rate,
                           p = params$weib_param,
                           ## time is relative to the origin
                           cens_time = params$cens - params$origin,
                           X_1 = X,
                           beta_1 = params$lhr,
                           X_2 = compet_exposure,
                           beta_2 = params$compet_lhr)
  ## time is relative to the origin:
  cohort$obs_time <- cohort$obs_time + params$origin
  cohort$f_time <- cohort$f_time + params$origin
  cohort$failed <- as.numeric(cohort$status==1)
  cohort <- cbind(cohort, X)
  # full cohort analysis
  st_full <- Surv(time    = cohort$obs_time, 
                  event   = cohort$failed, 
                  origin  = params$origin)
  exposure <- cohort$exposure
  matchvar <- cohort$matchvar
  full <- survreg(st_full ~ exposure + matchvar, dist = "weibull")
  full <- coef_se(full)
  coxfull <- coxph(st_full ~ exposure + matchvar)
  coxfull <- data.frame(intercept = rep(NA ,2), 
                        rbind(coef(coxfull), sqrt(diag(vcov(coxfull)))),
                        log_scale = rep(NA, 2),
                        param = c("coef", "se"),
                        stringsAsFactors=FALSE)
  
  # sample ncc
  ncc <- ncc_sample(entry = 0,
                    exit = obs_time, 
                    fail = failed,
                    match = matchvar,
                    controls = params$ncc_controls, 
                    include = list(exposure, obs_time),
                    data = cohort,
                    silent = TRUE)
  
  # fit conditional logistic regression
  cond <- clogit(ncc_fail ~ exposure + strata(ncc_set),
                 data=ncc)
  cond <- data.frame(intercept = rep(NA, 2),
                     exposure = c(coef(cond), sqrt(diag(vcov(cond)))),
                     matchvar = rep(NA, 2),
                     log_scale = rep(NA, 2),
                     param = c("coef", "se"),
                     stringsAsFactors = FALSE)
  # some people might have been sampled as controls more than once, 
  # keep only one record for the weighted analyses
  setkey(ncc, ncc_id, ncc_fail)
  ncc <- as.data.frame(ncc[, tail(.SD, 1), by = c("ncc_id")])
  
  # generate survival object
  st <- Surv(time   = ncc$obs_time, 
             event  = ncc$ncc_fail,
             origin = params$origin)
  weighted <- survreg(st ~ exposure + matchvar + cluster(ncc_id),  
                      weights = 1/ncc_pr,
                      robust = TRUE,
                      data = ncc) 
  weighted <- coef_se(weighted)
  coxweighted <- coxph(st ~ exposure + matchvar + cluster(ncc_id),  
                       weights = 1/ncc_pr,
                       robust = TRUE,
                       data = ncc) 
  coxweighted <- data.frame(intercept = rep(NA, 2),
                            rbind(coef(coxweighted), sqrt(diag(vcov(coxweighted)))),
                            log_scale = rep(NA, 2),
                            param = c("coef", "se"),
                            stringsAsFactors = FALSE)
  
  cond$model <- "conditional"
  full$model <- "full"
  weighted$model<- "weighted"
  coxfull$model <- "coxfull"
  coxweighted$model <- "coxweighted"
  res <- rbindlist(list(full, cond, weighted, coxfull, coxweighted))
  return(res)
}

## set up parameters that are constant accross scenarios
cohort_size <- 20000
cens <- 80
origin <- 30 # time at which people become at risk
baseline_rate_event <- 5e-6
baseline_rate_compet <- 10e-7
weib_param_event <- 2.5
weib_param_compet <- 3.4
ncc_controls <- 2
hr_exposure <- 2
hr_match <- 4
# odds ratio for assoc between matching factor and exposure
matchvar_exposure_or <- 1.5 
# proportion of people who have matchvar==1 in the population
prev_matchvar <- 0.5

# place sets of parameters in a list so we can use lapply
sim_settings <- list(
  # baseline rate 5-e7; weibull param 2.5 (rare disease)
  list(n = cohort_size, 
       cens = cens,
       baseline_rate = c(5e-7, baseline_rate_compet),
       weib_param = c(2.5, weib_param_compet),
       origin = origin, 
       ncc_controls = ncc_controls,
       compet_lhr = log(1),
       prev_matchvar = prev_matchvar,
       matchvar_exposure_or = matchvar_exposure_or,
       lhr = c(log(hr_exposure), log(hr_match))),
       
  # baseline rate 10-e5; weibull param 1.8 (common disease, e.g. PCa)
  list(n = cohort_size, 
       cens = cens,
       baseline_rate = c(10e-5, baseline_rate_compet),
       weib_param = c(1.8, weib_param_compet),
       origin = origin, 
       ncc_controls = ncc_controls,
       compet_lhr = log(1),
       prev_matchvar = prev_matchvar,
       matchvar_exposure_or = matchvar_exposure_or,
       lhr = c(log(hr_exposure), log(hr_match)))
)
  
##############################
# Run simulation
##############################
nsims <- 100
# seed chosen by sampling a random integer between 1 and 9999
# at <www.random.org> 
set.seed(3211)
output <- vector(mode = "list", length(sim_settings))
params <- vector(mode = "list", length(sim_settings))
for (i in 1:length(sim_settings)) {
  cat(paste0("\nSimulating for parameter set ", i, ":\n"))
  simres <- vector(mode = "list", nsims)
  for (j in 1:nsims) {
    simres[[j]] <- tryCatch(sim_strong_match(sim_settings[[i]]),
                       warning = function(e) {
                         warning("there was a warning")
                         return(NULL)
                       },
                       error = function(e) {
                         warning("there was an error")
                         return(NULL)
                       })
    cat(dots(j))
  }
  output[[i]] <- rbindlist(simres) 
  output[[i]]$settings <- i
  ## add simulation parameters
  params[[i]] <- as.data.frame(t(unlist(sim_settings[[i]])))
  params[[i]]$settings <- i
}
out <- rbindlist(output)
par <- rbindlist(params)
write.csv(out, file = "./analysis/output/o03_all_samples.csv")
write.csv(par, file = "./analysis/output/o03_params.csv")
