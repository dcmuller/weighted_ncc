# p01_sim_nomatch.r
# dcmuller
# Simulation study of weighted unconditional analysis of an NCC study
# with no matching factors
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
version_randsurv <- "0.1"
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
  colnames(res) <- c("intercept", "exposure", "log_scale")
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
##############################
## function: sim_no_match()
##############################
# mega function to sample the cohort and the NCC study, and return
# the results of each analysis. -params- is a list of parameters
# for the simulation.
sim_no_match <- function(params) {
  exposure <- as.matrix(rnorm(params$n))
  compet_exposure <- as.matrix(rnorm(params$n))
  cohort <- weibull_compet2(lambda0 = params$baseline_rate,
                            p = params$weib_param,
                            ## time is relative to the origin
                            cens_time = params$cens - params$origin,
                            X_1 = exposure,
                            beta_1 = params$lhr,
                            X_2 = compet_exposure,
                            beta_2 = params$compet_lhr)
  ## time is relative to the origin:
  cohort$obs_time <- cohort$obs_time + params$origin
  cohort$failed <- as.numeric(cohort$status==1)
  cohort$exposure <- exposure
  # full cohort analysis
  st_full <- Surv(time    = cohort$obs_time, 
                  event   = cohort$failed, 
                  origin  = params$origin)
  full <- survreg(st_full ~ cohort$exposure, dist = "weibull")
  full <- coef_se(full)
  coxfull <- coxph(st_full ~ cohort$exposure)
  coxfull <- data.frame(intercept = rep(NA ,2), 
                        exposure = c(coef(coxfull), sqrt(diag(vcov(coxfull)))),
                        log_scale = rep(NA, 2),
                        param = c("coef", "se"),
                        stringsAsFactors=FALSE)
  
  # sample ncc
  ncc <- ncc_sample(entry = 0, 
                    exit = obs_time, 
                    fail = failed, 
                    controls = params$ncc_controls, 
                    include = list(exposure, obs_time),
                    data = cohort,
                    silent = TRUE)
  # some people might have been sampled as controls more than once, 
  # keep only one record
  setkey(ncc, ncc_id, ncc_fail)
  ncc <- as.data.frame(ncc[, tail(.SD, 1), by = c("ncc_id")])
    
  # generate survival object
  st <- Surv(time   = ncc$obs_time, 
             event  = ncc$ncc_fail,
             origin = params$origin)
  weighted <- survreg(st ~ exposure + cluster(ncc_id),  
                      weights = 1/ncc_pr,
                      robust = TRUE,
                      data = ncc) 
  weighted <- coef_se(weighted)
  coxweighted <- coxph(st ~ exposure + cluster(ncc_id),  
                       weights = 1/ncc_pr,
                       robust = TRUE,
                      data = ncc) 
  coxweighted <- data.frame(intercept = rep(NA, 2),
                            exposure = c(coef(coxweighted), sqrt(diag(vcov(coxweighted)))),
                            log_scale = rep(NA, 2),
                            param = c("coef", "se"),
                            stringsAsFactors = FALSE)
  
  full$model <- "full"
  weighted$model<- "weighted"
  coxfull$model <- "coxfull"
  coxweighted$model <- "coxweighted"
  res <- rbindlist(list(full, weighted, coxfull, coxweighted))
  return(res)
}

# set up parameters that are constant accross scenarios
cohort_size <- 20000
cens <- 80
origin <- 30 # time at which people become at risk
baseline_rate_event <- 5e-6
baseline_rate_compet <- 10e-7
weib_param_event <- 2.5
weib_param_compet <- 3.4
ncc_controls <- 2

# place sets of parameters in a list so we can use lapply
sim_settings <- list(
  # hr of 2; baseline rate 5-e7; weibull param 2.5 (rare disease)
  list(n = cohort_size, 
       cens = cens,
       baseline_rate = c(5e-7, baseline_rate_compet),
       weib_param = c(2.5, weib_param_compet),
       origin = origin, 
       ncc_controls = ncc_controls,
       compet_lhr = log(1),
       lhr = c(log(2))),
  # hr of 2; baseline rate 10-e5; weibull param 1.8 (common disease, e.g. PCa)
  list(n = cohort_size, 
       cens = cens,
       baseline_rate = c(10e-5, baseline_rate_compet),
       weib_param = c(1.8, weib_param_compet),
       origin = origin, 
       ncc_controls = ncc_controls,
       compet_lhr = log(1),
       lhr = c(log(2)))
)
  
##############################
# Run simulation
##############################
nsims <- 500
# seed chosen by sampling a random integer between 1 and 9999
# at <www.random.org> 
set.seed(3211)
output <- vector(mode = "list", length(sim_settings))
params <- vector(mode = "list", length(sim_settings))
for (i in 1:length(sim_settings)) {
  cat(paste0("\nSimulating for parameter set ", i, ":\n"))
  simres <- vector(mode = "list", nsims)
  for (j in 1:nsims) {
    simres[[j]] <- tryCatch(sim_no_match(sim_settings[[i]]),
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
write.csv(out, file = "./analysis/output/o01_all_samples.csv")
write.csv(par, file = "./analysis/output/o01_params.csv")
