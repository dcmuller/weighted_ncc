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
  colnames(res) <- c(names(mod$coefficients), "Log(scale)")
  res$param <- c("coef", "se")
  return(res)
}

##############################
## function: dots()
##############################
# function to display dots to show progress of, and provide
# entertainment during, simulations
dots <- function(iter) {
  di_str <- ifelse((iter %% 10) == 0,
                   paste0(iter, "\n"),
                   ifelse((iter %% 50) == 0,
                          paste(iter),
                          "."))
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
                            cens_time = params$cens,
                            X_1 = exposure,
                            beta_1 = params$lhr,
                            X_2 = compet_exposure,
                            beta_2 = params$compet_lhr)
  cohort$failed <- as.numeric(cohort$status==1)
  cohort$exposure <- exposure
  # full cohort analysis
  full <- survreg(Surv(time = cohort$obs_time, event = cohort$failed) ~ exposure, 
                  dist = "weibull")
  full <- coef_se(full)
  
  # sample ncc
  ncc <- ncc_sample(entry = 0, 
                    exit = obs_time, 
                    fail = failed, 
                    controls = params$ncc_controls, 
                    include = list(exposure, obs_time),
                    data = cohort,
                    silent = TRUE)
  ncc <- as.data.frame(ncc)
  # generate survival object
  st <- Surv(time = ncc$obs_time, event = ncc$ncc_fail)
  weighted <- survreg(st ~ exposure,  
                      weights = 1/ncc_pr,
                      data = ncc) 
  weighted <- coef_se(weighted)
  
  full$model <- "full"
  weighted$model<- "weighted"
  res <- rbindlist(list(full, weighted))
  iter <- ifelse(exists("iter", inherits=FALSE), iter+1, 1)
  cat(dots(iter))
  return(res)
}

# set up parameters that are constant accross scenarios
cohort_size <- 10000
cens <- 80
baseline_rate_event <- 10e-10
baseline_rate_compet <- 10e-9
weib_param_event <- 4.3
weib_param_compet <- 4.2
ncc_controls <- 2

# place sets of parameters in a list so we can use lapply
sim_settings <- list(
  # hr of 2
  list(n = cohort_size, 
       cens = cens,
       baseline_rate = c(baseline_rate_event, baseline_rate_compet),
       weib_param = c(weib_param_event, weib_param_compet),
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
for(i in 1:length(sim_settings)) {
  cat(paste0("Simulating for parameter set ", i, ":\n"))
  output[[i]] <- replicate(nsims, 
                           sim_no_match(sim_settings[[i]]), 
                           simplify=FALSE)
}
out <- rbindlist(output[[1]])
write.csv(out, file="./analysis/output/o01_sim_nomatch.csv")
