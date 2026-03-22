#!/usr/bin/env Rscript

# 1. Library paths (FIRST)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

# 2. Load libraries
library(coda)
library(R2jags)
library(doParallel)
library(doRNG)


# 3. Define functions

#parallel in Trillium
run_one_sim <- function(seed,nburn = 5000, niter = 10000,
                        N = 1000, t0 = 4) {
  set.seed(seed)
  # Your simulation code here
  #nburn:Burn-in iterations; niter:number of posterior draws after burn-in
  #N: sample size; M:inner MC size; t0: predicted time
  
  #------------------------------------------------------#
  #     JAGS models for all three approaches.            #
  #------------------------------------------------------#
  
  jags_model_latU <- "
  model {
  for (i in 1:N) {
    
    # ----- Linear predictors: Event (D) and Censor (C)
    muD[i] <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*A[i] + eta_u*U[i]
    log_bD[i] <- -k_D * muD[i]
    bD[i] <- exp(log_bD[i])
    
    muC[i] <- zeta0 + zeta_x1*X1[i] + zeta_x2*X2[i] + zeta_a*A[i] + zeta_u*U[i]
    log_bC[i] <- -k_C * muC[i]
    bC[i] <- exp(log_bC[i])
    
    # ----- Log-likelihood parts at observed time T[i]
    # Event D
    logfD[i] <- log_bD[i] + log(k_D) + (k_D - 1) * log(T[i]) - bD[i] * pow(T[i], k_D)
    logSD[i] <- - bD[i] * pow(T[i], k_D)
    
    # Censor C
    logfC[i] <- log_bC[i] + log(k_C) + (k_C - 1) * log(T[i]) - bC[i] * pow(T[i], k_C)
    logSC[i] <- - bC[i] * pow(T[i], k_C)
    
    
    # Observed-data log-likelihood:
    # d=1:            log f_D(T) + log S_C(T)
    # d=0 & !admin:   log f_C(T) + log S_D(T)          # random censoring at T < tau
    # d=0 & admin:    log S_C(tau) + log S_D(tau)      # adminstrative censoring at tau
    logL[i] <- d[i] * (logfD[i] + logSC[i]) +
      (1 - d[i]) * ( (1 - admin[i]) * (logfC[i] + logSD[i]) + admin[i] * ( - bC[i] * pow(tau, k_C)  - bD[i] * pow(tau, k_D) ))
    
    
    # zeros trick
    lambda[i] <- -logL[i] + Cconst
    zeros[i] ~ dpois(lambda[i])
    
    # ----- Treatment model
    A[i] ~ dbern(pA[i])
    logit(pA[i]) <- alpha0 + alpha1*X1[i] + alpha2*X2[i] + alpha3*U[i]
    
    # ----- Latent U model
    U[i] ~ dbern(pU[i])
    logit(pU[i]) <- gamma0 + gamma1*X1[i] + gamma2*X2[i] 
  }
  
  # ----- Priors 
  # Event
  eta0   ~ dnorm(0, 0.1)
  eta_x1 ~ dnorm(0, 0.25)
  eta_x2 ~ dnorm(0, 0.25)
  eta_a  ~ dnorm(0, 0.25)
  eta_u  ~ dunif(-3, 3)
  
  # Censoring
  zeta0   ~ dnorm(0, 0.1)
  zeta_x1 ~ dnorm(0, 0.25)
  zeta_x2 ~ dnorm(0, 0.25)
  zeta_a  ~ dnorm(0, 0.25)
  zeta_u  ~ dunif(-3, 3)
  
  
  # U model
  gamma0 ~ dunif(-5, 5)
  gamma1 ~ dunif(-3, 3)
  gamma2 ~ dunif(-3, 3)
  
  # Treatment model
  alpha0 ~ dnorm(0, 0.1)
  alpha1 ~ dnorm(0, 0.25)
  alpha2 ~ dnorm(0, 0.25)
  alpha3 ~ dunif(-3, 3)
  
  # Shapes
  k_D ~ dgamma(9, 6)   # event shape
  k_C ~ dgamma(6, 5)  # censoring shape
  
  # Derived quantity;
  # Dirichlet weights for N observations
  
  for (i in 1:N) {
    alpha_w[i] <- 1.0
  }
  
  w[1:N] ~ ddirch(alpha_w[])
  # Survival at t0 under A=1 and A=0
  t0k <- pow(t0, k_D)
  
  for (i in 1:N) {
    # A=1
    muD1[i]    <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*1 + eta_u*U[i]
    log_bD1[i] <- -k_D * muD1[i]
    bD1[i]     <- exp(log_bD1[i])
    S1_ind[i]  <- exp(-bD1[i] * t0k)
  
    # A=0
    muD0[i]    <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*0 + eta_u*U[i]
    log_bD0[i] <- -k_D * muD0[i]
    bD0[i]     <- exp(log_bD0[i])
    S0_ind[i]  <- exp(-bD0[i] * t0k)
  }
 
  # Marginal survival and SPCE
  S1_marg <- inprod(w[], S1_ind[])
  S0_marg <- inprod(w[], S0_ind[])
  SPCE    <- S1_marg - S0_marg
}"
  
  jags_model_naiveU <- "
  model {
  for (i in 1:N) {
    
    # ----- Linear predictors: Event (D) and Censor (C)
    muD[i] <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*A[i] 
    log_bD[i] <- -k_D * muD[i]
    bD[i] <- exp(log_bD[i])
    
    muC[i] <- zeta0 + zeta_x1*X1[i] + zeta_x2*X2[i] + zeta_a*A[i] 
    log_bC[i] <- -k_C * muC[i]
    bC[i] <- exp(log_bC[i])
    
    # ----- Log-likelihood parts at observed time T[i]
    # Event D
    logfD[i] <- log_bD[i] + log(k_D) + (k_D - 1) * log(T[i]) - bD[i] * pow(T[i], k_D)
    logSD[i] <- - bD[i] * pow(T[i], k_D)
    
    # Censor C
    logfC[i] <- log_bC[i] + log(k_C) + (k_C - 1) * log(T[i]) - bC[i] * pow(T[i], k_C)
    logSC[i] <- - bC[i] * pow(T[i], k_C)
    
    
    # Observed-data log-likelihood:
    # d=1:            log f_D(T) + log S_C(T)
    # d=0 & !admin:   log f_C(T) + log S_D(T)          # random censoring at T < tau
    # d=0 & admin:    log S_C(tau) + log S_D(tau)      # adminstrative censoring at tau
    logL[i] <- d[i] * (logfD[i] + logSC[i]) +
      (1 - d[i]) * ( (1 - admin[i]) * (logfC[i] + logSD[i]) + admin[i] * ( - bC[i] * pow(tau, k_C)  - bD[i] * pow(tau, k_D) ))
    
    
    # zeros trick
    lambda[i] <- -logL[i] + Cconst
    zeros[i] ~ dpois(lambda[i])

  }
  
  # ----- Priors 
  # Event
  eta0   ~ dnorm(0, 0.1)
  eta_x1 ~ dnorm(0, 0.25)
  eta_x2 ~ dnorm(0, 0.25)
  eta_a  ~ dnorm(0, 0.25)

  # Censoring
  zeta0   ~ dnorm(0, 0.1)
  zeta_x1 ~ dnorm(0, 0.25)
  zeta_x2 ~ dnorm(0, 0.25)
  zeta_a  ~ dnorm(0, 0.25)

  # Shapes
  k_D ~ dgamma(9, 6)   # event shape
  k_C ~ dgamma(6, 5)  # censoring shape
  
  # Derived quantity;
  # Dirichlet weights for N observations
  
  for (i in 1:N) {
    alpha_w[i] <- 1.0
  }
  
  w[1:N] ~ ddirch(alpha_w[])
  # Survival at t0 under A=1 and A=0
  t0k <- pow(t0, k_D)
  
  for (i in 1:N) {
    # A=1
    muD1[i]    <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*1
    log_bD1[i] <- -k_D * muD1[i]
    bD1[i]     <- exp(log_bD1[i])
    S1_ind[i]  <- exp(-bD1[i] * t0k)
    
    # A=0
    muD0[i]    <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*0 
    log_bD0[i] <- -k_D * muD0[i]
    bD0[i]     <- exp(log_bD0[i])
    S0_ind[i]  <- exp(-bD0[i] * t0k)
  }
  
  # Marginal survival and SPCE
  S1_marg <- inprod(w[], S1_ind[])
  S0_marg <- inprod(w[], S0_ind[])
  SPCE    <- S1_marg - S0_marg
}"
  
  jags_model_obsU <- "
  model {
  for (i in 1:N) {
    
    # ----- Linear predictors: Event (D) and Censor (C)
    muD[i] <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*A[i] + eta_u*U[i]
    log_bD[i] <- -k_D * muD[i]
    bD[i] <- exp(log_bD[i])
    
    muC[i] <- zeta0 + zeta_x1*X1[i] + zeta_x2*X2[i] + zeta_a*A[i] + zeta_u*U[i]
    log_bC[i] <- -k_C * muC[i]
    bC[i] <- exp(log_bC[i])
    
    # ----- Log-likelihood parts at observed time T[i]
    # Event D
    logfD[i] <- log_bD[i] + log(k_D) + (k_D - 1) * log(T[i]) - bD[i] * pow(T[i], k_D)
    logSD[i] <- - bD[i] * pow(T[i], k_D)
    
    # Censor C
    logfC[i] <- log_bC[i] + log(k_C) + (k_C - 1) * log(T[i]) - bC[i] * pow(T[i], k_C)
    logSC[i] <- - bC[i] * pow(T[i], k_C)
    
    
    # Observed-data log-likelihood:
    # d=1:            log f_D(T) + log S_C(T)
    # d=0 & !admin:   log f_C(T) + log S_D(T)          # random censoring at T < tau
    # d=0 & admin:    log S_C(tau) + log S_D(tau)      # adminstrative censoring at tau
    logL[i] <- d[i] * (logfD[i] + logSC[i]) +
      (1 - d[i]) * ( (1 - admin[i]) * (logfC[i] + logSD[i]) + admin[i] * ( - bC[i] * pow(tau, k_C)  - bD[i] * pow(tau, k_D) ))
    
    
    # zeros trick
    lambda[i] <- -logL[i] + Cconst
    zeros[i] ~ dpois(lambda[i])

  }
  
  # ----- Priors 
  # Event
  eta0   ~ dnorm(0, 0.1)
  eta_x1 ~ dnorm(0, 0.25)
  eta_x2 ~ dnorm(0, 0.25)
  eta_a  ~ dnorm(0, 0.25)
  eta_u  ~ dnorm(0, 0.25)
  
  # Censoring
  zeta0   ~ dnorm(0, 0.1)
  zeta_x1 ~ dnorm(0, 0.25)
  zeta_x2 ~ dnorm(0, 0.25)
  zeta_a  ~ dnorm(0, 0.25)
  zeta_u  ~ dnorm(0, 0.25) 
  
  # Shapes
  k_D ~ dgamma(9, 6)   # event shape
  k_C ~ dgamma(6, 5)  # censoring shape
  
  # Derived quantity;
  # Dirichlet weights for N observations
  
  for (i in 1:N) {
    alpha_w[i] <- 1.0
  }
  
  w[1:N] ~ ddirch(alpha_w[])
  # Survival at t0 under A=1 and A=0
  t0k <- pow(t0, k_D)
  
  for (i in 1:N) {
    # A=1
    muD1[i]    <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*1 + eta_u*U[i]
    log_bD1[i] <- -k_D * muD1[i]
    bD1[i]     <- exp(log_bD1[i])
    S1_ind[i]  <- exp(-bD1[i] * t0k)
    
    # A=0
    muD0[i]    <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*0 + eta_u*U[i]
    log_bD0[i] <- -k_D * muD0[i]
    bD0[i]     <- exp(log_bD0[i])
    S0_ind[i]  <- exp(-bD0[i] * t0k)
  }
  
  # Marginal survival and SPCE
  S1_marg <- inprod(w[], S1_ind[])
  S0_marg <- inprod(w[], S0_ind[])
  SPCE    <- S1_marg - S0_marg
}"
  
  X1 <- rbinom(N, size = 1, prob = 0.5)
  X2 <- rnorm(N, mean = 0, sd = 1)
  
  #unmeasured confounder
  U <- rbinom(N, size = 1, prob = 0.4)
  
  # ---- Treatment assignment ----
  linpred_A <- -0.4 + 0.2 * X1 + 0.3 * X2 + 0.6 * U
  pi_A      <- plogis(linpred_A)
  A         <- rbinom(N, size = 1, prob = pi_A)
  
  # ---- Event times (Weibull) ----
  #generate via inverse CDF
  #T~weibull(k, b), F(t)=1-exp(- b *t^k), where b = exp(-k* lp)
  #T=(-log(U)/b)^(1/k)
  mu1 <- 0.9 + 0.1 * X1 - 0.25 * X2 + 0.6 * U + 0.8 
  mu0 <- 0.9 + 0.1 * X1 - 0.25 * X2 + 0.6 * U
  b1 <- exp(-1.5 * mu1)
  b0 <- exp(-1.5 * mu0)
  
  D_a0 <- (-log(runif(N)) / b0)^(1 / 1.5)
  D_a1 <- (-log(runif(N)) / b1)^(1 / 1.5)
  
  # ---- informative censoring: ONE censoring time per subject ----
  bC0 <- exp(-1.2 * (1 + 0.4 * X1 + 0.3 * X2 + 0.6 * U))
  bC1 <- exp(-1.2 * (1 + 0.4 * X1 + 0.3 * X2 + 0.6 * U + 0.3))
  
  C_a0 <- (-log(runif(N)) / bC0)^(1 / 1.2) #A=0
  C_a1 <- (-log(runif(N)) / bC1)^(1 / 1.2) #A=1
  
  # ---- Observed time & event indicator ----
  Time <- ifelse(A == 1, pmin(10, C_a1, D_a1),
                 pmin(10, C_a0, D_a0))
  d <- ifelse(A == 1,
              as.numeric(D_a1 <= pmin(10, C_a1)),
              as.numeric(D_a0 <= pmin(10, C_a0)))
  
  
  # --- prepare JAGS data  ---
  tau  <- 10
  T <- pmax(Time, 1e-8) #time used in log-terms (strictly > 0)
  eps   <- 1e-10 
  admin <- as.integer(d == 0 & abs(Time - 10) < eps)
  
  # --- Zeros trick components ---
  zeros <- rep(0L, N)     # all zeros
  Cconst <- 1000           # positive constant to ensure Poisson parameter > 0
  
  jags_data_latU <- list(
    N = N,
    T = T,       
    d = d,
    admin = admin,
    A = A,
    X1 = X1,
    X2 = X2,
    zeros = zeros,
    Cconst = Cconst,
    tau = 10,
    t0 = t0
  )
  
  jags_data_naiveU <- list(
    N = N,
    T = T,
    d = d,
    admin = admin,
    A = A,
    X1 = X1,
    X2 = X2,
    zeros = zeros,
    Cconst = Cconst,
    tau = 10,
    t0 = t0
  )
  
  jags_data_obsU <- list(
    N = N,
    T = T,
    d = d,
    admin = admin,
    A = A,
    X1 = X1,
    X2 = X2,
    U  = U,
    zeros = zeros,
    Cconst = Cconst,
    tau = 10,
    t0 = t0
  )
  
  # --- monitor parameters  ---
  params_latU <- c("SPCE")
  params_naiveU <- c("SPCE")
  params_obsU <- c("SPCE")
  
  #---  initial values --- 
  inits_latU <- function() {list(
    eta0 = 0, eta_x1 = 0, eta_x2 = 0, eta_a = 0, eta_u= 0,
    zeta0 = 0, zeta_x1 = 0, zeta_x2 = 0, zeta_a = 0, zeta_u= 0,
    gamma0 = 0, gamma1 = 0, gamma2 = 0,
    alpha0 = 0, alpha1 = 0, alpha2 = 0, alpha3 = 0,
    k_D = 1, k_C = 1)}
  inits_naiveU <- function() {list(
    eta0 = 0, eta_x1 = 0, eta_x2 = 0, eta_a = 0, 
    zeta0 = 0, zeta_x1 = 0, zeta_x2 = 0, zeta_a = 0,
    k_D = 1, k_C = 1)}
  inits_obsU <- function() {list(
    eta0 = 0, eta_x1 = 0, eta_x2 = 0, eta_a = 0, eta_u= 0, 
    zeta0 = 0, zeta_x1 = 0, zeta_x2 = 0, zeta_a = 0, zeta_u= 0,
    k_D = 1, k_C = 1)}
  
  # --- fit JAGS once per replicate  ---
  sim1 <- jags(data = jags_data_latU, inits = inits_latU, parameters.to.save = params_latU,
               model.file  = textConnection(jags_model_latU),
               n.chains = 1, n.iter = nburn + niter, n.burnin = nburn, n.thin = 1, jags.seed = seed)
  
  sim2 <- jags(data = jags_data_naiveU, inits = inits_naiveU, parameters.to.save = params_naiveU,
               model.file  = textConnection(jags_model_naiveU),
               n.chains = 1, n.iter = nburn + niter, n.burnin = nburn, n.thin = 1, jags.seed = seed+1000)
  
  sim3 <- jags(data = jags_data_obsU, inits = inits_obsU, parameters.to.save = params_obsU,
               model.file  = textConnection(jags_model_obsU),
               n.chains = 1, n.iter = nburn + niter, n.burnin = nburn, n.thin = 1, jags.seed = seed+2000)
  
  #--------------------------------------------------
  post_lat <- as.matrix(as.mcmc(sim1))
  SPCE_draws_latU <- post_lat[, "SPCE"]
  SPCE_sd_latU    <- sd(SPCE_draws_latU)
  SPCE_mean_latU  <- mean(SPCE_draws_latU)
  SPCE_CI_latU    <- quantile(SPCE_draws_latU, c(0.025, 0.975), names = FALSE)
  
  post_naive <- as.matrix(as.mcmc(sim2))
  SPCE_draws_naiveU <- post_naive[, "SPCE"]
  SPCE_sd_naiveU    <- sd(SPCE_draws_naiveU)
  SPCE_mean_naiveU  <- mean(SPCE_draws_naiveU)
  SPCE_CI_naiveU   <- quantile(SPCE_draws_naiveU, c(0.025, 0.975), names = FALSE)
  
  post_obs <- as.matrix(as.mcmc(sim3))
  SPCE_draws_obsU <- post_obs[, "SPCE"]
  SPCE_sd_obsU    <- sd(SPCE_draws_obsU)
  SPCE_mean_obsU  <- mean(SPCE_draws_obsU)
  SPCE_CI_obsU   <- quantile(SPCE_draws_obsU, c(0.025, 0.975), names = FALSE)
  
  results <- list(
    spce_draws_latU = SPCE_draws_latU,
    spce_mean_latU  = SPCE_mean_latU,
    spce_sd_latU    = SPCE_sd_latU,
    spce_ci_latU    = SPCE_CI_latU,
    spce_draws_naiveU = SPCE_draws_naiveU,
    spce_mean_naiveU  = SPCE_mean_naiveU,
    spce_sd_naiveU    = SPCE_sd_naiveU,
    spce_ci_naiveU    = SPCE_CI_naiveU,
    spce_draws_obsU = SPCE_draws_obsU,
    spce_mean_obsU = SPCE_mean_obsU,
    spce_sd_obsU    = SPCE_sd_obsU,
    spce_ci_obsU    = SPCE_CI_obsU
  )
  
  return(results)
}

# Detect number of cores allocated by Slurm
nsim <- 1000
ncores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

# Safety cap 
ncores <- min(ncores, 120, nsim)

cat("Detected cores:", ncores, "\n")

cl <- makeCluster(ncores)
registerDoParallel(cl)

# ---- Define seeds handled by this node ----
# For now: one node runs seeds 1:ncores
seeds <- 1:nsim

# ---- Run simulations in parallel ----
results <- foreach(seed = seeds,
                   .packages = c("R2jags", "coda"),
                   .export = c("run_one_sim"),
                   .options.RNG = 123) %dorng% {
                     run_one_sim(seed)
                   }

stopCluster(cl)

# 6. Save results

out_dir <- file.path(Sys.getenv("SCRATCH"), "sensitivity", "out")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

outfile <- file.path(out_dir,
                     sprintf("info_BB_pf1_ss1000.rds", ncores))

saveRDS(results, outfile)

cat("Completed", length(results), "simulations\n")
cat("Saved to:", outfile, "\n")
