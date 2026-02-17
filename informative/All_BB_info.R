t_star <- Sys.time()
library(doParallel)
library(doRNG) 
library(R2jags)
library(coda)

#-----------------------------------------------------------------#
# Data generator function (AFT Weibull) - informative censoring.  #
#-----------------------------------------------------------------#
gen_data_sim <- function(n, tau, alpha_0, alpha_x1, alpha_x2, alpha_u, #treatment parameters
                         eta_intercept, eta_x1, eta_x2, eta_u, eta_a, #potential event time parameters
                         zeta_intercept, zeta_a, zeta_x1, zeta_x2, zeta_u, #potential censored time parameters
                         k, pU) { #lambda_C censoring rate
  # ---- Covariate ----
  X1 <- rbinom(n, size = 1, prob = 0.5)
  X2 <- rnorm(n, mean = 0, sd = 1)
  
  #unmeasured confounder
  U <- rbinom(n, size = 1, prob = pU)
  
  # ---- Treatment assignment ----
  linpred_A <- alpha_0 + alpha_x1 * X1 + alpha_x2 * X2 + alpha_u * U
  pi_A      <- plogis(linpred_A)
  A         <- rbinom(n, size = 1, prob = pi_A)
  
  # ---- Event times (Weibull) ----
  #generate via inverse CDF
  #T~weibull(k, b), F(t)=1-exp(- b *t^k), where b = exp(-k* lp)
  #T=(-log(U)/b)^(1/k)
  mu1 <- eta_intercept + eta_x1 * X1 + eta_x2 * X2 + eta_u * U + eta_a 
  mu0 <- eta_intercept + eta_x1 * X1 + eta_x2 * X2 + eta_u * U
  b1 <- exp(-k * mu1)
  b0 <- exp(-k * mu0)
  
  D_a0 <- (-log(runif(n)) / b0)^(1 / k)
  D_a1 <- (-log(runif(n)) / b1)^(1 / k)
  
  # ---- informative censoring: ONE censoring time per subject ----
  linpred_C0 <- zeta_intercept + zeta_x1 * X1 + zeta_x2 * X2 + zeta_u  * U 
  linpred_C1 <- zeta_intercept + zeta_x1 * X1 + zeta_x2 * X2 + zeta_u  * U + zeta_a 
  C_a0 <- -log(runif(n)) / (exp(linpred_C0)) #A=0
  C_a1 <- -log(runif(n)) / (exp(linpred_C1)) #A=1
  
  # ---- Observed time & event indicator ----
  Time <- ifelse(A == 1, pmin(tau, C_a1, D_a1),
                 pmin(tau, C_a0, D_a0))
  d <- ifelse(A == 1,
              as.numeric(D_a1 <= pmin(tau, C_a1)),
              as.numeric(D_a0 <= pmin(tau, C_a0)))
  
  data.frame(
    D_a0 = D_a0, D_a1 = D_a1,
    C_a0 = C_a0, C_a1 = C_a1,
    Time = Time, d = d,
    A = A, X1 = X1, X2 = X2, U = U
  )
}

run_one_rep <- function(r,
                        nburn = 5000, niter = 10000, 
                        N = 1000, t0 = 4) {
  
  
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
  eta_u  ~ dnorm(0, 0.25)
  
  # Censoring
  zeta0   ~ dnorm(0, 0.1)
  zeta_x1 ~ dnorm(0, 0.25)
  zeta_x2 ~ dnorm(0, 0.25)
  zeta_a  ~ dnorm(0, 0.25)
  zeta_u  ~ dnorm(0, 0.25)
  
  
  # U model
  gamma0 ~ dunif(-5, 5)
  gamma1 ~ dunif(-2, 2)
  gamma2 ~ dunif(-2, 2)
  
  # Treatment model
  alpha0 ~ dnorm(0, 0.1)
  alpha1 ~ dnorm(0, 0.25)
  alpha2 ~ dnorm(0, 0.25)
  alpha3 ~ dnorm(0, 0.25)
  
  # Shapes
  k_D ~ dgamma(9, 6)   # event shape
  k_C ~ dgamma(5.76, 4.8)  # censoring shape
  
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
  k_C ~ dgamma(5.76, 4.8)  # censoring shape
  
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
  k_C ~ dgamma(5.76, 4.8)  # censoring shape
  
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

  
  set.seed(1000 + r)
  
  dat <- gen_data_sim(
    n = N, tau = 10,
    alpha_0 = -0.4, alpha_x1 = 0.2, alpha_x2 = 0.3, alpha_u = 0.4,
    eta_intercept = 0.9, eta_a = 0.8, eta_x1 = 0.1, eta_x2 = -0.25, eta_u = 0.4,
    zeta_intercept = -2, zeta_a = -0.25, zeta_x1 = 0.15, zeta_x2 = -0.15, zeta_u = 0.25,
    k = 1.5, pU = 0.4
  )
  
  #mean(dat$D_a1>4)-mean(dat$D_a0>4),0.37963
  #censor rate, mean(dat$d == 0)
  
  # --- prepare JAGS data  ---
  tau  <- 10
  N <- nrow(dat)
  Time_obs <- dat$Time
  T <- pmax(Time_obs, 1e-8) #time used in log-terms (strictly > 0)
  d   <- as.integer(dat$d)
  eps   <- 1e-10 
  admin <- as.integer(d == 0 & abs(Time_obs - tau) < eps)
  A  <- as.integer(dat$A)
  X1 <- as.integer(dat$X1)
  X2 <- as.numeric(dat$X2)
  U <- as.numeric(dat$U)
  
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
    tau = tau,
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
    tau = tau,
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
    tau = tau,
    t0 = t0
  )
  
  # --- monitor parameters  ---
  params_latU <- c("SPCE")
  params_naiveU <- c("SPCE")
  params_obsU <- c("SPCE")
  
  # --- fit JAGS once per replicate  ---
  sim1 <- jags(data = jags_data_latU, parameters.to.save = params_latU, 
               model.file  = textConnection(jags_model_latU),
               n.chains = 1, n.iter = nburn + niter, n.burnin = nburn, n.thin = 1, jags.seed = 5000 + r)
  
  sim2 <- jags(data = jags_data_naiveU, parameters.to.save = params_naiveU, 
               model.file  = textConnection(jags_model_naiveU),
               n.chains = 1, n.iter = nburn + niter, n.burnin = nburn, n.thin = 1, jags.seed = 5000 + r)
  
  sim3 <- jags(data = jags_data_obsU, parameters.to.save = params_obsU, 
               model.file  = textConnection(jags_model_obsU),
               n.chains = 1, n.iter = nburn + niter, n.burnin = nburn, n.thin = 1, jags.seed = 5000 + r)
  
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
  
  list(
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
}

# -------- parallel across replicates r (no inner parallel) --------
R <- 10

n.cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl)

set.seed(123)  # master seed for %dorng%

results <- foreach(r = 1:R,
                   .packages = c("R2jags", "MCMCprecision"),
                   .export   = c("gen_data_sim",
                                 "run_one_rep")) %dorng% { run_one_rep(r) }

parallel::stopCluster(cl)

t_end <- Sys.time()
runtime_min <- as.numeric(difftime(t_end, t_star, units = "mins"))
runtime_min
#===============================
t_star <- Sys.time()
res1 <- run_one_rep(1)
t_end <- Sys.time()
runtime_min <- as.numeric(difftime(t_end, t_star, units = "mins"))
runtime_min
#res1$spce_mean_naiveU;res1$spce_mean_latU;res1$spce_mean_obsU
#res1$spce_mean_naiveU-true_spce_aft;res1$spce_mean_latU-true_spce_aft;res1$spce_mean_obsU-true_spce_aft
#=========================
library(dplyr)
library(tidyr)
true_spce_aft <- 0.37963 #alpha_u = 0.4, eta_u = 0.4, zeta_u = 0.25

comb_long <- do.call(bind_rows, lapply(results, function(x) {
  if (is.null(x)) return(NULL)
  tibble(
    method = c("naiveU", "latU", "obsU"),
    mean_est   = c(x$spce_mean_naiveU, x$spce_mean_latU, x$spce_mean_obsU),
    sd     = c(x$spce_sd_naiveU, x$spce_sd_latU, x$spce_sd_obsU),
    ci_lo  = c(x$spce_ci_naiveU[1], x$spce_ci_latU[1], x$spce_ci_obsU[1]),
    ci_hi  = c(x$spce_ci_naiveU[2], x$spce_ci_latU[2], x$spce_ci_obsU[2])
  )
}))

comb_long <- do.call(bind_rows, lapply(results, function(x) {
  if (is.null(x)) return(NULL)
  tibble(
    method = c("naiveU", "latU", "obsU"),
    mean_est   = c(x$spce_mean_naiveU, x$spce_mean_latU, x$spce_mean_obsU),
    sd     = c(x$spce_sd_naiveU, x$spce_sd_latU, x$spce_sd_obsU),
    ci_lo  = c(x$spce_ci_naiveU[1], x$spce_ci_latU[1], x$spce_ci_obsU[1]),
    ci_hi  = c(x$spce_ci_naiveU[2], x$spce_ci_latU[2], x$spce_ci_obsU[2])
  )
}))
#dat1<- subset(comb_long, method == "latU", select = mean_est)
#boxplot(dat1)
#---------------------------------------------------------#
# Compute coverage, bias, and summary stats all together. #
#---------------------------------------------------------#
summary_comb <- comb_long %>%
  mutate(
    method = factor(method, levels = c("naiveU", "latU", "obsU")), 
    bias   = mean_est - true_spce_aft,
    cover  = (true_spce_aft >= ci_lo & true_spce_aft <= ci_hi)
  ) %>%
  group_by(method) %>%
  summarise(
    coverage  = mean(cover),
    post_mean  = mean(mean_est),
    mean_bias = mean(bias),
    ASE   = mean(sd),
    ESE   = sd(mean_est),
    .groups = "drop"
  ) %>%
  arrange(method)

summary_comb


  