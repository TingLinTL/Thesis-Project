t_star <- Sys.time()
library(R2jags)
library(coda)
library(doParallel)
library(doRNG)            # reproducible foreach RNG streams


#-----------------------------------------------------------------#
# Data generator function (AFT Weibull) - non-informative censoring.  #
#-----------------------------------------------------------------#
gen_data_sim <- function(n, tau, alpha_0, alpha_x1, alpha_x2, alpha_u, #treatment parameters
                         eta_intercept, eta_x1, eta_x2, eta_u, eta_a, #event time parameters
                         k, pU, lambda_C) { #lambda_C censoring rate
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
  
  # ---- Non-informative censoring: ONE censoring time per subject ----
  C <- rexp(n, rate = lambda_C)
  
  # ---- Observed time & event indicator ----
  Time <- ifelse(A == 1, pmin(D_a1, C, tau),
                 pmin(D_a0, C, tau))
  
  d <- ifelse(A == 1,
              as.integer(D_a1 <= pmin(C, tau)),
              as.integer(D_a0 <= pmin(C, tau)))
  data.frame(Time = Time, d = d, A = A, X1 = X1, X2 = X2, U = U,
             C = C, D_a0 = D_a0, D_a1 = D_a1)
}


run_one_rep <- function(r,
                        nburn = 5000, niter = 10000, 
                        N = 1000, M = 5000, t0 = 4) {
  #nburn:Burn-in iterations; niter:number of posterior draws after burn-in
  #N: sample size; M:inner MC size; t0: predicted time
  
  #------------------------------------------------------#
  #     JAGS models for all three approaches.            #
  #------------------------------------------------------#
  jags_model_latU <- "
  model {
    for (i in 1:N) {
      # Right-censoring via interval indicator
      is_cens[i] ~ dinterval(T[i], C[i,1])
      
      # Outcome model: Weibull distribution, dweib(shape = k, rate = b[i])
      T[i] ~ dweib(k_D, bD[i])
      
      # AFT linear predictor with latent U 
      # log(lambda_i) = eta0 + eta_x * X[i] + eta_a * A[i] 
      muD[i] <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*A[i] + eta_u*U[i]
      
      # rate parameter b[i] = exp(-k * mu[i])
      # log b[i] = -k * mu[i]
      log_bD[i] <- -k_D * muD[i]
      bD[i] <- exp(log_bD[i])
      
      # ----- Treatment model
      A[i] ~ dbern(pA[i])
      logit(pA[i]) <- alpha0 + alpha1*X1[i] + alpha2*X2[i] + alpha3*U[i]
      
      # ----- Latent U model
      logit(pU[i]) <- gamma0 + gamma1 * X1[i] + gamma2 * X2[i]
      U[i] ~ dbern(pU[i])
      
      # ----- Covariate X1 model; true X1 <- rbinom(n, size = 1, prob = 0.5)
      X1[i] ~ dbern(pX1)
      
      # ----- Covariate X2 model; true X2 <- rnorm(n, mean = 0, sd = 1)
      X2[i] ~ dnorm(muX2, tauX2) #mean and precision = 1/variance
    }
    
    # ----- Priors 
    # Event process
    eta0   ~ dnorm(0, 0.1) #precision=0.1, variance = 10
    eta_x1 ~ dnorm(0, 0.25) #precision=0.25, variance = 4
    eta_x2 ~ dnorm(0, 0.25)
    eta_a  ~ dnorm(0, 0.25)
    #eta_u  ~ dunif(-3, 3)
    eta_u  ~ dnorm(0, 0.25)
    
    # U model
    gamma0 ~ dunif(-3, 3)
    gamma1 ~ dunif(-2, 2)
    gamma2 ~ dunif(-2, 2)
    
    #X1 model
    pX1 ~ dbeta(1, 1) #uniform(0,1)
    #X2 model
    muX2 ~ dnorm(0, 0.1)
    tauX2 ~ dgamma(2, 2)
    
    # Treatment model
    alpha0 ~ dnorm(0, 0.1)
    alpha1 ~ dnorm(0, 0.25)
    alpha2 ~ dnorm(0, 0.25)
    #alpha3 ~ dunif(-3, 3)
    alpha3  ~ dnorm(0, 0.25)
    
    # Shapes
    k_D ~ dgamma(9, 6)   # event shape
    
    # Survival at t0 under A=1 and A=0
    t0k <- pow(t0, k_D)
    
    
    for (m in 1:M) {
      
      # --- draw new covariates from parametric models
      X1_new[m] ~ dbern(pX1)
      X2_new[m] ~ dnorm(muX2, tauX2)
      
      # --- draw latent U from its model
      logit(pU_new[m]) <- gamma0 + gamma1 * X1_new[m] + gamma2 * X2_new[m]
      U_new[m] ~ dbern(pU_new[m])
      
      # --- A = 1
      muD1_new[m] <- eta0 + eta_x1 * X1_new[m] + eta_x2 * X2_new[m] + eta_a  * 1 + eta_u  * U_new[m]
      log_bD1_new[m] <- -k_D * muD1_new[m]
      bD1_new[m]     <- exp(log_bD1_new[m])
      S1_new[m]      <- exp(-bD1_new[m] * t0k)
      
      # --- A = 0
      muD0_new[m] <- eta0 + eta_x1 * X1_new[m] + eta_x2 * X2_new[m] + eta_a  * 0 + eta_u  * U_new[m]
      
      log_bD0_new[m] <- -k_D * muD0_new[m]
      bD0_new[m]     <- exp(log_bD0_new[m])
      S0_new[m]      <- exp(-bD0_new[m] * t0k)
    }
    
    
    # Marginal survival and SPCE
    S1_marg <- mean(S1_new[])
    S0_marg <- mean(S0_new[])
    SPCE    <- S1_marg - S0_marg
  }"
  
  jags_model_naiveU <- "
  model {
    for (i in 1:N) {
      # Right-censoring via interval indicator
      is_cens[i] ~ dinterval(T[i], C[i,1])
      
      # Outcome model: Weibull distribution, dweib(shape = k, rate = b[i])
      T[i] ~ dweib(k_D, bD[i])
      
      # AFT linear predictor without latent U 
      # log(lambda_i) = eta0 + eta_x * X[i] + eta_a * A[i] 
      muD[i] <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*A[i]
      
      # rate parameter b[i] = exp(-k * mu[i])
      # log b[i] = -k * mu[i]
      log_bD[i] <- -k_D * muD[i]
      bD[i] <- exp(log_bD[i])
      
      # ----- Covariate X1 model; true X1 <- rbinom(n, size = 1, prob = 0.5)
      X1[i] ~ dbern(pX1)
      
      # ----- Covariate X2 model; true X2 <- rnorm(n, mean = 0, sd = 1)
      X2[i] ~ dnorm(muX2, tauX2) #mean and precision = 1/variance
    }
    
    # ----- Priors 
    # Event process
    eta0   ~ dnorm(0, 0.1)
    eta_x1 ~ dnorm(0, 0.25)
    eta_x2 ~ dnorm(0, 0.25)
    eta_a  ~ dnorm(0, 0.25)
    
    #X1 model
    pX1 ~ dbeta(1, 1)
    #X2 model
    muX2 ~ dnorm(0, 0.1)
    tauX2 ~ dgamma(2, 2)

    # Shapes
    k_D ~ dgamma(9, 6)   # event shape
    
    # Survival at t0 under A=1 and A=0
    t0k <- pow(t0, k_D)
    
    for (m in 1:M) {
      
      # --- draw new covariates from parametric models
      X1_new[m] ~ dbern(pX1)
      X2_new[m] ~ dnorm(muX2, tauX2)
      
      # --- A = 1
      muD1_new[m] <- eta0 + eta_x1 * X1_new[m] + eta_x2 * X2_new[m] + eta_a  * 1
      log_bD1_new[m] <- -k_D * muD1_new[m]
      bD1_new[m]     <- exp(log_bD1_new[m])
      S1_new[m]      <- exp(-bD1_new[m] * t0k)
      
      # --- A = 0
      muD0_new[m] <- eta0 + eta_x1 * X1_new[m] + eta_x2 * X2_new[m] + eta_a  * 0 
      
      log_bD0_new[m] <- -k_D * muD0_new[m]
      bD0_new[m]     <- exp(log_bD0_new[m])
      S0_new[m]      <- exp(-bD0_new[m] * t0k)
    }
    
    
    # Marginal survival and SPCE
    S1_marg <- mean(S1_new[])
    S0_marg <- mean(S0_new[])
    SPCE    <- S1_marg - S0_marg
  }"
  
  jags_model_obsU <- "
  model {
    for (i in 1:N) {
      # Right-censoring via interval indicator
      is_cens[i] ~ dinterval(T[i], C[i,1])
      
      # Outcome model: Weibull distribution, dweib(shape = k, rate = b[i])
      T[i] ~ dweib(k_D, bD[i])
      
      # AFT linear predictor without latent U 
      # log(lambda_i) = eta0 + eta_x * X[i] + eta_a * A[i] 
      muD[i] <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*A[i] + eta_u *U[i]
      
      # rate parameter b[i] = exp(-k * mu[i])
      # log b[i] = -k * mu[i]
      log_bD[i] <- -k_D * muD[i]
      bD[i] <- exp(log_bD[i])
      
      # ---- if U is observed
      U[i] ~ dbern(pU)
      
      # ----- Covariate X1 model; true X1 <- rbinom(n, size = 1, prob = 0.5)
      X1[i] ~ dbern(pX1)
      
      # ----- Covariate X2 model; true X2 <- rnorm(n, mean = 0, sd = 1)
      X2[i] ~ dnorm(muX2, tauX2) #mean and precision = 1/variance
    }
    
    # ----- Priors 
    # Event process
    eta0   ~ dnorm(0, 0.1)
    eta_x1 ~ dnorm(0, 0.25)
    eta_x2 ~ dnorm(0, 0.25)
    eta_a  ~ dnorm(0, 0.25)
    eta_u  ~ dnorm(0, 0.25)
    
    #X1 model
    pX1 ~ dbeta(1, 1)
    #X2 model
    muX2 ~ dnorm(0, 0.1)
    tauX2 ~ dgamma(2, 2)
    #U model
    pU ~ dbeta(1, 1)

    # Shapes
    k_D ~ dgamma(9, 6)   # event shape
    
    # Survival at t0 under A=1 and A=0
    t0k <- pow(t0, k_D)
    
    for (m in 1:M) {
      
      # --- draw new covariates from parametric models
      X1_new[m] ~ dbern(pX1)
      X2_new[m] ~ dnorm(muX2, tauX2)
      U_new[m] ~ dbern(pU)
      
      # --- A = 1
      muD1_new[m] <- eta0 + eta_x1 * X1_new[m] + eta_x2 * X2_new[m] + eta_a  * 1 + eta_u * U_new[m]
      log_bD1_new[m] <- -k_D * muD1_new[m]
      bD1_new[m]     <- exp(log_bD1_new[m])
      S1_new[m]      <- exp(-bD1_new[m] * t0k)
      
      # --- A = 0
      muD0_new[m] <- eta0 + eta_x1 * X1_new[m] + eta_x2 * X2_new[m] + eta_a  * 0 + eta_u * U_new[m]
      
      log_bD0_new[m] <- -k_D * muD0_new[m]
      bD0_new[m]     <- exp(log_bD0_new[m])
      S0_new[m]      <- exp(-bD0_new[m] * t0k)
    }
    
    
    # Marginal survival and SPCE
    S1_marg <- mean(S1_new[])
    S0_marg <- mean(S0_new[])
    SPCE    <- S1_marg - S0_marg
  }"
  # --- simulate one dataset  ---
  dat <- gen_data_sim(
    n = N, tau =10,
    alpha_0 = -0.4, alpha_x1 = 0.2, alpha_x2 = 0.3, alpha_u = 0.6,
    eta_intercept = 0.9, eta_a = 0.8, eta_x1 = 0.1, eta_x2 = -0.25, eta_u = 0.6,
    k = 1.5, pU = 0.4, lambda_C = 0.001
  )
  
  #0.1 35% censoring #0.001 13%
  #censor rate, mean(dat$d == 0)
  #predict time t0 = 4
  #mean(dat$D_a1>4)-mean(dat$D_a0>4), true = 0.3797853
  # --- prepare JAGS data  ---
  tau  <- 10
  N <- nrow(dat)          #number of observations
  cens <- matrix(c(dat$Time, rep(NA, length(dat$Time))),
                 nrow = length(dat$Time), ncol = 2)
  dat$Time[dat$d == 0] <- NA # Censored
  is.censored <- as.numeric(is.na(dat$Time))
  
  jags_data_latU <- list(
    N       = nrow(dat),
    T       = dat$Time,                # observed time
    is_cens = is.censored,                 # 1 = censored, 0 = event
    C       = cens,        # effective censoring time
    A       = dat$A,
    X1      = dat$X1,
    X2      = dat$X2,
    t0      = t0, # time point for SPCE
    M       = M# Monte Carlo size for parametric g-computation
  )
  
  jags_data_naiveU <- list(
    N       = nrow(dat),
    T       = dat$Time,                # observed time
    is_cens = is.censored,                 # 1 = censored, 0 = event
    C       = cens,        # effective censoring time
    A       = dat$A,
    X1      = dat$X1,
    X2      = dat$X2,
    t0      = t0, # time point for SPCE
    M       = M# Monte Carlo size for parametric g-computation
  )
  
  jags_data_obsU <- list(
    N       = nrow(dat),
    T       = dat$Time,                # observed time
    is_cens = is.censored,                 # 1 = censored, 0 = event
    C       = cens,        # effective censoring time
    A       = dat$A,
    X1      = dat$X1,
    X2      = dat$X2,
    U       = dat$U,
    t0      = t0, # time point for SPCE
    M       = M# Monte Carlo size for parametric g-computation
  )
  
  
  # --- monitor parameters  ---
  #params_latU <- c("eta0","eta_x1","eta_x2","eta_a","eta_u","pX1","muX2","tauX2","SPCE") 
  params_latU <- c("SPCE")
  params_naiveU <- c("SPCE")
  params_obsU <- c("SPCE")
  # --- fit JAGS once per replicate  ---
  #sim1 - latentU approach
  #sim1 <- jags.model(textConnection(jags_model_latU), data = jags_data_latU, n.chains = 1)
  #update(sim1, nburn)
  #post1 <- coda.samples(sim1, params_latU, n.iter = niter, thin = 1)
  
  sim1 <- jags(data = jags_data_latU, parameters.to.save = params_latU, 
               model.file  = textConnection(jags_model_latU),
               n.chains = 1, n.iter = nburn + niter, n.burnin = nburn, n.thin = 1)
  
  sim2 <- jags(data = jags_data_naiveU, parameters.to.save = params_naiveU, 
               model.file  = textConnection(jags_model_naiveU),
               n.chains = 1, n.iter = nburn + niter, n.burnin = nburn, n.thin = 1)
  
  sim3 <- jags(data = jags_data_obsU, parameters.to.save = params_obsU, 
               model.file  = textConnection(jags_model_obsU),
               n.chains = 1, n.iter = nburn + niter, n.burnin = nburn, n.thin = 1)
  
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
true_spce_aft <- 0.3797853 #0.4
true_spce_aft <- 0.3442317 #0.8
true_spce_aft <- 0.3638207 #0.6


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

