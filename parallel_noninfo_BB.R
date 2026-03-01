t_star <- Sys.time()
library(R2jags)
library(coda)
library(doParallel)
library(doRNG)            # reproducible foreach RNG streams


#parallel in Trillium
run_one_sim <- function(seed,nburn = 5000, niter = 10000,
                        N = 2000, t0 = 4) {
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

    }

    # ----- Priors 
    # Event process
    eta0   ~ dnorm(0, 0.1) #precision=0.1, variance = 10
    eta_x1 ~ dnorm(0, 0.25) #precision=0.25, variance = 4
    eta_x2 ~ dnorm(0, 0.25)
    eta_a  ~ dnorm(0, 0.25)
    #informative prior, center at true value
    #if OR is from (log(1/3),log(3)), then the mean = 0, sd = 0.5605
    eta_u ~ dnorm(0, 3.183)

    # U model
    #gamma0 ~ dnorm(-0.4055, 84)
    #gamma1 ~ dnorm(0, 15.36)
    #gamma2 ~ dnorm(0, 15.36)
    gamma0 ~ dnorm(0, 2.25)
    gamma1 ~ dnorm(0,1)
    gamma2 ~ dnorm(0,1)

    # Treatment model
    alpha0 ~ dnorm(0, 0.1)
    alpha1 ~ dnorm(0, 0.25)
    alpha2 ~ dnorm(0, 0.25)
    #alpha3 ~ dunif(-2, 2)
    #if OR is from (log(1/3),log(3)), then the mean = 0, sd = 0.5605
    alpha3 ~ dnorm(0, 3.183)

    # Shapes
    k_D ~ dgamma(9, 6)   # event shape

    # Survival at t0 under A=1 and A=0
    t0k <- pow(t0, k_D)
    
     # Derived quantity;
    # Dirichlet weights for N observations

    for (i in 1:N) {
     alpha_w[i] <- 1.0
     }

    w[1:N] ~ ddirch(alpha_w[])
    # Survival at t0 under A=1 and A=0

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
      # Right-censoring via interval indicator
      is_cens[i] ~ dinterval(T[i], C[i,1])

      # Outcome model: Weibull distribution, dweib(shape = k, rate = b[i])
      T[i] ~ dweib(k_D, bD[i])

      # AFT linear predictor without latent U 
      # log(lambda_i) = eta0 + eta_x1 * X1[i] + eta_x2 * X2[i] + eta_a * A[i] 
      muD[i] <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*A[i]

      # rate parameter b[i] = exp(-k * mu[i])
      # log b[i] = -k * mu[i]
      log_bD[i] <- -k_D * muD[i]
      bD[i] <- exp(log_bD[i])

    }

    # ----- Priors 
    # Event process
    eta0   ~ dnorm(0, 0.1)
    eta_x1 ~ dnorm(0, 0.25)
    eta_x2 ~ dnorm(0, 0.25)
    eta_a  ~ dnorm(0, 0.25)


    # Shapes
    k_D ~ dgamma(9, 6)   # event shape

    # Survival at t0 under A=1 and A=0
    t0k <- pow(t0, k_D)

    # Derived quantity;
    # Dirichlet weights for N observations

    for (i in 1:N) {
     alpha_w[i] <- 1.0
     }

    w[1:N] ~ ddirch(alpha_w[])
    # Survival at t0 under A=1 and A=0

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
      # Right-censoring via interval indicator
      is_cens[i] ~ dinterval(T[i], C[i,1])

      # Outcome model: Weibull distribution, dweib(shape = k, rate = b[i])
      T[i] ~ dweib(k_D, bD[i])

      # AFT linear predictor without latent U 
      # log(lambda_i) = eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a * A[i] 
      muD[i] <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*A[i] + eta_u *U[i]

      # rate parameter b[i] = exp(-k * mu[i])
      # log b[i] = -k * mu[i]
      log_bD[i] <- -k_D * muD[i]
      bD[i] <- exp(log_bD[i])

    }

    # ----- Priors 
    # Event process
    eta0   ~ dnorm(0, 0.1)
    eta_x1 ~ dnorm(0, 0.25)
    eta_x2 ~ dnorm(0, 0.25)
    eta_a  ~ dnorm(0, 0.25)
    #eta_u  ~ dnorm(0.6, 1/0.0225)
    eta_u ~ dnorm(0, 3.183)

    # Shapes
    k_D ~ dgamma(9, 6)   # event shape

    # Survival at t0 under A=1 and A=0
    t0k <- pow(t0, k_D)

    # Derived quantity;
    # Dirichlet weights for N observations

    for (i in 1:N) {
     alpha_w[i] <- 1.0
     }

    w[1:N] ~ ddirch(alpha_w[])
    # Survival at t0 under A=1 and A=0

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
  
  #Generating data
  # ---- Covariate ----
  X1 <- rbinom(N, size = 1, prob = 0.5)
  X2 <- rnorm(N, mean = 0, sd = 1)
  #unmeasured confounder
  U <- rbinom(N, size = 1, prob = 0.4)
  # ---- Treatment assignment ----
  pi_A      <- plogis(-0.4 + 0.2 * X1 + 0.3 * X2 + 0.6 * U)
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
  # ---- Non-informative censoring: ONE censoring time per subject ----
  C <- rexp(N, rate = 0.1)
  
  # ---- Observed time & event indicator ----
  Time <- ifelse(A == 1, pmin(D_a1, C, 10),
                 pmin(D_a0, C, 10))
  
  d <- ifelse(A == 1,
              as.integer(D_a1 <= pmin(C, 10)),
              as.integer(D_a0 <= pmin(C, 10)))
  cens <- matrix(c(Time, rep(NA, length(Time))),
                 nrow = length(Time), ncol = 2)
  Time[d == 0] <- NA # Censored
  is.censored <- as.numeric(is.na(Time))
  
  #prepare JAGS data
  jags_data_latU <- list(
    N       = N,
    T       = Time,                # observed time
    is_cens = is.censored,                 # 1 = censored, 0 = event
    C       = cens,        # effective censoring time
    A       = A,
    X1      = X1,
    X2      = X2,
    t0      = t0 # time point for SPCE
  )
  
  jags_data_naiveU <- list(
    N       = N,
    T       = Time,                # observed time
    is_cens = is.censored,                 # 1 = censored, 0 = event
    C       = cens,        # effective censoring time
    A       = A,
    X1      = X1,
    X2      = X2,
    t0      = t0# time point for SPCE
  )
  
  jags_data_obsU <- list(
    N       = N,
    T       = Time,                # observed time
    is_cens = is.censored,                 # 1 = censored, 0 = event
    C       = cens,        # effective censoring time
    A       = A,
    X1      = X1,
    X2      = X2,
    U       = U,
    t0      = t0          # time point for SPCE
  )
  
  # --- monitor parameters  ---
  #params_latU <- c("eta0","eta_x1","eta_x2","eta_a","eta_u","pX1","muX2","tauX2","SPCE") 
  params_latU <- c("SPCE")
  params_naiveU <- c("SPCE")
  params_obsU <- c("SPCE")
  
  #---  initial values --- 
  inits_latU <- function() {list(
    eta0 = 0, eta_x1 = 0, eta_x2 = 0, eta_a = 0, eta_u= 0,
    gamma0 = 0, gamma1 = 0, gamma2 = 0,
    alpha0 = 0, alpha1 = 0, alpha2 = 0, alpha3 = 0,
    k_D = 1)}
  inits_naiveU <- function() {list(
    eta0 = 0, eta_x1 = 0, eta_x2 = 0, eta_a = 0, k_D = 1)}
  inits_obsU <- function() {list(
    eta0 = 0, eta_x1 = 0, eta_x2 = 0, eta_a = 0, eta_u= 0, k_D = 1)}
    
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
t_end <- Sys.time()
run_one_sim(1)
