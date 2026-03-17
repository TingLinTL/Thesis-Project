dat <- gen_data_sim(
  n = N, tau = 10,
  alpha_0 = -0.4, alpha_x1 = 0.2, alpha_x2 = 0.3, alpha_u = 0.4,
  eta_intercept = 0.9, eta_a = 0.8, eta_x1 = 0.1, eta_x2 = -0.25, eta_u = 0.4,
  zeta_intercept = -2, zeta_a = -0.25, zeta_x1 = 0.15, zeta_x2 = -0.15, zeta_u = 0.25,
  k = 1.5, pU = 0.4
)

n <- 3000000
tau <- 10
X1 <- rbinom(n, size = 1, prob = 0.5)
X2 <- rnorm(n, mean = 0, sd = 1)

#unmeasured confounder
U <- rbinom(n, size = 1, prob = 0.4)

# ---- Treatment assignment ----
linpred_A <- -0.4 + 0.2 * X1 + 0.3 * X2 + 0.6 * U
pi_A      <- plogis(linpred_A)
A         <- rbinom(n, size = 1, prob = pi_A)

# ---- Event times (Weibull) ----
#generate via inverse CDF
#T~weibull(k, b), F(t)=1-exp(- b *t^k), where b = exp(-k* lp)
#T=(-log(U)/b)^(1/k)
mu1 <- 0.9 + 0.1 * X1 - 0.25 * X2 + 0.6 * U + 0.8 
mu0 <- 0.9 + 0.1 * X1 - 0.25 * X2 + 0.6 * U
b1 <- exp(-1.5 * mu1)
b0 <- exp(-1.5 * mu0)

D_a0 <- (-log(runif(n)) / b0)^(1 / 1.5)
D_a1 <- (-log(runif(n)) / b1)^(1 / 1.5)

# ---- informative censoring: ONE censoring time per subject ----
bC0 <- exp(-1.2 * (1+ 0.4 * X1 + 0.3 * X2 + 0.6 * U))
bC1 <- exp(-1.2 * (1+ 0.4 * X1 + 0.3 * X2 + 0.6 * U + 0.3))

C_a0 <- (-log(runif(n)) / bC0)^(1 / 1.2) #A=0
C_a1 <- (-log(runif(n)) / bC1)^(1 / 1.2) #A=1


# ---- Observed time & event indicator ----
Time <- ifelse(A == 1, pmin(tau, C_a1, D_a1),
               pmin(tau, C_a0, D_a0))
d <- ifelse(A == 1,
            as.numeric(D_a1 <= pmin(tau, C_a1)),
            as.numeric(D_a0 <= pmin(tau, C_a0)))

dat <- data.frame(
  D_a0 = D_a0, D_a1 = D_a1,
  C_a0 = C_a0, C_a1 = C_a1,
  Time = Time, d = d,
  A = A, X1 = X1, X2 = X2, U = U
)
#mean(dat$D_a1>4)-mean(dat$D_a0>4), 0.363073
#censor rate, mean(dat$d == 0) ~53%
#censoring intercept 5, ~ 13.5%
#censoring intercept 1.6, ~ 41%