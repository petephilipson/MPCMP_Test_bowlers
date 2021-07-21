# MPCMP model using direct polynomial solving via polysolve
# This has separate blocks for runs and opposition (both splines)
# Separate update for game-specifics
# Additional update for toss x first innings interaction
#library(mvnfast) # For multivariate normal RNG
#library(splines)
# Can also load covariance proposal matrix (from mpcmp package fit for instance)
MPCMPpoly_toss_fit <- function(
  data, niter, thin, burnin, alphapriors = NULL, nupriors = NULL,
  xipriors = NULL, phipriors = NULL, zetapriors = NULL, thetapriors = NULL, 
  gammapriors = NULL, inits = NULL, OppMat, Sigma1, runsknots = 4)
{
  # Set up design matrix
  runs_knots <- as.numeric(c(quantile(log(data$Runs), probs = seq(0, 1, length = runsknots + 2))))
  BRuns <- bs(log(data$Runs), knots = runs_knots[2:(length(runs_knots) - 1)], degree = 3, intercept = TRUE)
  pruns <- ncol(BRuns)
  pplayers <- nlevels(data$ID) 
  XOpp <- as.matrix(OppMat)
  popps <- ncol(XOpp)
  OppInd <- XOpp != 0
  X1 <- model.matrix(~ BRuns, data = data)
  X <- X1[, -1]
  # Set up result matrix
  npara <- ncol(X) + 2*pplayers + popps + 9
  matres <- matrix(0, ncol = 1 + npara, nrow = niter)
  matres[, 1] <- 1:niter
  # Extract data summaries
  summ <- summarydata(data)
  # Prior values
  if(is.null(alphapriors)){alphapriors <- rep(c(0, 0.5), pruns)}
  if(is.null(nupriors)){nupriors <- rep(c(0, 0.5*log(3)), pplayers)} 
  if(is.null(thetapriors)){thetapriors <- rep(c(0, 0.25), pplayers)}
  if(is.null(xipriors)){xipriors <- c(rep(c(0, 0.5*log(2)), 4))}
  if(is.null(phipriors)){phipriors <- c(rep(c(0, 0.5*log(2)), 2))}
  if(is.null(gammapriors)){gammapriors <- c(rep(c(0, 0.5*log(2)), 2))}
  if(is.null(zetapriors)){zetapriors <- rep(c(0, 0.5), popps)}
  beta_m <- c(alphapriors[seq(1, length(alphapriors), 2)])
  beta_s <- c(alphapriors[seq(2, length(alphapriors), 2)])
  nu_m <- nupriors[seq(1, length(nupriors), 2)]
  nu_s <- nupriors[seq(2, length(nupriors), 2)]
  theta_m <- thetapriors[seq(1, length(thetapriors), 2)]
  theta_s <- thetapriors[seq(2, length(thetapriors), 2)]
  zeta_m <- zetapriors[seq(1, length(zetapriors), 2)]
  zeta_s <- zetapriors[seq(2, length(zetapriors), 2)]
  xi_m <- xipriors[seq(1, length(xipriors), 2)]
  xi_s <- xipriors[seq(2, length(xipriors), 2)]
  phi_m <- phipriors[seq(1, length(phipriors), 2)]
  phi_s <- phipriors[seq(2, length(phipriors), 2)]
  gamma_m <- gammapriors[seq(1, length(gammapriors), 2)]
  gamma_s <- gammapriors[seq(2, length(gammapriors), 2)]
  paras <- initialise(data, X, inits, alphapriors, nupriors, zetapriors, 
                      thetapriors, xipriors, phipriors, gammapriors)
  res <- NULL
  for (i in 1 : burnin){
    mupara <- muUpdate(data, X, XOpp, beta_m, beta_s, summ$unique, summ$lfac, 
                       paras$beta, paras$nu, paras$theta, paras$zeta, paras$xi, 
                       paras$phi, paras$gamma, sig_prop = Sigma1)
    nupara <- nuUpdate(data, X, XOpp, nu_m, nu_s, pplayers, summ$unique, summ$lfac, summ$lfacWkts, 
                       mupara$beta, paras$nu, paras$theta, paras$zeta, paras$xi, 
                       paras$phi, paras$gamma)
    thetapara <- thetaUpdate(data, X, XOpp, theta_m, theta_s, pplayers, summ$unique, 
                             summ$lfac, mupara$beta, nupara$nu, paras$theta, paras$zeta, 
                             paras$xi, paras$phi, paras$gamma)
    zetapara <- zetaUpdate(data, X, XOpp, OppInd, zeta_m, zeta_s, popps, summ$unique, 
                           summ$lfac, mupara$beta, nupara$nu, thetapara$theta, 
                           paras$zeta, paras$xi, paras$phi, paras$gamma)
    xipara <- xiUpdate(data, X, XOpp, xi_m, xi_s, summ$unique, summ$lfac, 
                       mupara$beta, nupara$nu, thetapara$theta, zetapara$zeta, 
                       paras$xi, paras$phi, paras$gamma)
    phipara <- phiUpdate(data, X, XOpp, phi_m, phi_s, summ$unique, summ$lfac, 
                         mupara$beta, nupara$nu, thetapara$theta, zetapara$zeta, 
                         xipara$xi, paras$phi, paras$gamma)
    gammapara <- gammaUpdate(data, X, XOpp, gamma_m, gamma_s, summ$unique, summ$lfac, 
                             mupara$beta, nupara$nu, thetapara$theta, zetapara$zeta, 
                             xipara$xi, phipara$phi, paras$gamma)
    paras <- list(beta = mupara$beta, nu = nupara$nu, theta = thetapara$theta, zeta = zetapara$zeta, 
                  xi = xipara$xi, phi = phipara$phi, gamma = gammapara$gamma)
  }
  for (j in 1 : niter){
    for (k in 1 : thin){
      mupara <- muUpdate(data, X, XOpp, beta_m, beta_s, summ$unique, summ$lfac, 
                         paras$beta, paras$nu, paras$theta, paras$zeta, paras$xi, 
                         paras$phi, paras$gamma, sig_prop = Sigma1)
      nupara <- nuUpdate(data, X, XOpp, nu_m, nu_s, pplayers, summ$unique, summ$lfac, summ$lfacWkts, 
                         mupara$beta, paras$nu, paras$theta, paras$zeta, paras$xi, 
                         paras$phi, paras$gamma)
      thetapara <- thetaUpdate(data, X, XOpp, theta_m, theta_s, pplayers, summ$unique, 
                               summ$lfac, mupara$beta, nupara$nu, paras$theta, paras$zeta, 
                               paras$xi, paras$phi, paras$gamma)
      zetapara <- zetaUpdate(data, X, XOpp, OppInd, zeta_m, zeta_s, popps, summ$unique, 
                             summ$lfac, mupara$beta, nupara$nu, thetapara$theta, 
                             paras$zeta, paras$xi, paras$phi, paras$gamma)
      xipara <- xiUpdate(data, X, XOpp, xi_m, xi_s, summ$unique, summ$lfac, 
                         mupara$beta, nupara$nu, thetapara$theta, zetapara$zeta, 
                         paras$xi, paras$phi, paras$gamma)
      phipara <- phiUpdate(data, X, XOpp, phi_m, phi_s, summ$unique, summ$lfac, 
                           mupara$beta, nupara$nu, thetapara$theta, zetapara$zeta, 
                           xipara$xi, paras$phi, paras$gamma)
      gammapara <- gammaUpdate(data, X, XOpp, gamma_m, gamma_s, summ$unique, summ$lfac, 
                               mupara$beta, nupara$nu, thetapara$theta, zetapara$zeta, 
                               xipara$xi, phipara$phi, paras$gamma)
      paras <- list(beta = mupara$beta, nu = nupara$nu, theta = thetapara$theta, zeta = zetapara$zeta, 
                    xi = xipara$xi, phi = phipara$phi, gamma = gammapara$gamma)
    }
    post_ll <- llhood(data, X, XOpp, summ$unique, summ$lfac, summ$lfacWkts, 
                      paras$beta, paras$nu, paras$theta, paras$zeta, 
                      paras$xi, paras$phi, paras$gamma)
    matres[j, -1] <- c(paras$beta, paras$nu, paras$theta, paras$zeta, paras$xi, paras$phi,
                       paras$gamma, post_ll$ll)
  }
  list(res = matres, dataset = data)
}

# Some data summaries needed in updates
summarydata <- function(data)
{
  N <- nrow(data)
  unique <- 0:10
  lfac <- lfactorial(unique)
  lfacWkts <- lfactorial(data$Wkts)
  list(N = N, unique = unique, lfac = lfac, lfacWkts = lfacWkts)
}

initialise <- function(data, X, inits, alphapriors, nupriors, zetapriors, 
                       thetapriors, xipriors, phipriors, gammapriors)
{
  # Initialise using supplied starting values or sampling from priors
  if(is.null(inits)){
    alpha <- rnorm(length(alphapriors)/2, alphapriors[seq(1, length(alphapriors), 2)], 
                   alphapriors[seq(2, length(alphapriors), 2)])
    nu <- rnorm(length(nupriors)/2, nupriors[seq(1, length(nupriors), 2)], 
                nupriors[seq(2, length(nupriors), 2)])
    xi <- rnorm(4, xipriors[c(1, 3, 5, 7)], xipriors[c(2, 4, 6, 8)])
    phi <- rnorm(2, phipriors[c(1, 3)], phipriors[c(2, 4)])
    zeta <- rnorm(length(zetapriors)/2, zetapriors[seq(1, length(zetapriors), 2)], 
                  zetapriors[seq(2, length(zetapriors), 2)])
    theta <- rnorm(length(thetapriors)/2, thetapriors[seq(1, length(thetapriors), 2)], 
                   thetapriors[seq(2, length(thetapriors), 2)])
    gamma <- rnorm(2, gammapriors[c(1, 3)], gammapriors[c(2, 4)])
  }
  else {
    alpha <- inits$alpha0
    nu <- inits$nu0
    xi <- c(0, inits$xi0[1:3])
    phi <- c(0, inits$xi0[4])
    zeta <- c(inits$zeta0)
    theta <- c(inits$theta0)
    gamma <- c(0, inits$xi0[5])
  }
  list(beta = alpha, nu = nu, theta = theta, zeta = zeta, xi = xi, phi = phi,
       gamma = gamma)
}

# Function to solve polynomial
findlam <- function(mu, nu){
  lam <- polyroot(((0:10) - mu)/factorial(0:10)^nu)
  lam <- Re(lam[Re(lam) >= 0 & round(Im(lam), 6) ==  0])
  if(length(lam) == 0){
    print(mu)
    print(nu)
  }
  if(lam < 0){
    print(mu)
    print(nu)
  }
  return (lam)
}

# Update for runs spline and game-specific covariates
muUpdate <- function(data, X, XOpp, beta_m ,beta_s, unique, lfac, 
                     beta, nu, theta, zeta, xi, phi, gamma, sig_prop){
  mu <- with(data, exp(X%*%beta + theta[ID] + XOpp%*%zeta +
                         xi[Inns] + phi[HA] + gamma[TossInns1]))
  if (max(mu) >= 10){mu <- pmin(9.9, mu)}
  if (min(mu) < 0.001){mu <- pmax(0.001, mu)}
  nu_lp <- exp(nu[data$ID])
  if (max(nu_lp >= 10)){nu_lp <- pmin(9.9, nu_lp)}
  if (min(nu_lp) < 0.001){nu_lp <- pmax(0.001, nu_lp)}
  #log_lambda <- log(lam_grid_poly_10K[cbind(round(mu*1000), round(nu_lp*1000))])
  log_lambda <- log(mapply(findlam, mu, nu_lp))
  # Proposal
  beta_p <- as.vector(rmvn(1, beta, 0.2*sig_prop)) 
  mu_p <- with(data, exp(X%*%beta_p + theta[ID] + XOpp%*%zeta + 
                           xi[Inns] + phi[HA] + gamma[TossInns1]))
  if (max(mu_p) >= 10){mu_p <- pmin(9.9, mu_p)}
  if (min(mu_p) < 0.001){mu_p <- pmax(0.001, mu_p)}
  #log_lambda_p <- log(lam_grid_poly_10K[cbind(round(mu_p*1000), round(nu_lp*1000))])
  log_lambda_p <- log(mapply(findlam, mu_p, nu_lp))
  lG <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(nu_lp, lfac)))) 
  lG_p <- log(rowSums(exp(tcrossprod(log_lambda_p, unique) - tcrossprod(nu_lp, lfac)))) 
  # Block update
  denrat <- sum(-0.5*((beta_p - beta_m)^2 - (beta - beta_m)^2)/beta_s^2)
  denrat <- denrat + sum(data$Wkts*(log_lambda_p - log_lambda) - (lG_p - lG))
  laccept <- min(0, denrat)
  accept <- (log(runif(1)) < laccept)
  #print(accept)
  if (accept){beta <- beta_p}
  list(beta = beta)
} 

# Dispersion parameters update (log-scale)
nuUpdate <- function(data, X, XOpp, nu_m, nu_s, pplayers, unique, lfac, 
                     lfacWkts, beta, nu, theta, zeta, xi, phi, gamma, rwsd){
  # One nu parameter for each player
  mu <- with(data, exp(X%*%beta + theta[ID] + XOpp%*%zeta + 
                         xi[Inns] + phi[HA] + gamma[TossInns1]))
  if (max(mu) >= 10){mu <- pmin(9.9, mu)}
  if (min(mu) < 0.001){mu <- pmax(0.001, mu)}
  nu_lp <- exp(nu[data$ID])
  if (max(nu_lp >= 10)){nu_lp <- pmin(9.9, nu_lp)}
  if (min(nu_lp) < 0.001){nu_lp <- pmax(0.001, nu_lp)}
  log_lambda <- log(mapply(findlam, mu, nu_lp))
  # Proposal 
  nu_p <- rnorm(pplayers, nu, 0.75)
  nu_p_lp <- exp(nu_p[data$ID])
  if (max(nu_p_lp >= 10)){nu_p_lp <- pmin(9.9, nu_p_lp)}
  if (min(nu_p_lp) < 0.001){nu_p_lp <- pmax(0.001, nu_p_lp)}
  log_lambda_p <- log(mapply(findlam, mu, nu_p_lp))
  # Calculate ln(G) terms
  lG <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(nu_lp, lfac)))) 
  lG_p <- log(rowSums(exp(tcrossprod(log_lambda_p, unique) - tcrossprod(nu_p_lp, lfac)))) 
  # Component-wise update 
  denrat <- with(data, tapply(Wkts*(log_lambda_p - log_lambda) + lG - lG_p + (exp(nu[ID]) - exp(nu_p[ID]))*lfacWkts, ID, sum))
  denrat <- denrat -0.5*((nu_p - nu_m)^2 - (nu - nu_m)^2)/nu_s^2
  laccept <- pmin(0, denrat)
  accept <- (log(runif(pplayers)) < laccept)
  nu[accept] <- nu_p[accept]
  list(nu = nu)	
}

# Player block
thetaUpdate <- function(data, X, XOpp, theta_m, theta_s, pplayers, unique, lfac, 
                        beta, nu, theta, zeta, xi, phi, gamma){
  mu <- with(data, exp(X%*%beta + theta[ID] + XOpp%*%zeta + 
                         xi[Inns] + phi[HA] + gamma[TossInns1]))
  if (max(mu) >= 10){mu <- pmin(9.9, mu)}
  if (min(mu) < 0.001){mu <- pmax(0.001, mu)}
  nu_lp <- exp(nu[data$ID])
  if (max(nu_lp >= 10)){nu_lp <- pmin(9.9, nu_lp)}
  if (min(nu_lp) < 0.001){nu_lp <- pmax(0.001, nu_lp)}
  log_lambda <- log(mapply(findlam, mu, nu_lp))
  # Proposal
  theta_p <- rnorm(pplayers, theta, 0.25) # 0.1*sig_prop
  mu_p <- with(data, exp(X%*%beta + theta_p[ID]+ XOpp%*%zeta + 
                           xi[Inns] + phi[HA] + gamma[TossInns1]))
  if (max(mu_p) >= 10){mu_p <- pmin(9.9, mu_p)}
  if (min(mu_p) < 0.001){mu_p <- pmax(0.001, mu_p)}
  log_lambda_p <- log(mapply(findlam, mu_p, nu_lp))
  lG <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(nu_lp, lfac)))) 
  lG_p <- log(rowSums(exp(tcrossprod(log_lambda_p, unique) - tcrossprod(nu_lp, lfac)))) 
  # Component-wise update (vectorised) 
  denrat <- with(data, tapply(Wkts*(log_lambda_p - log_lambda) + lG - lG_p, ID, sum))
  denrat <- denrat - 0.5*((theta_p - theta_m)^2 - (theta - theta_m)^2)/theta_s^2
  if (sum(is.na(denrat)) != 0){
    denrat[is.na(denrat)] <- 0
  }
  laccept <- pmin(0, denrat)
  accept <- (log(runif(pplayers)) < laccept)
  theta[accept] <- theta_p[accept]
  theta <- theta - mean(theta)
  list(theta = theta)
} 

# Opposition effects (time-varying due to spline structure)
zetaUpdate <- function(data, X, XOpp, OppInd, zeta_m, zeta_s, popps, unique, 
                       lfac, beta, nu, theta, zeta, xi, phi, gamma){
  mu <- with(data, exp(X%*%beta + theta[ID] + XOpp%*%zeta + 
                         xi[Inns] + phi[HA] + gamma[TossInns1]))
  if (max(mu) >= 10){mu <- pmin(9.9, mu)}
  if (min(mu) < 0.001){mu <- pmax(0.001, mu)}
  nu_lp <- exp(nu[data$ID])
  if (max(nu_lp >= 10)){nu_lp <- pmin(9.9, nu_lp)}
  if (min(nu_lp) < 0.001){nu_lp <- pmax(0.001, nu_lp)}
  log_lambda <- log(mapply(findlam, mu, nu_lp))
  # Proposal
  zeta_p <- rnorm(popps, zeta, 0.1) 
  #zeta_p <- as.vector(rmvnorm(1, zeta, 0.01*sig_prop)) 
  mu_p <- with(data, exp(X%*%beta + theta[ID] + XOpp%*%zeta_p + 
                           xi[Inns] + phi[HA] + gamma[TossInns1]))
  if (max(mu_p) >= 10){mu_p <- pmin(9.9, mu_p)}
  if (min(mu_p) < 0.001){mu_p <- pmax(0.001, mu_p)}
  log_lambda_p <- log(mapply(findlam, mu_p, nu_lp))
  lG <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(nu_lp, lfac)))) 
  lG_p <- log(rowSums(exp(tcrossprod(log_lambda_p, unique) - tcrossprod(nu_lp, lfac)))) 
  # Component-wise update (vectorised)
  denrat <- colSums((data$Wkts*(log_lambda_p - log_lambda) - (lG_p - lG))*OppInd)
  denrat <- denrat - 0.5*((zeta_p - zeta_m)^2 - (zeta - zeta_m)^2)/zeta_s^2
  laccept <- pmin(0, denrat)
  accept <- (log(runif(popps)) < laccept)
  zeta[accept] <- zeta_p[accept]
  list(zeta = zeta)
} 

# Innings effects
xiUpdate <- function(data, X, XOpp, xi_m, xi_s, unique, lfac, 
                     beta, nu, theta, zeta, xi, phi, gamma){
  mu <- with(data, exp(X%*%beta + theta[ID] + XOpp%*%zeta + 
                         xi[Inns] + phi[HA] + gamma[TossInns1]))
  if (max(mu) >= 10){mu <- pmin(9.9, mu)}
  if (min(mu) < 0.001){mu <- pmax(0.001, mu)}
  nu_lp <- exp(nu[data$ID])
  if (max(nu_lp >= 10)){nu_lp <- pmin(9.9, nu_lp)}
  if (min(nu_lp) < 0.001){nu_lp <- pmax(0.001, nu_lp)}
  log_lambda <- log(mapply(findlam, mu, nu_lp))
  # Proposal
  xi_p <- rnorm(4, xi, 0.03)
  mu_p <- with(data, exp(X%*%beta + theta[ID]+ XOpp%*%zeta + 
                           xi_p[Inns] + phi[HA] + gamma[TossInns1]))
  if (max(mu_p) >= 10){mu_p <- pmin(9.9, mu_p)}
  if (min(mu_p) < 0.001){mu_p <- pmax(0.001, mu_p)}
  log_lambda_p <- log(mapply(findlam, mu_p, nu_lp))
  lG <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(nu_lp, lfac)))) 
  lG_p <- log(rowSums(exp(tcrossprod(log_lambda_p, unique) - tcrossprod(nu_lp, lfac)))) 
  # Component-wise update (vectorised)  - check the below
  denrat <- with(data, tapply(Wkts*(log_lambda_p - log_lambda) + lG - lG_p, Inns, sum))
  denrat <- denrat - 0.5*((xi_p - xi_m)^2 - (xi - xi_m)^2)/xi_s^2
  laccept <- pmin(0, denrat)
  accept <- (log(runif(4)) < laccept)
  xi[accept] <- xi_p[accept]
  xi[1] <- 0
  list(xi = xi)
} 

# Home/away effects
phiUpdate <- function(data, X, XOpp, phi_m, phi_s, unique, lfac, 
                      beta, nu, theta, zeta, xi, phi, gamma){
  mu <- with(data, exp(X%*%beta + theta[ID] + XOpp%*%zeta + 
                         xi[Inns] + phi[HA] + gamma[TossInns1]))
  if (max(mu) >= 10){mu <- pmin(9.9, mu)}
  if (min(mu) < 0.001){mu <- pmax(0.001, mu)}
  nu_lp <- exp(nu[data$ID])
  if (max(nu_lp >= 10)){nu_lp <- pmin(9.9, nu_lp)}
  if (min(nu_lp) < 0.001){nu_lp <- pmax(0.001, nu_lp)}
  log_lambda <- log(mapply(findlam, mu, nu_lp))
  # Proposal
  phi_p <- rnorm(2, phi, 0.02)
  mu_p <- with(data, exp(X%*%beta + theta[ID]+ XOpp%*%zeta + 
                           xi[Inns] + phi_p[HA] + gamma[TossInns1]))
  if (max(mu_p) >= 10){mu_p <- pmin(9.9, mu_p)}
  if (min(mu_p) < 0.001){mu_p <- pmax(0.001, mu_p)}
  log_lambda_p <- log(mapply(findlam, mu_p, nu_lp))
  lG <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(nu_lp, lfac)))) 
  lG_p <- log(rowSums(exp(tcrossprod(log_lambda_p, unique) - tcrossprod(nu_lp, lfac)))) 
  # Component-wise update (vectorised) 
  denrat <- with(data, tapply(Wkts*(log_lambda_p - log_lambda) + lG - lG_p, HA, sum))
  denrat <- denrat - 0.5*((phi_p - phi_m)^2 - (phi - phi_m)^2)/phi_s^2
  laccept <- pmin(0, denrat)
  accept <- (log(runif(2)) < laccept)
  phi[accept] <- phi_p[accept]
  phi[1] <- 0
  list(phi = phi)
} 

# Toss x first innings interaction update
gammaUpdate <- function(data, X, XOpp, gamma_m, gamma_s, unique, lfac, 
                        beta, nu, theta, zeta, xi, phi, gamma){
  mu <- with(data, exp(X%*%beta + theta[ID] + XOpp%*%zeta + 
                         xi[Inns] + phi[HA] + gamma[TossInns1]))
  if (max(mu) >= 10){mu <- pmin(9.9, mu)}
  if (min(mu) < 0.001){mu <- pmax(0.001, mu)}
  nu_lp <- exp(nu[data$ID])
  if (max(nu_lp >= 10)){nu_lp <- pmin(9.9, nu_lp)}
  if (min(nu_lp) < 0.001){nu_lp <- pmax(0.001, nu_lp)}
  log_lambda <- log(mapply(findlam, mu, nu_lp))
  # Proposal
  gamma_p <- rnorm(2, gamma, 0.02)
  mu_p <- with(data, exp(X%*%beta + theta[ID]+ XOpp%*%zeta + 
                           xi[Inns] + phi[HA] + gamma_p[TossInns1]))
  if (max(mu_p) >= 10){mu_p <- pmin(9.9, mu_p)}
  if (min(mu_p) < 0.001){mu_p <- pmax(0.001, mu_p)}
  log_lambda_p <- log(mapply(findlam, mu_p, nu_lp))
  lG <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(nu_lp, lfac)))) 
  lG_p <- log(rowSums(exp(tcrossprod(log_lambda_p, unique) - tcrossprod(nu_lp, lfac)))) 
  # Component-wise update (vectorised) 
  denrat <- with(data, tapply(Wkts*(log_lambda_p - log_lambda) + lG - lG_p, TossInns1, sum))
  denrat <- denrat - 0.5*((gamma_p - gamma_m)^2 - (gamma - gamma_m)^2)/gamma_s^2
  laccept <- pmin(0, denrat)
  accept <- (log(runif(2)) < laccept)
  gamma[accept] <- gamma_p[accept]
  gamma[1] <- 0
  #print(gamma)
  list(gamma = gamma)
} 

# Function to calculate log-likelihood
llhood <- function(data, X, XOpp, unique, lfac, lfacWkts, 
                   beta, nu, theta, zeta, xi, phi, gamma){
  mu <- with(data, exp(X%*%beta + theta[ID] + XOpp%*%zeta + xi[Inns] + 
                         phi[HA] + gamma[TossInns1]))
  if (max(mu) >= 10){mu <- pmin(9.9, mu)}
  if (min(mu) < 0.001){mu <- pmax(0.001, mu)}
  nu_lp <- exp(nu[data$ID])
  if (max(nu_lp > 10)){nu_lp <- pmin(9.9, nu_lp)}
  if (min(nu_lp) < 0.001){nu_lp <- pmax(0.001, nu_lp)}
  log_lambda <- log(lam_grid_poly_10K[cbind(round(mu*1000), round(nu_lp*1000))])
  lG <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(nu_lp, lfac)))) 
  ll <- sum(data$Wkts*log_lambda - nu_lp*lfacWkts - lG)
  list(ll = ll)
}

