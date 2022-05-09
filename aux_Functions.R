# jonashaslbeck@gmail.com; May 9, 2022

# ------------------------------------------------------
# -------- What is this? -------------------------------
# ------------------------------------------------------

# Auxiliary functions

# ------------------------------------------------------
# -------- Load Packages -------------------------------
# ------------------------------------------------------

library(multimode)

# ------------------------------------------------------
# -------- Modality Detection via Excess Mass ----------
# ------------------------------------------------------

ExcessMass <- function(x, pcrit=0.1, B=20, maxM=10) {

  # Apply Excess Mass method
  pval <- 0
  counter <- 1
  l_out <- list()
  while(pval < pcrit) {
    l_out[[counter]] <- modetest(x, method="ACR", mod0=counter, B=B)
    pval <- l_out[[counter]]$p.value
    counter <- counter + 1
    if(counter>maxM) pval <- 1 # stop after 10
  } # end while

  # Get M-estimate
  M <- counter-1

  outlist <- list("objects"=l_out,
                  "M"=M)
  return(outlist)

} # eoF


# ------------------------------------------------------
# -------- Modality Detection via Silverman ------------
# ------------------------------------------------------

# Input: data
# Output: number of modes + density estimate

Silverman <- function(x, pcrit=0.1, B=500, maxM=10) {

  # Apply Silverman's method
  pval <- 0
  counter <- 1
  l_out <- list()
  while(pval < pcrit) {
    l_out[[counter]] <- modetest(x, method="SI", mod0=counter, B=B)
    pval <- l_out[[counter]]$p.value
    counter <- counter + 1
    if(counter>maxM) pval <- 1 # stop after 10
  } # end while

  # Get M-estimate
  M <- counter-1

  # Detect scale
  x_s <- seq(1, length(unique(x)), length=100)
  if(max(x) > 8) x_s <- seq(0, 100, length=100) # better would be to input and hardcode all of them

  # Get density estimate between M and M+1 [just for viz]
  Xden <- matrix(NA, length(x), 100)

  if(M>1) {
    sd_dis <- (l_out[[counter-2]]$statistic + l_out[[counter-1]]$statistic) / 2
  } else {
    sd_dis <- l_out[[counter-1]]$statistic * 2
  }

  for(s in 1:length(x)) Xden[s,] <- dnorm(x_s, mean=x[s], sd = sd_dis)
  XdenS <- colMeans(Xden)

  # Return
  outlist <- list("M"=M,
                  "dens_x"=x_s,
                  "dens_y"=XdenS)
  return(outlist)

} # eoF





# ------------------------------------------------------
# -------- Plotting Aux Functions ----------------------
# ------------------------------------------------------

# Plot Labels
PL <- function(tex, cex = 1.5, srt=0, x=0.45, y=0.5) {

  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(x, y, labels = tex, srt=srt, cex=cex, adj=0.25)

}

# ------------------------------------------------------
# -------- MM Detection with density E -----------------
# ------------------------------------------------------

GMM_MMdet <- function(X, noise=0.4, n=100) {

  # Add noise
  Xn <- X + rnorm(n=length(X), 0, noise)

  ## -- GMM + compute densities --
  out <- Mclust(Xn, modelNames = "V", verbose=FALSE)
  # compute mixture density
  x <- seq(min(Xn), max(Xn), length=n)
  l_dens_GMM <- list()
  for(k in 1:out$G) l_dens_GMM[[k]] <- dnorm(x, mean=out$parameters$mean[k],
                                             sd = sqrt(out$parameters$variance$sigmasq[k]))
  dens_mix <- colSums(do.call(rbind, l_dens_GMM))
  # Count roots to get to modality
  ni <- length(dens_mix)
  diff <- dens_mix[-ni] - dens_mix[-1]
  sign_reversals <- sign(diff[-(ni-1)]) != sign(diff[-1])
  Nrev <- sum(sign_reversals)
  M <- ceiling(Nrev/2)

  # Outlist
  outlist <- list("M"=M,
                  "den_x" = x,
                  "den_y" = dens_mix)

} # eoF

# ------------------------------------------------------
# -------- MM Detection via GMM + roots ----------------
# ------------------------------------------------------

DensMMdet <- function(X, adjust=3, noise=0.4, n = 100, plot=FALSE) {

  # Add noise
  Xn <- X + rnorm(length(X), 0, noise)

  # Density Estimation
  den <- density(Xn, bw="SJ", adjust=adjust, n=n)
  if(plot) plot(den)

  # Compute number of reversals
  ni <- length(den$y) # usually 512

  diff <- den$y[-ni] - den$y[-1]

  sign_reversals <- sign(diff[-(ni-1)]) != sign(diff[-1])
  Nrev <- sum(sign_reversals)

  Modality <- ceiling(Nrev/2) # since between each mode there is another reversal

  outlist <- list("M" = Modality,
                  "den_x" = den$x,
                  "den_y" = den$y)

  return(outlist)

} # eoF



# ------------------------------------------------------
# -------- Data Generation: Continuous -----------------
# ------------------------------------------------------

DataGen_cont <- function(type, N) {

  # unimodal
  if(type == 1) {
    data <- rnorm(n = N, mean=0, sd=1)
  }
  # unimodal skewed
  if(type == 2) {
    data2.1 <- rnorm(n = round(N*1/2), mean=0, sd=1)
    data2.2 <- rnorm(n = round(N*1/2), mean=3.5, sd=4)
    data <- c(data2.1, data2.2)
  }
  ## Bimodal
  if(type == 3) {
    data3.1 <- rnorm(n = round(N/2), mean=0, sd=1)
    data3.2 <- rnorm(n = round(N/2), mean=4, sd=1)
    data <- c(data3.1, data3.2)
  }
  ## 3 Modes
  if(type == 4) {
    data4.1 <- rnorm(n = round(N/3), mean = -1, sd=1)
    data4.2 <- rnorm(n = round(N/3), mean = 3, sd=1)
    data4.3 <- rnorm(n = round(N/3), mean = 7, sd=1)
    data <- c(data4.1, data4.2, data4.3)
  }

  return(data)

} # eoF


# ------------------------------------------------------
# -------- Data Generation: 100-scale ------------------
# ------------------------------------------------------


DataGen_100 <- function(type, N, l_d) {

  J <- 100

  # Convert GMM densities in multinomial probabilities
  v_p <- l_d[[type]][seq(1, length(l_d[[type]]), length=J)]

  # sample from multinomial
  out <- t(rmultinom(N, size=1, prob=v_p))
  data <- apply(out, 1, function(x) which(x==1))

  return(data)

} # eoF


# ------------------------------------------------------
# -------- Data Generation: Likert-Scale ---------------
# ------------------------------------------------------

DataGen_ord <- function(type, N) {

  # unimodal
  if(type == 1) {
    out <- t(rmultinom(N, size=1, prob=c(.1, .2, .3, .4, .3, .2, .1)))
    data <- apply(out, 1, function(x) which(x==1))
  }
  # unimodal skewed
  if(type == 2) {
    out <- t(rmultinom(N, size=1, prob=c(.3, .6, .3, .2, .1, .05, .01)))
    data <- apply(out, 1, function(x) which(x==1))
  }
  ## Bimodal
  if(type == 3) {
    out <- t(rmultinom(N, size=1, prob=c(.1, .6, .3, .1, .3, .6, .1)))
    data <- apply(out, 1, function(x) which(x==1))
  }
  ## 3 Modes
  if(type == 4) {
    out <- t(rmultinom(N, size=1, prob=c(.7, .2, .3, .6, .3, .2, .7)))
    data <- apply(out, 1, function(x) which(x==1))
  }

  return(data)

}

