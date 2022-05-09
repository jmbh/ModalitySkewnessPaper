# jonashaslbeck@gmail.com; May 9, 2022

# ------------------------------------------------------
# -------- Load Packages -------------------------------
# ------------------------------------------------------

# Plotting
library(scales)
library(RColorBrewer)

# Data Generation
library(MASS)

# Modality detection
library(mclust)
library(mousetrap)
library(diptest)
library(moments)
library(multimode)
source("aux_Functions.R")

# ------------------------------------------------------
# -------- Load Data -----------------------------------
# ------------------------------------------------------

baseDir <- "Data/"

data_files <- c("data_Rowland2020.RDS",
                "data_Bringmann2016.RDS",
                "data_Vrijen2018.RDS",
                "data_Fisher2017.RDS",
                "data_Bringmann2013.RDS",
                "data_Fried2021.RDS",
                "data_Wendt2020.RDS")

# Subsetting of emotion variables in each dataset
l_ind_vars <- list("Rowland2020" = 5:12, # 0-100 scale
                   "Bringmann2016" = 2:7, # 0-100
                   "Vrijen2018" = c(5,6,7,8), # 0-100
                   "Fisher2017" = c(2:16, 19:20), # 0-100
                   "Bringmann2013" = c(6,8:11), # 1-7
                   "Fried2021" = c(2:8,10:11), # 1-5
                   "Wendt2020" = 9:39) # 1-5 scale

# Get Rowland dataset
v <- 1
data_Rowland <- readRDS(paste0(baseDir, names(l_ind_vars)[v], "/", data_files[v]))
var_set <- l_ind_vars[[v]]
data_Rowland_ss <- data_Rowland[, var_set]
head(data_Rowland_ss)

# Get Bringmann2013 / MM dataset
v <- 5
data_MM <- readRDS(paste0(baseDir, names(l_ind_vars)[v], "/", data_files[v]))
var_set <- l_ind_vars[[v]]
data_MM_ss <- data_MM[, var_set]
head(data_MM_ss)


# ------------------------------------------------------
# -------- Validation 1: Empirical Examples ------------
# ------------------------------------------------------

# ---------- A) Rowland Data ---------

# Select 9 representative examples from Rowland
u_ptp <- unique(data_Rowland$subj_id)

# unimodal OK
R3 <- data_Rowland_ss[data_Rowland$subj_id==u_ptp[3], 3] # unimodal, OK
R5 <- data_Rowland_ss[data_Rowland$subj_id==u_ptp[3], 7] # unimodal skew, OK

# bimodal OK
R1 <- data_Rowland_ss[data_Rowland$subj_id==u_ptp[1], 1] # bimodal, OK
R9  <- data_Rowland_ss[data_Rowland$subj_id==u_ptp[1], 8] # multi-modal, OKish

# bimodal bad
R2 <- data_Rowland_ss[data_Rowland$subj_id==u_ptp[2], 3] # bimodal, bad
R7  <- data_Rowland_ss[data_Rowland$subj_id==u_ptp[4], 2] # bimodal, bad
R8  <- data_Rowland_ss[data_Rowland$subj_id==u_ptp[6], 2] # multi-modal, bad

# not so clear
R4 <- data_Rowland_ss[data_Rowland$subj_id==u_ptp[21], 1] # not so clear
R6  <- data_Rowland_ss[data_Rowland$subj_id==u_ptp[14], 5] # not so clear


l_R <- list(R3, R5, R6,  # unimodal / skewed
            R7, R4, R1,
            R8, R2, R9)

# --------- Estimation --------

# Storage: M-estimates
m_estM <- matrix(NA, 9, 6)
l_den <- list()

# Reproducibility
set.seed(1)

for(i in 1:9) {

  print(i)
  X_i <- na.omit(l_R[[i]])

  ## ----- Bimodality coefficient -----
  m_estM[i, 1] <- bimodality_coefficient(X_i) > 5/9

  ## ----- Hardigan's dip ------
  dip_out <- dip.test(X_i)
  m_estM[i, 2] <- dip_out$p.value < 0.05

  ## ----- GMM + compute densities ------
  out <- Mclust(X_i, modelNames = "V")
  # compute mixture density
  x <- seq(0, 100, length=100)
  l_dens_GMM <- list()
  for(k in 1:out$G) l_dens_GMM[[k]] <- dnorm(x, mean=out$parameters$mean[k],
                                             sd = sqrt(out$parameters$variance$sigmasq[k]))
  dens_mix <- colSums(do.call(rbind, l_dens_GMM))
  # Count roots to get to modality
  ni <- length(dens_mix)
  diff <- dens_mix[-ni] - dens_mix[-1]
  sign_reversals <- sign(diff[-(ni-1)]) != sign(diff[-1])
  Nrev <- sum(sign_reversals)
  m_estM[i, 3] <- ceiling(Nrev/2)

  ## ----- Our method -----
  prop_modecat <- (max(table(X_i))/length(X_i))
  base_noise <- 0.035
  noise_adj0 <- base_noise + prop_modecat / 5
  noise <- noise_adj0 * 100

  # RERUN 10 times to reduce arbitrariness of specific noise draws
  v_M <- rep(NA, 10)
  l_M <- list()
  for(r in 1:10) { # do 10 reruns
    l_M[[r]] <- DensMMdet(X_i, n=100, noise=noise, adjust=2)
    v_M[r] <- l_M[[r]]$M
  }
  tb <- table(v_M)
  m_estM[i, 4] <- as.numeric(names(tb)[which.max(tb)]) # if equal, take smaller one
  l_den[[i]] <- l_M[[which(v_M == m_estM[i, 4])[1] ]] # save density of first one with selected #modes

  ## ----- Silvermann ------
  out_S <- Silverman(X_i, pcrit = 0.1, B=500)
  m_estM[i, 5] <- out_S$M

  ## ----- Excess Mass ------
  out_EM <- ExcessMass(X_i, pcrit = 0.05, B=10, maxM = 10)
  m_estM[i, 6] <- out_EM$M
}

## ----- Plotting ------
library(scales)
sc <- 0.9
pdf("Figures/Validation_P1_D1.pdf", width = 12*sc, height = 7.5*sc)

# Setup layout
lmat <- matrix(1:18, nrow=3, byrow = TRUE)
w_l <- .5
lo <- layout(lmat, widths = c(1, w_l, 1, w_l, 1, w_l))
# layout.show(lo)

for(i in 1:9) {

  # Plot Histogram
  par(mar=c(3,5,2,0))
  hout <- hist(l_R[[i]], breaks=seq(0,100, length=30), axes=FALSE, xlab="", main="")
  axis(1, seq(0,100, length=5))
  axis(2, las=2)
  den_out <- l_den[[i]]
  lines(den_out$den_x, (den_out$den_y/ max(den_out$den_y)) * max(hout$counts),
        col="black", lwd=2)
  title(main=paste0("(", letters[i],") "), font.main=1)

  # Plot Classification Results
  methods_label <- c("Bim-coef", "Hardigan's d", "GMM", "Our method", "Silverman", "Excess-mass")
  cex <- 1
  col <- "black"
  y_locs <- seq(0.7, .35, length=6)
  x_loc <- 0
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  for(s in 1:2) {
    if(m_estM[i, s]==0) text_d <- paste0(methods_label[s], ": M = 1")
    if(m_estM[i, s]==1) text_d <- paste0(methods_label[s], ": M > 1")
    text(x_loc, y_locs[s], text_d, adj=0, cex=cex, col=col)
  }
  for(s in 3:6) text(x_loc, y_locs[s], paste0(methods_label[s], ": M = ",  m_estM[i, s]),
                     adj=0, cex=cex, col=col)

} # end for: 9 examples

dev.off()


# ---------- B) Bringmann2013 Data ---------

# Select 9 representative examples from Bringmann2013

u_ptp <- unique(data_MM$subj_id)

# Clear unimodal (1)
R1 <- data_MM_ss[data_MM$subj_id==u_ptp[31], 2]

# Unimodal, skew (2)
R2 <- data_MM_ss[data_MM$subj_id==u_ptp[6], 2]
R3 <- data_MM_ss[data_MM$subj_id==u_ptp[1], 1]

# Bi/multimodal (incorrectly classified) (2)
R4 <- data_MM_ss[data_MM$subj_id==u_ptp[9], 2]
R5 <- data_MM_ss[data_MM$subj_id==u_ptp[19], 2]

# Bi/multimodal (correctly classified) (3)
R6 <- data_MM_ss[data_MM$subj_id==u_ptp[1], 2]
R7 <- data_MM_ss[data_MM$subj_id==u_ptp[41], 2]
R8 <- data_MM_ss[data_MM$subj_id==u_ptp[17], 2]

# not so clear (1)
R9 <- data_MM_ss[data_MM$subj_id==u_ptp[34], 4]

l_R2 <- list(R1, R3,R2, R9, # unimodal / skewed / unclear
             R4, R5, # multimodal (incorrect)
             R6, R7, R8) # multimodal (correct)

# --------- Estimation --------

# Storage: M-estimates
m_estM2 <- matrix(NA, 9, 6)
l_den2 <- list()

# Reproducibility
set.seed(1)

for(i in 1:9) {

  print(i)
  X_i <- na.omit(l_R2[[i]])

  ## ----- Bimodality coefficient -----
  m_estM2[i, 1] <- bimodality_coefficient(X_i) > 5/9

  ## ----- Hardigan's dip ------
  dip_out <- dip.test(X_i)
  m_estM2[i, 2] <- dip_out$p.value < 0.05

  ## ----- GMM + compute densities ------
  out <- Mclust(X_i, modelNames = "V")
  # compute mixture density
  x <- seq(0, 100, length=100)
  l_dens_GMM <- list()
  for(k in 1:out$G) l_dens_GMM[[k]] <- dnorm(x, mean=out$parameters$mean[k],
                                             sd = sqrt(out$parameters$variance$sigmasq[k]))
  dens_mix <- colSums(do.call(rbind, l_dens_GMM))
  print(out$G)
  # Count roots to get to modality
  ni <- length(dens_mix)
  diff <- dens_mix[-ni] - dens_mix[-1]
  sign_reversals <- sign(diff[-(ni-1)]) != sign(diff[-1])
  Nrev <- sum(sign_reversals)
  m_estM2[i, 3] <- ceiling(Nrev/2)

  ## ----- Our method -----
  prop_modecat <- (max(table(X_i))/length(X_i))
  noise <- prop_modecat

  # RERUN 10 times to reduce arbitrariness of specific noise draws
  v_M <- rep(NA, 10)
  l_M <- list()
  for(r in 1:10) { # do 10 reruns
    l_M[[r]] <- DensMMdet(X_i, n=100, noise=noise, adjust=2)
    v_M[r] <- l_M[[r]]$M
  }
  tb <- table(v_M)
  m_estM2[i, 4] <- as.numeric(names(tb)[which.max(tb)]) # if equal, take smaller one
  l_den2[[i]] <- l_M[[which(v_M == m_estM2[i, 4])[1] ]] # save density of first one with selected #modes

  ## ----- Silvermann ------
  out_S <- Silverman(X_i, pcrit = 0.1, B=500)
  m_estM2[i, 5] <- out_S$M

  ## ----- Excess Mass ------
  out_EM <- ExcessMass(X_i, pcrit = 0.05, B=10, maxM = 10)
  m_estM2[i, 6] <- out_EM$M

} # end loop: over 9 examples


sc <- 0.9
pdf("Figures/Validation_P1_D2.pdf", width = 12*sc, height = 7.5*sc)

# Setup layout
lmat <- matrix(1:18, nrow=3, byrow = TRUE)
w_l <- .5
lo <- layout(lmat, widths = c(1, w_l, 1, w_l, 1, w_l))
# layout.show(lo)

for(i in 1:9) {

  # Plot Histogram
  par(mar=c(3,5,2,0))
  hout <- hist(l_R2[[i]], breaks=seq(1,7, length=30), axes=FALSE, xlab="", main="")
  axis(1, seq(1, 7, length=7))
  axis(2, las=2)
  den_out <- l_den2[[i]]
  lines(den_out$den_x, (den_out$den_y/ max(den_out$den_y)) * max(hout$counts),
        col="black", lwd=2)
  title(main=paste0("(", letters[i],") "), font.main=1)

  # Plot Classification Results
  methods_label <- c("Bim-coef", "Hardigan's d", "GMM", "Our method", "Silverman", "Excess-mass")
  cex <- 1
  col <- "black"
  y_locs <- seq(0.7, .35, length=6)
  x_loc <- 0
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  for(s in 1:2) {
    if(m_estM2[i, s]==0) text_d <- paste0(methods_label[s], ": M = 1")
    if(m_estM2[i, s]==1) text_d <- paste0(methods_label[s], ": M > 1")
    text(x_loc, y_locs[s], text_d, adj=0, cex=cex, col=col)
  }
  for(s in 3:6) text(x_loc, y_locs[s], paste0(methods_label[s], ": M = ",  m_estM2[i, s]),
                     adj=0, cex=cex, col=col)

} # end for: 9 examples

dev.off()




# ------------------------------------------------------
# -------- Validation 2: Simulation Study --------------
# ------------------------------------------------------

# Vary:
# - modality: unimodal, unimodal skewed, bimodal, 3 modes
# - sample size: 20 - 200

# -------- Select Continuous Distributions --------

N <- 1000
x <- seq(-4,10, length=1000)

## Unimodal
data1 <- rnorm(n = N, mean=0, sd=1)
d1 <- dnorm(x, mean = 0, sd=1)

## Unimodal skewed
data2.1 <- rnorm(n = N, mean=0, sd=1)
data2.2 <- rnorm(n = N/3, mean=3.5, sd=4)
data2 <- c(data2.1, data2.2)
d2.1 <- dnorm(x, mean = 0, sd=1)
d2.2 <- dnorm(x, mean = 3.5, sd=4)
d2 <- d2.1*1/2 + d2.2*1/2
plot(d2, type="l")

## Bimodal
data3.1 <- rnorm(n = N, mean=0, sd=1)
data3.2 <- rnorm(n = N, mean=3.5, sd=1)
data3 <- c(data3.1, data3.2)
d3.1 <- dnorm(x, mean = 0, sd=1)
d3.2 <- dnorm(x, mean = 4.5, sd=1)
d3 <- d3.1*0.5 + d3.2*0.5
plot(d3)

## 3 Modes
data4.1 <- rnorm(n = N, mean=0, sd=1)
data4.2 <- rnorm(n = N, mean=3.5, sd=1)
data4.3 <- rnorm(n = N, mean=7, sd=1)
data4 <- c(data4.1, data4.2, data4.3)
d4.1 <- dnorm(x, mean = -1.5, sd=1)
d4.2 <- dnorm(x, mean = 3, sd=1)
d4.3 <- dnorm(x, mean = 7.5, sd=1)
d4 <- d4.1*1/3 + d4.2*1/3 + d4.3*1/3
plot(d4)

l_d <- list(d1, d2, d3, d4)


# -------- Simulation Settings --------

# Basics
n_seq <- seq(20, 200, length=10)
n_n_seq <- length(n_seq)
nIter <- 500

# ---- Storage ----
# M-estimation
a_out <- array(NA, dim=c(nIter, 4, n_n_seq, 2, 4)) # nIter, 4 types, n-var, domain, 4 methods
# Skewness
a_out_skew <- array(NA, dim=c(nIter, 2, n_n_seq, 2)) # only two types
# -------- Run Simulation --------

set.seed(1)

for(i in 1:nIter) { # iterations
  for(type in 1:4) { # type of modality
    for(n in 1:n_n_seq) { # n- variations
      for(dom in 1:2) { # domain: continuous vs. ordinal

        # ---- Generate data ----

        N <- n_seq[n]

        if(dom==1) {
          data <- DataGen_100(type = type, N = N, l_d = l_d)
        } else {
          data <- DataGen_ord(type = type, N = N)
        }

        # ---- Estimate Modality ----

        ## -- Our Method --
        # Set Noise hyper-parameter based on scale
        u_cat <- length(unique(data))
        range_i_j <- max(data) - min(data)
        max_th_range <- 14 # roughly the range that includes almost all data
        prop_modecat <- (max(table(data))/length(data))

        if(dom==1) { # continuous case
          base_noise <- 0.035
          noise_adj0 <- base_noise + prop_modecat / 5
          noise <- noise_adj0 * max_th_range
        } else { # ordinal case
          noise <- prop_modecat
        }

        # Set bandwidth hyper-parameter conditional on scale
        adjust <- ifelse(dom==1, 2, 3) # lower adjust (i.e., lower bandwidth) for continuous data

        den_out <- DensMMdet(X = data, noise = noise, adjust=adjust) # playaround
        a_out[i, type, n, dom, 1] <- ifelse(den_out$M>1, 1, 0)

        ## -- Bimodality coefficient --
        a_out[i, type, n, dom, 2] <- bimodality_coefficient(data) > 5/9

        ## -- Hardigan's dip --
        dip_out <- dip.test(data)
        a_out[i, type, n, dom, 3] <- dip_out$p.value < 0.05

        ## -- GMM + compute densities --
        GMM_out <- GMM_MMdet(X=data, noise=noise, n=100)
        a_out[i, type, n, dom, 4] <- ifelse(GMM_out$M>1, 1, 0)

        # For types 1&2: do skewness detection
        if(type %in% 1:2) {
          a_out_skew[i, type, n, dom] <- skewness(data)
        }
      }
    }
  }
  print(i)
} # end: loops


# SAVE
saveRDS(a_out, file="Files/SimResults_Validation.RDS")
saveRDS(a_out_skew, file="Files/SimResults_Validation_Skew.RDS")

# -------- Evaluate / Make Figure --------

titles <- c("Unimodal", "Unimodal skewed", "Bimodal", "3 Modes")
cols <- brewer.pal(4, "Set1")

# --- Pre-processing ---
# Turn everything into unimodal vs. bimodal
a_out_bin <- a_out
v_correct <- c(0,0,1,1)
for(type in 1:4) a_out_bin[, type, , , ] <- a_out[, type, , , ] == v_correct[type]
# Aggregate over iterations
a_out_bin_agg <- apply(a_out_bin, 2:5, mean)

a_out_bin_agg[type, , 1, 1]

# ----- Figure: continuous -----
sc <- .8
pdf("figures/Validation_P2_D1.pdf", width=9*sc, height=4*sc)

dom <- 1 # select results for continuous domain

lmat <- rbind(1:4, 5:8)
layout(mat = lmat, heights = c(.6,1))

par(mar=c(1,3,1,1))
## Plot true probability distribution
for(type in 1:4) {
  v_p <- l_d[[type]][seq(1, length(l_d[[type]]), length=30)]
  barplot(v_p, axes=FALSE, border="black")
  title(main = titles[type], font.main=1, cex.main=1.3, line=1.5)
}

## Plot data
par(mar=c(4,3,1,1))

for(type in 1:4) {
  plot.new()
  plot.window(xlim=c(1,10), ylim=c(0,1))
  axis(1, labels = n_seq, at=1:10, las=2)
  axis(2, las=2)
  title(ylab="Proportion Correct")
  title(xlab="Sample size", line=3)

  # data / methods
  for(j in 1:4) lines(a_out_bin_agg[type, , dom, j], col=cols[j], lwd=2)

  # legend
  if(type==1) legend("right", legend = c("Density-based", "BC", "H's Dip", "GMM"),
                     col=cols, lwd=rep(2,4), bty="n")

}

dev.off()


# ----- Figure: ordinal -----
dom <- 2 # select results for ordinal domain

prob_table <- rbind(c(.1, .2, .3, .4, .3, .2, .1),
                    c(.3, .6, .3, .2, .1, .05, .01),
                    c(.1, .6, .3, .1, .3, .6, .1),
                    c(.7, .2, .3, .6, .3, .2, .7))

sc <- .7
pdf("figures/Validation_P2_D2.pdf", width=9*sc, height=4*sc)

lmat <- rbind(1:4, 5:8)
layout(mat = lmat, heights = c(.6,1))

par(mar=c(0,3,3,1))
## Plot true densities
for(type in 1:4) {
  barplot(prob_table[type, ], axes=FALSE)
  title(main = titles[type], font.main=1, cex.main=1.3, line=1.5)
}

## Plot data
par(mar=c(4,3,1,1))

for(type in 1:4) {
  plot.new()
  plot.window(xlim=c(1,10), ylim=c(0,1))
  axis(1, labels = n_seq, at=1:10, las=2)
  axis(2, las=2)
  title(ylab="Proportion Correct")
  title(xlab="Sample size", line=3)

  # data
  for(j in 1:4) lines(a_out_bin_agg[type, , dom, j], col=cols[j], lwd=2)

  # legend
  if(type==1) legend("right", legend = c("Density-based", "BC", "H's Dip", "GMM"),
                     col=cols, lwd=rep(2,4), bty="n")

}

dev.off()

# -------- Evaluate / Make Figure [Skewness] --------

n_seq <- seq(20, 200, length=10)

cutoff <- 2/3 # Defined cutoff

a_out_skew_p <- apply(a_out_skew, 1:4, function(x) x>cutoff)
a_out_skew_p[, 1, , ] <- a_out_skew_p[, 1, , ] == 0
a_out_skew_p[, 2, , ] <- a_out_skew_p[, 2, , ] == 1
a_out_skew_agg <- apply(a_out_skew_p, 2:4, mean)

# ----- Make Figure -----
titles <- c("Unimodal", "Unimodal skewed", "Bimodal", "3 Modes")

sc <- 0.9
pdf("figures/Validation_Skew.pdf", width=9*sc, height=4*sc)

lmat <- matrix(1:8, nrow=2, byrow = T)
layout(lmat, heights = c(.8, 1))

## Top labels
par(mar=c(1,4,3,1))
# Continuous
for(type in 1:2) {
  v_p <- l_d[[type]][seq(1, length(l_d[[type]]), length=30)]
  barplot(v_p, axes=FALSE, border="black")
  title(main = titles[type], font.main=1, cex.main=1.3, line=1.5)
}
# Ordinal
for(type in 1:2) {
  barplot(prob_table[type, ], axes=FALSE)
  title(main = titles[type], font.main=1, cex.main=1.3, line=1.5)
}

## Data
par(mar=c(4,4,0,1))
for(dom in 1:2) {
  for(type in 1:2) {

    plot.new()
    plot.window(xlim=c(1,10), ylim=c(0,1))
    axis(1, labels = n_seq, at=1:10, las=2)
    axis(2, las=2)
    title(ylab="Proportion Correct")
    title(xlab="Sample size", line=3)

    lines(a_out_skew_agg[type, , dom], lwd=2)

  }
}

dev.off()




